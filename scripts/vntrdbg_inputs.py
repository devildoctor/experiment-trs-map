"""Read existing anchored evidence plus conservative alignment-local fragments.

Only exact slices of one stored BAM record are exported. Unanchored alignment
fragments are explicitly local-only evidence, never complete allele paths.
"""
from collections import Counter,defaultdict
from hashlib import sha256
from pathlib import Path
import csv
import json
import re
import subprocess

import vntranchors
from vntrcore import compact,decompose
from vntrgraph import IndexedFastaReference,discover_bams,read_catalog
from vntrdbg import digest,save_json,tsv


def rows(path):
    with Path(path).open(encoding='utf-8',newline='') as handle:
        return list(csv.DictReader(handle,delimiter='\t'))


def intervals(record,start,end):
    """Query intervals bounded by actually aligned bases; N splits the path.

    Internal insertions remain real query bases; terminal CIGAR-only boundaries
    are uncertain. Deletions add no bases; supplementary records are not joined.
    """
    ref,q=record.pos,0
    left=right=None
    found=[]
    for n,op in record.ops:
        if op in 'M=X':
            lo,hi=max(ref,start),min(ref+n,end)
            if lo<hi:
                if left is None:
                    left=q+lo-ref
                right=q+hi-ref
            ref+=n
            q+=n
        elif op=='N':
            if left is not None:
                found.append((left,right))
            left=right=None
            ref+=n
        elif op=='D':
            ref+=n
        elif op in 'IS':
            q+=n
    if left is not None:
        found.append((left,right))
    return found


def prepare(args):
    csv.field_size_limit(64*1024*1024)
    models=rows(args.models)
    match=[m for m in models if m['locus']==args.locus]
    if len(match)!=1 or match[0]['model']!='compact_dbg':
        raise ValueError('Locus must explicitly select compact_dbg in the model catalog')
    loci=[x for x in read_catalog(args.catalog) if x.name==args.locus]
    if len(loci)!=1 or loci[0].str_structure:
        raise ValueError('Requires a single non-STR locus')
    locus=loci[0]
    params={r['parameter']:r['value'] for r in rows(args.evidence_dir/'parameters.tsv')}
    if params['catalog_sha256']!=sha256(args.catalog.read_bytes()).hexdigest():
        raise ValueError('Evidence catalog does not match requested catalog')
    if params['reference_fai_sha256']!=sha256(Path(str(args.reference)+'.fai').read_bytes()).hexdigest():
        raise ValueError('Evidence reference index does not match')
    if int(params['reference_size'])!=args.reference.stat().st_size or int(params['reference_mtime_ns'])!=args.reference.stat().st_mtime_ns:
        raise ValueError('Reference changed since evidence extraction')
    evidence_bams={str(Path(r['bam']).resolve()) for r in rows(args.evidence_dir/'bam_audit.tsv')}
    bams=discover_bams(args.bam_dir,args.locus)
    if evidence_bams!={str(b.resolve()) for _,_,b in bams}:
        raise ValueError('BAM inventory differs from the anchor extraction')
    decisions={r['molecule_key']:r for r in map(json.loads,(args.evidence_dir/'anchor_decisions.jsonl').read_text().splitlines())}
    complete={r['molecule_key']:r for r in rows(args.evidence_dir/'molecules.tsv')}
    records_by_mol=defaultdict(list)
    audit=[]
    namespaces={}
    source_table=params.get('sources','None')
    if source_table not in ('None',''):
        source_path=Path(source_table)
        for r in rows(source_path):
            p=Path(r['bam'])
            namespaces[str((p if p.is_absolute() else source_path.parent/p).resolve())]=r['namespace']
    with IndexedFastaReference(args.reference) as reference:
        refseq=reference.fetch(locus.chrom,locus.start,locus.end)
        for sample,platform,bam in bams:
            subprocess.run([args.samtools,'quickcheck',str(bam)],check=True)
            text=subprocess.check_output([args.samtools,'view','-h',str(bam)],text=True)
            rr,stats=vntranchors.parse_records(text.splitlines(),sample,reference,locus,str(bam),namespaces.get(str(bam.resolve())),platform)
            for r in rr:
                records_by_mol[compact([sample,r.namespace,r.name])].append(r)
            audit.append(dict(bam=str(bam),sample=sample,platform=platform,records=len(rr),stats=dict(stats),
                              bam_sha256=sha256(bam.read_bytes()).hexdigest()))
    output=[dict(id='REF_GRCh38_'+args.locus,kind='reference',sample='',sequence=refseq,
                 sequence_sha256=digest(refseq),coverage='catalog_interval',platforms=[],source_records=[])]
    provenance=[]
    for key,group in sorted(records_by_mol.items()):
        sample=json.loads(key)[0]
        decision=decisions.get(key,{}).get('status','not_in_anchor_audit')
        if not vntranchors.consistent_sequences(group):
            provenance.append(dict(molecule_key=key,status='excluded_inconsistent_sequences',anchor_status=decision,paths=[]))
            continue
        emitted=[]
        if key in complete:
            c=complete[key]
            seq=c['sequence']
            a,b=int(c['qL']),int(c['qR'])
            sources=json.loads(c['records'])
            if not any(r.source in sources and r.sequence[a:b]==seq for r in group):
                raise AssertionError('Anchored sequence cannot be traced to current BAM SEQ')
            identifier='OBS_'+sample+'_'+sha256(key.encode()).hexdigest()[:14]+'_A'
            item=dict(id=identifier,kind='observed',sample=sample,molecule_key=key,sequence=seq,
                      sequence_sha256=digest(seq),coverage='two_flank_anchored',
                      platforms=json.loads(c['platforms']) or [c['platform']],source_records=sources,
                      stored_query_interval=[a,b],original_query_interval=json.loads(c['original_query_interval']),
                      anchor_status=decision,boundary_status=locus.boundary_status)
            output.append(item)
            emitted.append(identifier)
        elif decision in {'unresolved_molecule_conflict','conflicting_complete_observations','inconsistent_record_sequences_or_lengths'}:
            provenance.append(dict(molecule_key=key,status='excluded_anchor_conflict',anchor_status=decision,paths=[]))
            continue
        else:
            fragments={}
            for r in group:
                if r.sequence=='*' or r.chrom!=locus.chrom or r.flag&(4|256|512|1024) or r.mapq<20 or r.mapq==255:
                    continue
                for a,b in intervals(r,locus.start,locus.end):
                    for part in re.finditer('[ACGT]+',r.sequence[a:b]):
                        lo,hi=a+part.start(),a+part.end()
                        if hi-lo<2:
                            continue
                        seq=r.sequence[lo:hi]
                        original=r.original_interval(lo,hi)
                        identity=(original,seq)
                        item=fragments.setdefault(identity,dict(sequence=seq,platforms=set(),source_records=set(),stored_intervals=[]))
                        item['platforms'].add(r.platform)
                        item['source_records'].add(r.source)
                        item['stored_intervals'].append(dict(record=r.source,start=lo,end=hi))
            for (original,seq),item in sorted(fragments.items()):
                name='OBS_'+sample+'_'+sha256((key+str(original)+seq).encode()).hexdigest()[:14]+'_L'
                output.append(dict(item,id=name,kind='observed',sample=sample,molecule_key=key,
                                   sequence_sha256=digest(seq),platforms=sorted(item['platforms']),
                                   source_records=sorted(item['source_records']),original_query_interval=original,
                                   coverage='alignment_local_only',anchor_status=decision,
                                   boundary_status='CIGAR_projected_uncertain'))
                emitted.append(name)
        provenance.append(dict(molecule_key=key,status='included' if emitted else 'no_eligible_local_sequence',
                               anchor_status=decision,paths=emitted))
    annotations={}
    for c in rows(args.evidence_dir/'audit.tsv'):
        name='CAND_'+c['allele']
        output.append(dict(id=name,kind='candidate',sample=c['sample'],sequence=c['sequence'],
                           sequence_sha256=digest(c['sequence']),coverage='candidate_sequence',
                           platforms=[],source_records=[],uncertain=c['uncertain'],flags=json.loads(c['flags']),
                           candidate_molecule_ids=json.loads(c['molecule_ids'])))
    for t in rows(args.evidence_dir/'tokens.tsv'):
        item=dict(t,start=int(t['start']),end=int(t['end']))
        annotations.setdefault('CAND_'+t['allele'],[]).append(item)
    motifs=[r['motif'] for r in rows(args.evidence_dir/'motifs.tsv') if r['usable']=='True' and r['level']=='base']
    if motifs:
        tokens,_,_=decompose(refseq,tuple(motifs),return_details=True)
        annotations[output[0]['id']]=[dict(start=t.start,end=t.end,exact_bases=t.exact_bases,family=t.family,status=t.status) for t in tokens]
    if len({r['id'] for r in output})!=len(output):
        raise AssertionError('Path ID collision')
    for name,tt in annotations.items():
        seq=next(r['sequence'] for r in output if r['id']==name)
        assert ''.join(t['exact_bases'] for t in sorted(tt,key=lambda t:t['start']))==seq
    save_json(args.out/'input_paths.json',output)
    save_json(args.out/'motif_annotations.json',annotations)
    tsv(args.out/'input_audit.tsv',audit)
    tsv(args.out/'molecule_audit.tsv',provenance)
    tsv(args.out/'classification.tsv',models)
    summary=dict(locus=args.locus,bams=len(bams),input_records=sum(r['records'] for r in audit),
                 path_types=dict(Counter(r['kind'] for r in output)),
                 observation_scope=dict(Counter(r['coverage'] for r in output if r['kind']=='observed')),
                 independent_molecules=len({r['molecule_key'] for r in output if r['kind']=='observed'}),
                 molecule_decisions=dict(Counter(r['status'] for r in provenance)),
                 evidence_dir=str(args.evidence_dir),candidate_support_is_not_observed_edge_support=True,
                 alignment_fragments='CIGAR-local-only, MAPQ>=20; not flank-unique or complete-allele evidence',
                 source_evidence_hashes={p.name:sha256(p.read_bytes()).hexdigest() for p in args.evidence_dir.iterdir() if p.is_file() and p.suffix in {'.tsv','.json','.jsonl'}},
                 catalog_sha256=sha256(args.catalog.read_bytes()).hexdigest())
    save_json(args.out/'input_summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k!='source_evidence_hashes'}),flush=True)
    return output,annotations
