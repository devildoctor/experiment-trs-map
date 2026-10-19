#!/usr/bin/env python3
"""Evidence-preserving motif structure graphs (not a genotype/SV caller).

Catalog seeds guide annotation; repeated HiFi observations can add candidate
variants. Exact sequence graphs and abstract motif graphs are separate products.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
from hashlib import sha256
from html import escape
from pathlib import Path
import json
import math
import re

from vntrgraph import IndexedFastaReference, read_catalog, discover_bams
from vntrreads import read_bam, write_table, collect_reads, annotation_results
from vntrcore import edit_distance, revcomp, compact, canonical_motif
from vntrsummary import structural_summary, reference_report, exact_runs, summary_svg

from vntrmotifs import primitive_unit, definition, forms, context_support, TwoPaths, base_metrics, compound_annotations

VERSION='5.1-phase1'
PALETTE=['#3874ad','#bf722b','#927c29','#a34e83','#718441']


def support_key(row):
    # BAM copies and multiple alignment records must not multiply support.
    return compact([row['sample'],row['read_id']])


def seed_dictionary(motifs):
    if not motifs:
        raise ValueError('Structure mode requires MOTIFS in BED; seeds are not inferred silently')
    unique=[]
    for motif in motifs:
        if not motif or set(motif)-set('ACGT'):
            raise ValueError('MOTIFS must contain nonempty ACGT strings')
        if motif not in unique: unique.append(motif)
    return normalize_dictionary([dict(id=f'M{i+1}',base=f'M{i+1}',sequence=m,source='catalog_seed',support=0)
            for i,m in enumerate(unique)])


def normalize_dictionary(dictionary):
    for record in list(dictionary):
        record.update(definition(record['sequence']))
        record.setdefault('evidence_status','not_observed_in_repeat_context')
        primitive=record['primitive_sequence']
        if primitive != record['sequence'] and not any(canonical_motif(m['sequence'])==canonical_motif(primitive) for m in dictionary):
            name=f'M{len(dictionary)+1}'
            dictionary.append(dict(id=name,base=name,sequence=primitive,source='derived_primitive',
                                   parent_id=record['id'],support=0,**definition(primitive),evidence_status='not_observed_in_repeat_context'))
    for record in dictionary:
        unit=record['sequence']
        shorter={v for m in dictionary if len(m['sequence'])<len(unit) for v in forms(m['sequence'])}
        reachable={0}
        for i in range(len(unit)):
            if i in reachable:
                reachable.update(i+len(m) for m in shorter if unit.startswith(m,i))
        if len(unit) in reachable:record['level']='compound'
    return dictionary


def add_reference_units(dictionary, sequence):
    """Use the best seed phase to enumerate reference units, including interruptions."""
    keys={canonical_motif(m['sequence']) for m in dictionary}
    for seed in list(dictionary):
        length=len(seed['sequence'])
        if length>len(sequence):continue
        options=[]
        for offset in range(length):
            units=[sequence[i:i+length] for i in range(offset,len(sequence)-length+1,length)]
            if not units:continue
            similar=sum(edit_distance(u,seed['sequence'])<=max(1,length//4) for u in units)
            exact=sum(u==seed['sequence'] for u in units)
            options.append((similar,exact,-offset,units))
        if not options:continue
        units=max(options,key=lambda x:x[:3])[3]
        for unit in dict.fromkeys(units):
            if set(unit)-set('ACGT') or edit_distance(unit,seed['sequence'])>max(1,length//4):continue
            key=canonical_motif(unit)
            if key in keys:continue
            keys.add(key)
            name=f'M{len(dictionary)+1}'
            dictionary.append(dict(id=name,base=name,parent_id=seed['id'],sequence=unit,
                                   source='reference_unit',support=0))
    # Scan independently: approximate seed matches must not hide a new period.
    exact_covered=set()
    for token in segment((sequence,dictionary,'HiFi')):
        if token['status']=='unit' and token['edit_distance']==0:exact_covered.update(range(token['start'],token['end']))
    for hit in exact_runs(sequence):
        unit=hit['sequence'];key=canonical_motif(unit)
        if key in keys or sum(i in exact_covered for i in range(hit['start'],hit['end']))>=.8*(hit['end']-hit['start']):continue
        name=f'M{len(dictionary)+1}';keys.add(key)
        dictionary.append(dict(id=name,base=name,sequence=unit,source='reference_period',support=0))
    return normalize_dictionary(dictionary)


def exact_tandem_units(sequence,min_length=2,max_length=200):
    return {h['sequence'] for h in exact_runs(sequence,max_period=max_length,
            min_period=min_length,min_copies=2,min_span=2*min_length)}


def segment(task,return_details=False):
    """Seed-guided exact-slice DP. Gap/length penalties exceed UNKNOWN cost.

    Unknown=1/base; substitutions cost min(3, 0.6*motif length). Both motif
    orientations are tested, not confused with the BAM alignment strand.
    """
    sequence,dictionary,platform=task
    n=len(sequence)
    variants=[]
    seen=set()
    for m in dictionary:
        if m.get('level',definition(m['sequence'])['level']) != 'base':continue
        for orientation,s in (('+',m['sequence']),('-',revcomp(m['sequence']))):
            if (m['id'],s) not in seen:
                variants.append((m,s,orientation));seen.add((m['id'],s))
    candidates=defaultdict(set)
    for v,(m,s,orientation) in enumerate(variants):
        length=len(s); k=min(length,5,max(2,length//5))
        slip=max(1,min(8,math.ceil(length*.08)))
        for offset in sorted({0,max(0,length//3),max(0,2*length//3-k),length-k}):
            hit=sequence.find(s[offset:offset+k])
            while hit>=0:
                for delta in range(-slip,slip+1):
                    start=hit-offset+delta
                    if 0<=start<n:candidates[start].add(v)
                hit=sequence.find(s[offset:offset+k],hit+1)
    dp=TwoPaths(n);cost=dp.cost
    def update(end,value,step):
        dp.offer(step[0],end,value-cost[step[0]],(step[0],end)+step[1:])
    for i in range(n):
        update(i+1,cost[i]+1,(i,-1,'unknown',0))
        for v in sorted(candidates.get(i,())):
            m,s,orientation=variants[v];length=len(s)
            allowed=0 if length==1 else 1 if length<=6 else max(1,math.floor(length*(.32 if platform=='HiFi' else .38)))
            slip=0 if length==1 else max(1,min(8,math.ceil(length*(.08 if platform=='HiFi' else .16))))
            for width in range(max(1 if length==1 else 2,length-slip),min(n-i,length+slip)+1):
                d=edit_distance(s,sequence[i:i+width])
                if d<=allowed:
                    update(i+width,cost[i]+.1+min(3,length*.6)*d+2*abs(width-length),(i,v,'unit',d))
        # Very short terminal matches remain UNKNOWN instead of 1bp "copies".
        for v,(m,s,orientation) in enumerate(variants):
            minimum=max(2,min(8,len(s)//4))
            if i==0:
                for width in range(minimum,min(len(s),n+1)):
                    if sequence[:width]==s[-width:]:update(width,.2,(0,v,'partial',0))
            if minimum<=n-i<len(s) and sequence[i:]==s[:n-i]:
                update(n,cost[i]+.2,(i,v,'partial',0))
    audit=dp.audit();alternate=dp.trace(1);alternate_set=set(alternate)
    def convert(steps):
        result=[]
        for start,end,v,status,d in steps:
            m,s,orientation=variants[v] if v>=0 else ({'id':'U','base':'U'},'','?')
            family=canonical_motif(s) if s else 'UNKNOWN'
            from vntrcore import motif_orientation
            phase,canonical_orientation=motif_orientation(s,family) if s else (0,'?')
            t=dict(start=start,end=end,exact_bases=sequence[start:end],status=status,
                   motif=m['id'],base=m['base'],template=s,orientation=orientation,family=family,
                   phase=phase,phase_policy='fixed_template_phase',canonical_orientation=canonical_orientation,
                   orientation_reference='catalog_template',phase_reference='canonical_family',edit_distance=d,
                   annotation_ambiguity=audit['annotation_ambiguity'] and (start,end,v,status,d) not in alternate_set)
            if status=='unknown' and result and result[-1]['status']=='unknown':
                result[-1]['end']=end;result[-1]['exact_bases']+=t['exact_bases']
                result[-1]['annotation_ambiguity']|=t['annotation_ambiguity']
            else:result.append(t)
        assert ''.join(t['exact_bases'] for t in result)==sequence
        return result
    result=convert(dp.trace())
    audit['metrics']=base_metrics(result)
    audit['alternative']=convert(alternate) if alternate else []
    audit['alternative_metrics']=base_metrics(audit['alternative']) if alternate else None
    return (result,audit) if return_details else result


def segment_details(task):
    return segment(task,return_details=True)


def annotate_rows(rows,dictionary,workers):
    fragments=[(r,f) for r in rows for f in r['fragments']]
    tasks=((f['sequence'],dictionary,r['platform']) for r,f in fragments)
    for i,((r,f),(tokens,annotation)) in enumerate(zip(fragments,annotation_results(segment_details,tasks,workers)),1):
        f['tokens']=tokens;f['annotation']=annotation
        if workers>1 and i%100==0:print(f'Annotated {i}/{len(fragments)}',flush=True)


def learn_variants(rows,dictionary,minimum=3,limit=256,min_fraction=.10):
    """Discover primitive local repeats independently of current annotation.

    Exact recurring interruptions are also retained, with raw-sequence support.
    """
    existing={canonical_motif(m['sequence']) for m in dictionary}
    observed=defaultdict(set);by_sample=defaultdict(lambda:defaultdict(set))
    sources=defaultdict(set);informative=defaultdict(set)
    for r in rows:
        if r['platform']!='HiFi':continue
        for f in r['fragments']:
            if f['alignment_type']!='primary' or f['mapq']<20 or set(f.get('flags',[]))&{'duplicate_flag','qc_fail'}:continue
            informative[r['sample']].add(support_key(r))
            sequence=f.get('sequence',''.join(t['exact_bases'] for t in f['tokens']))
            runs=exact_runs(sequence,max_period=200,min_period=2,min_copies=2,min_span=12)
            candidates=[(h['sequence'],'','local_exact_tandem') for h in runs]
            for i,t in enumerate(f['tokens']):
                seq=t['exact_bases'];left=f['tokens'][i-1] if i else None;right=f['tokens'][i+1] if i+1<len(f['tokens']) else None
                neighbors=[x for x in (left,right) if x and x['status']=='unit']
                if any(h['start']<=t.get('start',-1) and t.get('end',len(sequence)+1)<=h['end'] for h in runs):continue
                if t['status']=='unit' and t['orientation']=='+':
                    if neighbors and abs(len(seq)-len(t['template']))<=max(1,math.ceil(len(t['template'])*.04)):
                        candidates.append((seq,t['base'],'matched_variant'))
                if t['status']=='unknown':
                    if 2<=len(seq)<=200 and len(neighbors)==2:
                        candidates.append((seq,left['base'],'unknown_interruption'))
            for candidate,base,source in candidates:
                candidate=primitive_unit(candidate)
                if not 2<=len(candidate)<=200 or set(candidate)-set('ACGT'):continue
                key=canonical_motif(candidate)
                if key in existing:continue
                # A candidate must have exact neighboring unit context in this molecule.
                context=context_support(sequence,[candidate]+[m['sequence'] for m in dictionary if m.get('level','base')=='base'])
                if not context[candidate]['context_occurrences']:continue
                item=(key,candidate,base);molecular=support_key(r)
                observed[item].add(molecular);by_sample[item][r['sample']].add(molecular);sources[item].add(source)
    grouped=defaultdict(list)
    for item,ids in observed.items():grouped[item[0]].append((item,ids))
    audit=[];admitted=[]
    for key,variants in grouped.items():
        item,ids=max(variants,key=lambda x:(len(x[1]),'local_exact_tandem' in sources[x[0]],x[0][1]))
        all_ids=set().union(*(v for _,v in variants));sample_counts=defaultdict(set)
        for variant,_ in variants:
            for sample,seen in by_sample[variant].items():sample_counts[sample].update(seen)
        fraction=max((len(seen)/len(informative[sample]) for sample,seen in sample_counts.items() if informative[sample]),default=0)
        source_set=set().union(*(sources[v] for v,_ in variants))
        reason='accepted' if len(all_ids)>=minimum and fraction>=min_fraction else 'low_HiFi_support' if len(all_ids)<minimum else 'low_sample_fraction'
        audit.append(dict(sequence=item[1],base=item[2],length=len(item[1]),support=len(all_ids),hifi_support=len(all_ids),sample_fraction=round(fraction,4),sources=','.join(sorted(source_set)),status=reason,phase_key=key))
        if reason=='accepted':admitted.append((item,len(all_ids),sorted(all_ids),source_set))
    admitted.sort(key=lambda x:(-x[1],x[0][1]))
    if len(dictionary)+len(admitted)>limit:raise ValueError(f'Motif candidate limit {limit} exceeded; no silent truncation')
    added=[]
    for item,support,ids,source_set in admitted:
        name=f'M{len(dictionary)+len(added)+1}'
        added.append(dict(id=name,base=name,parent_id=item[2],sequence=item[1],source='learned_'+','.join(sorted(source_set)),support=support,support_ids=ids,evidence_status='supported'))
    return normalize_dictionary(dictionary+added),sorted(audit,key=lambda x:(x['status']!='accepted',-x['support'],x['sequence']))


def run_blocks(tokens,family=False):
    blocks=[]
    for i,t in enumerate(tokens):
        key=(t['base'] if family else t['motif'])+('r' if t['orientation']=='-' else '')
        if t['status']=='unknown':key='U'
        if t['status']=='partial':key+='p'
        if blocks and blocks[-1]['label']==key and t['status']=='unit':
            blocks[-1]['count']+=1;blocks[-1]['end']=i+1;blocks[-1]['bp']+=len(t['exact_bases'])
        else:blocks.append(dict(label=key,count=1,start=i,end=i+1,bp=len(t['exact_bases']),status=t['status']))
    return blocks


def compound_blocks(tokens,catalog=()):
    return [dict(h,start=h['start_token'],end=h['end_token'],bp_start=h['start'],bp_end=h['end'],width=h['base_unit_count'])
            for h in compound_annotations(tokens,catalog)]


def summarize(rows,dictionary):
    exact={};abstract={};edges={};exact_edges={};patterns={}
    family_nodes={};family_edges={};base_lookup={m['id']:m['base'] for m in dictionary}
    def node(store,key,seq,status):
        if key not in store:store[key]=dict(id=key,sequence=seq,status=status,reads=set(),primary_reads=set(),occurrences=0,samples=defaultdict(set))
        return store[key]
    def add(n,row,fragment):
        key=support_key(row);n['reads'].add(key);n['occurrences']+=1;n['samples'][row['sample']].add(key)
        if fragment['alignment_type']=='primary' and fragment['mapq']>=20:n['primary_reads'].add(key)
    def edge(store,a,b,row,fragment):
        key=(a,b)
        if key not in store:store[key]=dict(source=a,target=b,reads=set(),primary_reads=set(),occurrences=0,samples=defaultdict(set))
        add(store[key],row,fragment)
    for row in rows:
        for f in row['fragments']:
            for t in f['tokens']:
                t['qstart']=f['qstart']+t['start'];t['qend']=f['qstart']+t['end']
                t['exact_node']='E'+sha256(t['exact_bases'].encode()).hexdigest()[:16]
                if t['status']=='unknown':sid='U'+sha256(t['exact_bases'].encode()).hexdigest()[:12]
                elif t['status']=='partial':sid=t['motif']+('r' if t['orientation']=='-' else '')+'p'+sha256(t['exact_bases'].encode()).hexdigest()[:8]
                else:sid=t['motif']+('r' if t['orientation']=='-' else '')
                t['structure_node']=sid
                add(node(exact,t['exact_node'],t['exact_bases'],'exact'),row,f)
                add(node(abstract,sid,t['template'] if t['status']=='unit' else t['exact_bases'],t['status']),row,f)
                family=t['base']+('r' if t['orientation']=='-' else '') if t['status']=='unit' else 'U' if t['status']=='unknown' else 'P'
                t['family_node']=family
                add(node(family_nodes,family,'','motif_family' if t['status']=='unit' else 'visual_bucket'),row,f)
            for a,b in zip(f['tokens'],f['tokens'][1:]):
                edge(edges,a['structure_node'],b['structure_node'],row,f)
                edge(exact_edges,a['exact_node'],b['exact_node'],row,f)
                edge(family_edges,a['family_node'],b['family_node'],row,f)
            f['blocks']=run_blocks(f['tokens']);f['family_blocks']=run_blocks(f['tokens'],True)
            f['compounds']=compound_blocks(f['tokens'],[m['sequence'] for m in dictionary if m['level']=='compound' and m['source']=='catalog_seed'])
            f['structure']=' → '.join(f"{b['label']}×{b['count']}" if b['status']=='unit' else f"{b['label']}({b['bp']}bp)" for b in f['blocks'])
            f['family_structure']=' → '.join(f"{b['label']}×{b['count']}" if b['status']=='unit' else f"{b['label']}({b['bp']}bp)" for b in f['family_blocks'])
            f['features']=[]
            if any(b['count']>=2 for b in f['blocks']):f['features'].append('tandem_run')
            if f['compounds']:f['features'].append('compound_candidate')
            if any(t['status']=='unknown' for t in f['tokens'][1:-1]):f['features'].append('unresolved_internal_segment')
            if any(t['orientation']=='-' for t in f['tokens']):f['features'].append('reverse_motif_candidate')
            signature=compact([f['coverage'],[t['structure_node'] for t in f['tokens']]])
            pid='T'+sha256(signature.encode()).hexdigest()[:12]
            f['pattern_id']=pid
            if pid not in patterns:patterns[pid]=dict(id=pid,structure=f['structure'],coverage=f['coverage'],reads=set(),primary_reads=set(),occurrences=0,samples=defaultdict(set))
            add(patterns[pid],row,f)
            assert ''.join(exact[t['exact_node']]['sequence'] for t in f['tokens'])==f['sequence']
    def finish(store):
        result=[]
        for key,r in sorted(store.items()):
            r['support']=len(r.pop('reads'));r['primary_support']=len(r.pop('primary_reads'))
            r['samples']={s:len(v) for s,v in r['samples'].items()};result.append(r)
        return result
    return dict(exact_nodes=finish(exact),structure_nodes=finish(abstract),edges=finish(edges),exact_edges=finish(exact_edges),
                patterns=sorted(finish(patterns),key=lambda p:(-p['support'],p['id'])),family_nodes=finish(family_nodes),family_edges=finish(family_edges))


def color(label):
    if label.startswith('U'):return '#87909b'
    if label.startswith('P'):return '#b8a87a'
    match=re.match(r'M(\d+)',label)
    base=PALETTE[(int(match[1])-1)%len(PALETTE)] if match else '#87909b'
    variant=re.search(r'v(\d+)',label)
    factor=1 if not variant else .58+.1*(int(variant[1])%4)
    return '#'+''.join(f'{round(int(base[i:i+2],16)*factor):02x}' for i in (1,3,5))


def graph_svg(graph,locus):
    """Small family-level topology; U/P are explicit display buckets only."""
    nodes=graph['family_nodes'];edges=graph['family_edges'];pos={n['id']:(130+i*245,155) for i,n in enumerate(nodes)}
    width=max(900,245*len(nodes)+40);height=360+len(edges)*22
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0 0 L8 4 L0 8" fill="#55677c"/></marker></defs>',
           '<rect width="100%" height="100%" fill="white"/>',f'<text x="24" y="28" font-size="20">{escape(locus)} motif-family graph / distinct-read support</text>',
           '<text x="24" y="51" font-size="13">U/P = visual buckets, not merged biological sequences. Only recorded read paths are supported.</text>']
    for e in edges:
        x,y=pos[e['source']];xx,yy=pos[e['target']]
        if x==xx:d=f'M {x-25} {y-26} C {x-100} {y-110}, {x+100} {y-110}, {x+25} {y-26}'
        elif x<xx:d=f'M {x+72} {y-10} Q {(x+xx)/2} {55-abs(xx-x)*.02} {xx-72} {yy-10}'
        else:d=f'M {x-72} {y+10} Q {(x+xx)/2} {250+abs(xx-x)*.02} {xx+72} {yy+10}'
        parts.append(f'<path d="{d}" fill="none" stroke="#55677c" stroke-width="1.5" marker-end="url(#arrow)"><title>{e["source"]} → {e["target"]}: {e["support"]} reads; {e["occurrences"]} traversals</title></path>')
    for n in nodes:
        x,y=pos[n['id']]
        parts.append(f'<rect x="{x-72}" y="{y-27}" width="144" height="54" fill="{color(n["id"])}" stroke="#26364a"/>')
        parts.append(f'<text x="{x}" y="{y-6}" text-anchor="middle" fill="white" font-size="17">{n["id"]}</text><text x="{x}" y="{y+15}" text-anchor="middle" fill="white" font-size="13">{n["support"]} reads</text>')
    for i,e in enumerate(sorted(edges,key=lambda e:(-e['support'],e['source'],e['target']))):
        parts.append(f'<text x="24" y="{320+i*22}" font-size="14">{e["source"]} → {e["target"]}: {e["support"]} reads / {e["primary_support"]} primary MAPQ≥20 / {e["occurrences"]} traversals</text>')
    return ''.join(parts)+'</svg>'


def rows_svg(rows,max_length,title):
    count=sum(max(1,len(r['fragments'])) for r in rows);height=100+count*65
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="1500" height="{height}" viewBox="0 0 1500 {height}" font-family="sans-serif"><rect width="100%" height="100%" fill="white"/>',f'<text x="20" y="25" font-size="18">{escape(title)}</text>',
           f'<text x="20" y="49" font-size="12">Observed read-relative bp, not reference alignment. Shared scale 0–{max_length} bp. All rows retained.</text>']
    scale=1090/max(1,max_length);y=83
    for row in rows:
        fragments=row['fragments'] or [None]
        for f in fragments:
            label=f'{row["sample"]}/{row["platform"]} {row["read_id"]}'
            parts.append(f'<text x="20" y="{y}" font-size="11"><title>{escape(label)}</title>{escape(label[:49])}</text>')
            if not f:
                parts.append(f'<text x="390" y="{y}" font-size="12">{escape(row["status"])}</text>');y+=65;continue
            parts.append(f'<text x="20" y="{y+18}" font-size="10">{escape(f["coverage"])} / {escape(f["anchor"])}</text>')
            for t in f['tokens']:
                x=390+t['start']*scale;w=(t['end']-t['start'])*scale
                label=t['motif']+('r' if t['orientation']=='-' else '')
                dashed=' stroke-dasharray="3 2"' if t['status']=='partial' else ''
                parts.append(f'<rect x="{x:.2f}" y="{y-12}" width="{w:.2f}" height="22" fill="{color(label)}" stroke="#fff" stroke-width=".5"{dashed}><title>{escape(label)}; {t["structure_node"]}; [{t["qstart"]},{t["qend"]}); edit={t["edit_distance"]}; {escape(t["exact_bases"])}</title></rect>')
                if w>=25:parts.append(f'<text x="{x+w/2:.2f}" y="{y+3}" text-anchor="middle" font-size="9" fill="white">{escape(label)}</text>')
            summary=f['family_structure']
            parts.append(f'<text x="390" y="{y+28}" font-size="11"><title>{escape(summary)}</title>{escape(summary[:140])}{"… (full path in HTML)" if len(summary)>140 else ""}</text>')
            y+=65
    return ''.join(parts)+'</svg>'


def write_gfa(path,nodes,edges,rows,abstract=False):
    with path.open('w',encoding='utf-8') as out:
        out.write('H\tVN:Z:1.0\tPG:Z:vntrstructure-4\tDS:Z:'+('abstract_motif_graph_not_base_reconstruction' if abstract else 'exact_sequence_graph')+'\n')
        for n in nodes:
            seq='*' if abstract else n['sequence']
            out.write(f'S\t{n["id"]}\t{seq}\tLN:i:{len(n["sequence"])}\tRC:i:{n["support"]}\tPC:i:{n["primary_support"]}\tOC:i:{n["occurrences"]}\tTS:Z:{n["sequence"]}\n')
        for e in edges:out.write(f'L\t{e["source"]}\t+\t{e["target"]}\t+\t0M\tRC:i:{e["support"]}\tOC:i:{e["occurrences"]}\n')
        key='structure_node' if abstract else 'exact_node'
        for r in rows:
            for f in r['fragments']:
                p=','.join(t[key]+'+' for t in f['tokens'])
                if p:out.write(f'P\t{f["id"]}\t{p}\t*\tRN:Z:{r["read_id"]}\tCV:Z:{f["coverage"]}\n')


def emit(destination,locus,rows,dictionary,audit,reference_structure=None,summary_minimum=3,summary_fraction=.10):
    graph=summarize(rows,dictionary)
    summary=structural_summary(rows,summary_minimum,summary_fraction)
    audit['structure_summary']=summary['audit']
    tokens=[t for r in rows for f in r['fragments'] for t in f['tokens']]
    fragments=[f for r in rows for f in r['fragments']]
    motif_reads=defaultdict(set)
    for r in rows:
        for f in r['fragments']:
            for t in f['tokens']:
                if t['status']=='unit':motif_reads[t['motif']].add(support_key(r))
    evidence=defaultdict(set);positions=defaultdict(list)
    for r in rows:
        if r['platform']!='HiFi':continue
        for f in r['fragments']:
            if f['alignment_type']!='primary' or f['mapq']<20 or set(f.get('flags',[]))&{'duplicate_flag','qc_fail'}:continue
            hits=context_support(f['sequence'],[m['sequence'] for m in dictionary])
            for motif in dictionary:
                for start,end in hits[motif['sequence']]['context_occurrences']:
                    evidence[motif['id']].add(support_key(r))
                    positions[motif['id']].append(dict(molecule_key=support_key(r),fragment_id=f['id'],start=start,end=end))
    for motif in dictionary:
        motif['final_read_support']=len(motif_reads[motif['id']])
        motif['evidence_read_support']=len(evidence[motif['id']])
        motif['support_positions']=positions[motif['id']]
        motif['evidence_status']='supported' if len(evidence[motif['id']])>=summary_minimum else 'observed_below_threshold' if evidence[motif['id']] else 'not_observed_in_repeat_context'
    for f in fragments:
        f['base_metrics']=base_metrics(f['tokens'])
    total=sum(len(t['exact_bases']) for t in tokens);unknown=sum(len(t['exact_bases']) for t in tokens if t['status']=='unknown')
    audit.update(version=VERSION,status='complete',graph_fragments=len(fragments),token_occurrences=len(tokens),
                 exact_nodes=len(graph['exact_nodes']),structure_nodes=len(graph['structure_nodes']),
                 shared_structure_nodes=sum(n['support']>1 for n in graph['structure_nodes']),
                 self_loop_edges=sum(e['source']==e['target'] for e in graph['edges']),
                 unknown_bases=unknown,total_observed_bases=total,unknown_fraction=unknown/max(1,total),
                 sequence_reconstruction_checks_passed=len(fragments),
                 feature_counts=dict(Counter(x for f in fragments for x in f['features'])),
                 support_definition='distinct sample/QNAME; repeated units, BAM copies, platform labels and alignments counted once per node/edge',
                 limitations='Exploratory seed-guided annotation. Candidate-window boundaries. Abstract graph is not allele sequence; U/P overview buckets do not assert sequence identity; merged graph combinations require observed read paths.')
    audit['platform_annotation']={}
    for platform in sorted({r['platform'] for r in rows}):
        ts=[t for r in rows if r['platform']==platform for f in r['fragments'] for t in f['tokens']]
        bases=sum(len(t['exact_bases']) for t in ts)
        audit['platform_annotation'][platform]=dict(bases=bases,unknown_fraction=sum(len(t['exact_bases']) for t in ts if t['status']=='unknown')/max(1,bases))
    write_gfa(destination/f'{locus.name}.reads.gfa',graph['exact_nodes'],graph['exact_edges'],rows)
    write_gfa(destination/f'{locus.name}.structure.gfa',graph['structure_nodes'],graph['edges'],rows,True)
    write_table(destination/'motifs.tsv',dictionary,['id','base','parent_id','sequence','source','level','primitive_sequence','repeat_factor','family','phase_policy','support','final_read_support','evidence_read_support','evidence_status','support_ids','support_positions'])
    write_table(destination/'motif_candidates.tsv',audit.get('motif_candidates',[]),['round','sequence','base','length','support','hifi_support','sample_fraction','sources','status','phase_key'])
    write_table(destination/'structure_classes.tsv',summary['classes'],['id','skeleton','status','support','hifi_support','samples','length_min','length_max','subtypes','representative'])
    write_table(destination/'local_modules.tsv',summary['modules'],['id','kind','structure','status','support','hifi_support','exact_support','approximate_support','samples','occurrences','instances'])
    write_table(destination/'read_assignments.tsv',summary['assignments'],['row_id','read_id','fragment_id','sample','platform','class_id','subtype_id','evidence','reason','length','module_ids'])
    (destination/'reference_structure.json').write_text(json.dumps(reference_structure,ensure_ascii=False,indent=2),encoding='utf-8')
    write_table(destination/'structure_nodes.tsv',graph['structure_nodes'],['id','sequence','status','support','primary_support','occurrences','samples'])
    write_table(destination/'structure_edges.tsv',graph['edges'],['source','target','support','primary_support','occurrences','samples'])
    write_table(destination/'patterns.tsv',graph['patterns'],['id','structure','coverage','support','primary_support','occurrences','samples'])
    write_table(destination/'reads.tsv',[dict(r,fragment_count=len(r['fragments'])) for r in rows],['id','read_id','sample','platform','bam','record_count','fragment_count','status'])
    write_table(destination/'read_nodes.tsv',[dict(t,read_id=r['read_id'],row_id=r['id'],fragment_id=f['id'],order=i+1) for r in rows for f in r['fragments'] for i,t in enumerate(f['tokens'])],['row_id','read_id','fragment_id','order','structure_node','exact_node','motif','base','orientation','phase','phase_policy','status','qstart','qend','edit_distance','annotation_ambiguity','exact_bases'])
    write_table(destination/'compound_blocks.tsv',[dict(h,fragment_id=f['id']) for f in fragments for h in f['compounds']],['fragment_id','start_token','end_token','bp_start','bp_end','pattern','count','units_per_block','base_unit_count','catalog_motif','basis','exact_bases'])
    (destination/'annotation_alternatives.json').write_text(json.dumps({f['id']:f.get('annotation',{}) for f in fragments},ensure_ascii=False,indent=2),encoding='utf-8')
    graph_image=summary_svg(summary,locus.name)
    (destination/f'{locus.name}.family-transitions.svg').write_text(graph_svg(graph,locus.name),encoding='utf-8')
    (destination/f'{locus.name}.overview.svg').write_text(graph_image,encoding='utf-8')
    svgdir=destination/'svg';svgdir.mkdir()
    max_length=max([len(f['sequence']) for f in fragments]+[1]);links=[]
    # Stable sample/platform/length ordering, with no read filtering.
    ordered=sorted(rows,key=lambda r:(r['sample'],r['platform'],-max([len(f['sequence']) for f in r['fragments']]+[0]),r['read_id']))
    for start in range(0,len(ordered),40):
        number=start//40+1;name=f'{locus.name}.structures.{number:04d}.svg'
        (svgdir/name).write_text(rows_svg(ordered[start:start+40],max_length,f'{locus.name} all reads | page {number} | rows {start+1}–{min(start+40,len(rows))} / {len(rows)}'),encoding='utf-8')
        links.append(f'<li><a href="{name}">Page {number}</a></li>')
    (svgdir/'index.html').write_text('<!doctype html><meta charset="utf-8"><h2>全部读段，无抽样</h2><ol>'+''.join(links)+'</ol>',encoding='utf-8')
    audit['svg_pages']=len(links)
    payload=dict(locus=asdict(locus),rows=ordered,motifs=dictionary,graph=graph,audit=audit,max_length=max_length,
                 summary=summary,reference_structure=reference_structure)
    (destination/f'{locus.name}.structures.json').write_text(compact(payload),encoding='utf-8')
    (destination/'status.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    html=Path(__file__).with_name('structure_summary_view.html').read_text(encoding='utf-8')
    serialized=compact(payload).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
    (destination/f'{locus.name}.structure.html').write_text(html.replace('__TITLE__',locus.name).replace('__GRAPH__',graph_image).replace('__DATA__',serialized),encoding='utf-8')
    return audit


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bam-dir',type=Path,required=True);p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--catalog',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--locus',required=True);p.add_argument('--samtools',default='samtools');p.add_argument('--workers',type=int,default=2)
    p.add_argument('--variant-support',type=int,default=3);p.add_argument('--max-variants',type=int,default=256)
    p.add_argument('--min-sample-fraction',type=float,default=.10);p.add_argument('--discovery-rounds',type=int,default=3)
    p.add_argument('--gene-strand',choices=['+','-','unknown'],default='unknown')
    args=p.parse_args()
    if args.workers<1 or args.variant_support<2 or args.max_variants<1 or not 0<args.min_sample_fraction<=1 or not 1<=args.discovery_rounds<=3:p.error('Invalid workers/support/variant/discovery limits')
    selected=[l for l in read_catalog(args.catalog) if l.name==args.locus]
    if len(selected)!=1:p.error('Unknown locus')
    locus=selected[0];dictionary=seed_dictionary(locus.motifs)
    destination=args.out/locus.name;destination.mkdir(parents=True,exist_ok=False)
    rows=[];audits=[]
    with IndexedFastaReference(args.reference) as reference:
        reference_sequence=reference.fetch(locus.chrom,locus.start,locus.end).upper()
        dictionary=add_reference_units(dictionary,reference_sequence)
        ref_start=max(0,(locus.search_start if locus.search_start is not None else locus.start)-100)
        ref_end=min(reference.contig_length(locus.chrom),(locus.search_end if locus.search_end is not None else locus.end)+100)
        ref_window=reference.fetch(locus.chrom,ref_start,ref_end)
        if len(dictionary)>args.max_variants:raise ValueError('Reference motif dictionary exceeds --max-variants')
        with (destination/'records.jsonl').open('w',encoding='utf-8') as records:
            bams=discover_bams(args.bam_dir,locus.name)
            if not bams:raise ValueError('No BAM input')
            rows,audits=collect_reads(bams,locus,reference,args.samtools,100,records)
        print(f'Discovery: {len(rows)} rows; {len(dictionary)} BED/reference motifs',flush=True)
        candidate_audit=[]
        hifi=[r for r in rows if r['platform']=='HiFi']
        for round_number in range(1,args.discovery_rounds+1):
            annotate_rows(hifi,dictionary,args.workers)
            updated,candidates=learn_variants(hifi,dictionary,args.variant_support,args.max_variants,args.min_sample_fraction)
            candidate_audit.extend(dict(item,round=round_number) for item in candidates)
            added=len(updated)-len(dictionary)
            dictionary=updated
            print(f'Discovery round {round_number}: {added} motifs added; {len(dictionary)} total',flush=True)
            if not added:break
        annotate_rows(rows,dictionary,args.workers)
        reference_structure=reference_report(locus,reference_sequence,ref_window,ref_start,dictionary,segment,args.gene_strand,args.reference)
    audit=dict(locus=asdict(locus),bam_count=len(bams),read_rows=len(rows),rows_with_repeat_sequence=sum(bool(r['fragments']) for r in rows),
               rows_without_extractable_repeat=sum(not r['fragments'] for r in rows),alignment_records=sum(r['record_count'] for r in rows),bam_audits=audits,
               motif_candidates=candidate_audit,discovery_rounds_run=round_number,
               motif_source_counts=dict(Counter(m['source'] for m in dictionary)))
    audit=emit(destination,locus,rows,dictionary,audit,reference_structure,args.variant_support,args.min_sample_fraction)
    hashes={f.name:sha256(f.read_bytes()).hexdigest() for f in Path(__file__).parent.glob('vntr*.py')}
    hashes['structure_summary_view.html']=sha256(Path(__file__).with_name('structure_summary_view.html').read_bytes()).hexdigest()
    (destination/'parameters.json').write_text(json.dumps(dict(args={k:str(v) for k,v in vars(args).items()},source_hashes=hashes,catalog_sha256=sha256(args.catalog.read_bytes()).hexdigest()),ensure_ascii=False,indent=2),encoding='utf-8')
    print(compact({k:v for k,v in audit.items() if k!='bam_audits'}),flush=True)


if __name__=='__main__':main()
