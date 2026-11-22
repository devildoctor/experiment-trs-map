"""Regression cases from the November 22 point-graph review."""
import json
import os
import shutil
import subprocess
from hashlib import sha256
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from run_efficient_family import run, member_identity
from vntrflankalign import build_context, pair_endpoints, PREVIOUS_VERSION, write_flank_outputs
from vntrgraph import Locus
from vntrhmm import segment_hmm
from vntrreads import extract_record
from vntrpathstructure import annotate_paths
from vntrpartition import load_model
from vntrstructure import annotate_rows


class ReviewRepairs(unittest.TestCase):
    def test_reference_identity_evidence_and_missing_fields(self):
        model=dict(locus='L',assembly='GRCh38',chrom='chr1',start=10,end=13)
        ref=dict(id='P0001',role='reference',sequence='CAG',sha256=sha256(b'CAG').hexdigest(),
                 locus='L',assembly='GRCh38',observations=[],molecules=[])
        for interval in (['chr1',10,13],dict(chrom='chr1',start=10,end=13)):
            row=dict(ref,reference_interval=interval)
            self.assertEqual(member_identity([row],model,'CAG')['status'],'verified')
            with self.assertRaisesRegex(ValueError,'reference bases'):
                member_identity([row],model,'CAA')
            with self.assertRaisesRegex(ValueError,'assembly'):
                member_identity([dict(row,assembly='GRCh37')],model)
        self.assertEqual(member_identity([ref],model)['status'],'partially_verified')
        self.assertEqual(member_identity([],model,'CAG')['status'],'partially_verified')

    def test_wrong_locus_fails_before_output(self):
        model = ROOT / 'catalogs/D1S80.structure.v1.json'
        ref = 'CAG' * 23
        row = dict(id='P0001', role='reference', sequence=ref, sha256=sha256(ref.encode()).hexdigest(),
                   reference_interval=['chr12', 111598950, 111599019], observations=[], molecules=[], samples={})
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            source = base / 'members.json'
            source.write_text(json.dumps([row]))
            with self.assertRaisesRegex(ValueError, 'identity.*(chrom|interval|length)'):
                run(source, model, base / 'out', render=False)
            self.assertFalse((base / 'out').exists())

    def test_zero_core_keeps_both_cuts(self):
        from types import SimpleNamespace
        rng = random.Random(17)
        left, right = [''.join(rng.choice('ACGT') for _ in range(250)) for _ in range(2)]
        dna = left + 'CAG'*10 + right
        locus = Locus('X','chr1',250,280,(),boundary_status='confirmed')
        reference=SimpleNamespace(fetch=lambda c,a,b:dna[a:b],contig_length=lambda c:len(dna))
        context=build_context(locus,reference)
        for middle, cigar, status in [('', '250M30D250M', 'complete_empty'), ('CAG','253M27D250M','success'),
                                     ('C','251M29D250M','success')]:
            fields=['x','0','chr1','1','60',cigar,'*','0','0',left+middle+right,'*']
            f=extract_record(fields,locus,context)
            self.assertEqual((f['qstart'],f['qend'],f['sequence'],f['core_status']),
                             (250,250+len(middle),middle,status))
        old=build_context(locus,reference,PREVIOUS_VERSION)
        fields[5]='250M30D250M';fields[9]=left+right
        self.assertEqual(len(extract_record(fields,locus,old)['sequence']),1)
        self.assertNotEqual(old.model_hash,context.model_hash)
        pairs=pair_endpoints([[250,0,250,250]],[[250,0,250,250],[251,0,250,250]],[250,250],allow_empty=True)
        fake=dict(status='success',pairs=pairs,representative=pairs[0],ambiguous=True,
                  left={'endpoints':[[250,0,250,250]]},right={'endpoints':[[250,0,250,250]]})
        with patch.object(context,'locate',return_value=fake):
            mixed=context.extract({},left+right,None,[], 'confirmed')
        self.assertEqual(mixed['core_status'],'empty_or_nonempty')
        with tempfile.TemporaryDirectory() as temp:
            rec=dict(mixed,record_id='r1')
            row=dict(id='R1',molecule_key='M',sample='S',platform='HiFi',read_id='Q',records=[rec,dict(rec,record_id='r2')])
            counts=write_flank_outputs(Path(temp),[row])
            self.assertEqual(counts['empty_or_nonempty'],dict(records=2,molecules=1))

    def test_iupac_in_both_hmm_passes(self):
        dictionary=[dict(id='U01',base='U01',sequence='ACGTACGT')]
        for assessed in (False, True):
            tokens,audit=segment_hmm(('ACNT',dictionary,'HiFi'),return_details=True,assess_structure=assessed)
            self.assertEqual((tokens[0]['exact_bases'],tokens[0]['status'],audit['decoder']),
                             ('ACNT','unknown','opaque_iupac'))
            with self.assertRaises(ValueError):
                segment_hmm(('AC!T',dictionary,'HiFi'),assess_structure=assessed)

    @unittest.skipUnless(os.environ.get('SAMTOOLS') or shutil.which('samtools'),'samtools is required for BAM integration')
    def test_all_empty_bam_exports_reference_and_counts(self):
        samtools=os.environ.get('SAMTOOLS') or shutil.which('samtools')
        rng=random.Random(42)
        left,right=[''.join(rng.choice('ACGT') for _ in range(250)) for _ in range(2)]
        model=json.loads((ROOT/'catalogs/D1S80.structure.v1.json').read_bytes())
        core=model['motifs'][0]['sequence']*2
        model.update(start=250,end=250+len(core))
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);ref=root/'ref.fa';bed=root/'locus.bed';mp=root/'model.json'
            ref.write_bytes(('>chr1\n'+left+core+right+'\n').encode())
            bed.write_text(f"chr1\t250\t{model['end']}\tID=D1S80;BOUNDARY=confirmed\n")
            mp.write_text(json.dumps(model))
            bam=root/'bams/D1S80/test/HiFi/test.D1S80.bam';bam.parent.mkdir(parents=True)
            sam=(f'@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{500+len(core)}\n'
                 f'z\t0\tchr1\t1\t60\t250M{len(core)}D250M\t*\t0\t0\t{left+right}\t*\n')
            subprocess.run([samtools,'faidx',str(ref)],check=True)
            subprocess.run([samtools,'view','-b','-o',str(bam),'-'],input=sam.encode(),check=True)
            subprocess.run([samtools,'index',str(bam)],check=True)
            cmd=[sys.executable,str(ROOT/'scripts/vntrpathgraph_cli.py'),'--bam-dir',str(root/'bams'),
                 '--reference',str(ref),'--catalog',str(bed),'--locus','D1S80','--locus-model',str(mp),
                 '--graph-mode','structure','--supernode-policy',str(ROOT/'catalogs/D1S80.supernodes.v1.json'),
                 '--samtools',samtools,'--out',str(root/'out')]
            subprocess.run(cmd,check=True,stdout=subprocess.PIPE)
            out=root/'out/D1S80';summary=json.loads((out/'summary.json').read_bytes())
            self.assertEqual((summary['paths'],summary['molecules']),(1,0))
            self.assertEqual(summary['core_length_evidence']['complete_empty'],dict(records=1,molecules=1))
            self.assertEqual(json.loads((out/'extraction_audit.jsonl').read_bytes())['core_status'],'complete_empty')
            self.assertEqual(json.loads((out/'source_members.json').read_bytes())[0]['sequence'],core)

    def test_export_q_matrix_preserves_unknown_and_binds_delivery(self):
        from test_supernode_export import fixture
        from vntrsupernode_export import export_supernodes, validate_bundle, digest
        model,records=fixture()
        locus=Locus(model['locus'],model['chrom'],model['start'],model['end'],(),boundary_status='confirmed')
        row=records[-1]
        row['sequence']=row['sequence'][:9]+'N'+row['sequence'][10:]
        row['sha256']=sha256(row['sequence'].encode()).hexdigest()
        with tempfile.TemporaryDirectory() as temp:
            for q in (None,0,.2):
                out=Path(temp)/str(q);out.mkdir()
                export_supernodes(out,locus,records,ROOT/'catalogs/D1S80.structure.v1.json',
                                  ROOT/'catalogs/D1S80.supernodes.v1.json',render=False,q=q)
                features=json.loads((out/'annotation_audit.json').read_bytes())
                self.assertEqual(features[row['id']]['annotation']['decoder'],'opaque_iupac')
                self.assertEqual(''.join(t['exact_bases'] for t in features[row['id']]['tokens']),row['sequence'])
                self.assertEqual(validate_bundle(out,ROOT)['paths'],3)
            txt=out/'D1S80.path_sequences.txt'
            txt.write_bytes(txt.read_bytes().replace(row['sequence'].encode(),b'CAG'))
            prov_path=out/'provenance.json'
            prov=json.loads(prov_path.read_bytes());prov['files'][txt.name]=digest(txt)
            prov_path.write_text(json.dumps(prov))
            with self.assertRaisesRegex(ValueError,'TXT path content'):
                validate_bundle(out)

    def test_imported_terminal_obeys_path_sample_scope(self):
        model_path=ROOT/'catalogs/D1S80.structure.v1.json'
        model=json.loads(model_path.read_bytes())
        locus=Locus(model['locus'],model['chrom'],model['start'],model['end'],(),boundary_status='confirmed')
        _,dictionary=load_model(model_path,locus)
        term=model['terminal_templates'][0]
        sample=term['samples'][0]
        seq=model['motifs'][0]['sequence']+term['sequence']
        frag=dict(id='F1',sequence=seq)
        rows=[dict(sample=sample,platform='HiFi',fragments=[frag])]
        annotate_rows(rows,dictionary,1)
        rec=dict(id='P0002',role='observation',sequence=seq,samples={sample:1},
                 observations=[dict(sample=sample,platform='HiFi',fragment='F1')])
        frozen=dict(model_sha256=sha256(model_path.read_bytes()).hexdigest())
        direct=annotate_paths([rec],model_path,locus)[2]['P0002']
        imported=annotate_paths([rec],model_path,locus,rows,frozen)[2]['P0002']
        self.assertEqual(direct['tokens'],imported['tokens'])
        self.assertEqual(imported['tokens'][-1].get('role'),'right_terminal')
        self.assertEqual(''.join(t['exact_bases'] for t in imported['tokens']),seq)
        rec['samples']['outside-scope']=1
        mixed=annotate_paths([rec],model_path,locus,rows,frozen)[2]['P0002']['tokens']
        self.assertFalse(any(t.get('role')=='right_terminal' for t in mixed))
        rec['samples']={sample:1}
        annotate_rows(rows,dictionary,1,model=model)
        self.assertEqual(frag['tokens'][-1]['role'],'right_terminal')
        self.assertEqual(frag['annotation']['ambiguity_scope'],'all_minimum_cost_paths_in_seed_candidate_DAG')
        self.assertEqual(annotate_paths([rec],model_path,locus,rows,frozen)[2]['P0002']['tokens'],direct['tokens'])
        rec['samples']['outside-scope']=1
        self.assertFalse(any(t.get('role')=='right_terminal'
                         for t in annotate_paths([rec],model_path,locus,rows,frozen)[2]['P0002']['tokens']))
