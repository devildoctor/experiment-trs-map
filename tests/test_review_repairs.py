"""Regression cases from the November 22 point-graph review."""
import json
from hashlib import sha256
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from run_efficient_family import run
from vntrflankalign import build_context
from vntrgraph import Locus
from vntrhmm import segment_hmm
from vntrreads import extract_record
from vntrpathstructure import annotate_paths
from vntrpartition import load_model
from vntrstructure import annotate_rows


class ReviewRepairs(unittest.TestCase):
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
        context = build_context(locus, SimpleNamespace(fetch=lambda c,a,b:dna[a:b],contig_length=lambda c:len(dna)))
        for middle, cigar, status in [('', '250M30D250M', 'complete_empty'), ('CAG','253M27D250M','success'),
                                     ('C','251M29D250M','success')]:
            fields=['x','0','chr1','1','60',cigar,'*','0','0',left+middle+right,'*']
            f=extract_record(fields,locus,context)
            self.assertEqual((f['qstart'],f['qend'],f['sequence'],f['core_status']),
                             (250,250+len(middle),middle,status))

    def test_iupac_in_both_hmm_passes(self):
        dictionary=[dict(id='U01',base='U01',sequence='ACGTACGT')]
        for assessed in (False, True):
            tokens,audit=segment_hmm(('ACNT',dictionary,'HiFi'),return_details=True,assess_structure=assessed)
            self.assertEqual((tokens[0]['exact_bases'],tokens[0]['status'],audit['decoder']),
                             ('ACNT','unknown','opaque_iupac'))
            with self.assertRaises(ValueError):
                segment_hmm(('AC!T',dictionary,'HiFi'),assess_structure=assessed)

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
