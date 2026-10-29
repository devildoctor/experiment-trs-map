"""Regression checks for empty evidence and explicit dBG support thresholds."""
from collections import Counter
from contextlib import redirect_stdout
import csv
from hashlib import sha256
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
import vntrdbg
import vntrdbg_inputs
import vntrdbg_report
from vntrdbg_analysis import JOINT_MOLECULE_FIELDS, MOLECULE_EVENT_FIELDS
from validate_vntrdbg_analysis import validate_joint


def record(name, sequence, kind='observed', sample='S'):
    row = dict(id=name, sequence=sequence, kind=kind, sample=sample,
               sequence_sha256=sha256(sequence.encode()).hexdigest(), platforms=['HiFi'],
               coverage='two_flank_anchored' if kind == 'observed' else 'candidate_sequence')
    if kind == 'observed':
        row['molecule_key'] = sample + '/' + name
    return row


def write_inputs(root, records, motif):
    dbg, point = root/'dbg', root/'point'
    dbg.mkdir()
    point.mkdir()
    (point/'motifs.tsv').write_text('motif\tusable\tlevel\n'+motif+'\tTrue\tbase\n', encoding='utf-8')
    (point/'tokens.tsv').write_text('allele\tstart\tend\tstatus\n', encoding='utf-8')
    (dbg/'input_paths.json').write_text(json.dumps(records), encoding='utf-8')
    summary = dict(locus='D1S80', path_types=dict(Counter(r['kind'] for r in records)),
                   source_evidence_hashes={p.name: sha256(p.read_bytes()).hexdigest() for p in point.iterdir()})
    (dbg/'input_summary.json').write_text(json.dumps(summary), encoding='utf-8')
    return dbg, point


def analysis_cli(dbg, point, out, *extra):
    return subprocess.run([sys.executable, '-X', 'utf8', str(SCRIPTS/'vntrdbg_analysis.py'),
                           '--dbg-dir', str(dbg), '--pointgraph-dir', str(point), '--out', str(out),
                           '--k', '3', *extra], capture_output=True, text=True, encoding='utf-8')


def read_rows(path):
    with path.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return reader.fieldnames, list(reader)


class EmptyEvidenceTests(unittest.TestCase):
    def test_cli_no_observations_for_candidate_sample_writes_complete_results(self):
        motif = 'CCTGGTGGTCCTCCTT'
        records = [record('REF', motif, 'reference', ''), record('CAND_C', motif, 'candidate'),
                   record('OBS_T', motif, sample='T')]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dbg, point = write_inputs(root, records, motif)
            result = analysis_cli(dbg, point, root/'out')
            self.assertEqual(result.returncode, 0, result.stderr)
            for name, expected in [('joint_molecules.tsv', JOINT_MOLECULE_FIELDS),
                                   ('molecule_events.tsv', MOLECULE_EVENT_FIELDS)]:
                fields, rows = read_rows(root/'out'/name)
                self.assertEqual(fields, expected)
                self.assertEqual(rows, [])
            _, rows = read_rows(root/'out/joint_support.tsv')
            self.assertTrue(rows)
            for row in rows:
                self.assertEqual(row['tested_molecules'], '0')
                self.assertEqual(row['support_molecules'], '0')
                self.assertEqual(row['conflict_molecules'], '0')
                self.assertEqual(row['uninformative_molecules'], '0')
                self.assertEqual(row['support_fraction'], 'NA')
                self.assertEqual(row['status'], 'no_observed_molecules')
            manifest = json.loads((root/'out/manifest.json').read_text())
            self.assertEqual(manifest['decomposed_path_count'], 3)
            self.assertEqual(manifest['branch_min_support'], 3)
            self.assertEqual(set(p.name for p in (root/'out').iterdir()),
                             {'joint_queries.tsv', 'events.tsv', 'joint_support.tsv', 'joint_molecules.tsv',
                              'molecule_events.tsv', 'primitive_units.tsv', 'unit_summary.tsv',
                              'motif_counts.tsv', 'manifest.json'})

    def test_cli_mixed_candidate_samples_keeps_supported_and_empty_summaries(self):
        motif = 'CCTGGTGGTCCTCCTT'
        records = [record('REF', motif, 'reference', ''), record('CAND_C', motif, 'candidate'),
                   record('CAND_T', motif, 'candidate', 'T'), record('OBS_S', motif)]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dbg, point = write_inputs(root, records, motif)
            result = analysis_cli(dbg, point, root/'out')
            self.assertEqual(result.returncode, 0, result.stderr)
            fields, molecules = read_rows(root/'out/joint_molecules.tsv')
            self.assertEqual(fields, JOINT_MOLECULE_FIELDS)
            self.assertTrue(molecules)
            self.assertEqual({r['sample'] for r in molecules}, {'S'})
            _, summaries = read_rows(root/'out/joint_support.tsv')
            full = {r['sample']: r for r in summaries if r['kind'] == 'full_candidate_exact'}
            self.assertEqual(full['S']['support_molecules'], '1')
            self.assertEqual(full['S']['status'], 'computed')
            self.assertEqual(full['T']['status'], 'no_observed_molecules')


class ThresholdTests(unittest.TestCase):
    def test_joint_cli_threshold_changes_observed_branch_definition(self):
        target, alternate = 'AAACCGTATGCTAGTC', 'AAATCGTATGCTAGTC'
        records = [record('CAND_C', target, 'candidate'),
                   record('OBS_1', alternate), record('OBS_2', alternate)]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dbg, point = write_inputs(root, records, target)
            for threshold in (2, 3):
                out = root/f'out{threshold}'
                result = analysis_cli(dbg, point, out, '--branch-min-support', str(threshold))
                self.assertEqual(result.returncode, 0, result.stderr)
                manifest = json.loads((out/'manifest.json').read_text())
                self.assertEqual(manifest['branch_min_support'], threshold)
                self.assertEqual(manifest['branch_support_scope'], 'max_unique_molecules_within_one_sample')
                _, queries = read_rows(out/'joint_queries.tsv')
                branches = [q for q in queries if q['kind'] == 'branch_block']
                self.assertEqual(bool(branches), threshold == 2)
                self.assertEqual(validate_joint(records, manifest, out)[0], len(queries))
                legacy = {key: value for key, value in manifest.items() if key != 'branch_min_support'}
                if threshold == 3:
                    self.assertEqual(validate_joint(records, legacy, out)[0], len(queries))
                else:
                    with self.assertRaisesRegex(AssertionError, 'query set is incomplete'):
                        validate_joint(records, legacy, out)
                invalid = dict(manifest, branch_min_support=0)
                with self.assertRaisesRegex(AssertionError, 'invalid branch support threshold'):
                    validate_joint(records, invalid, out)

    def test_invalid_joint_threshold_rejected_before_output_creation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = analysis_cli(root/'missing', root/'missing', root/'out', '--branch-min-support', '0')
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('--branch-min-support must be >= 1', result.stderr)
            self.assertFalse((root/'out').exists())

    def test_edge_color_uses_per_sample_maximum_and_selected_threshold(self):
        records = [record('REF', 'AAC', 'reference', ''), record('OBS_A', 'AACGT', sample='A'),
                   record('OBS_B1', 'AACGT', sample='B'), record('OBS_B2', 'AACGT', sample='B')]
        with tempfile.TemporaryDirectory() as temporary:
            for threshold, expected_color in [(3, '#b8c0ca'), (2, '#47566e')]:
                graph = vntrdbg.compact_graph(records, vntrdbg.raw_graph(records, 3), 3,
                                              'TEST', 'all_evidence', threshold)
                self.assertTrue(graph['edges'])
                self.assertTrue(all(e['molecule_support'] == 3 and e['max_sample_support'] == 2
                                    for e in graph['edges']))
                target = Path(temporary)/f'k3_{threshold}.svg'
                with patch.object(vntrdbg.subprocess, 'run'):
                    vntrdbg.topology(target, graph)
                dot = target.with_suffix('.dot').read_text()
                edge_lines = [line for line in dot.splitlines() if ' -> ' in line]
                self.assertTrue(all('color="'+expected_color+'"' in line for line in edge_lines))
                self.assertIn(f'min_support={threshold}', dot)
                self.assertIn('max_sample_support=2', dot)

    def test_graph_cli_propagates_threshold_to_both_views_and_metrics(self):
        records = [record('REF', 'AAC', 'reference', ''), record('OBS_1', 'AACGT'), record('OBS_2', 'AACGT')]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            argv = ['vntrdbg.py', '--evidence-dir', str(root), '--bam-dir', str(root),
                    '--reference', str(root/'reference.fa'), '--catalog', str(root/'catalog'),
                    '--models', str(root/'models'), '--locus', 'TEST', '--samtools', 'samtools',
                    '--out', str(root/'out'), '--k', '3', '4', '--min-support', '2']
            source = SimpleNamespace(prepare=lambda args: (records, {}))
            with patch.object(sys, 'argv', argv), patch.dict(sys.modules, {'vntrdbg_inputs': source}), \
                    patch.object(vntrdbg, 'topology'), redirect_stdout(io.StringIO()):
                vntrdbg.main()
            for k in (3, 4):
                folder = root/'out'/f'k{k}'
                for suffix in ('', '.supported'):
                    graph = json.loads((folder/f'TEST.k{k}{suffix}.json').read_text())
                    self.assertEqual(graph['min_support'], 2)
                    self.assertEqual(graph['support_scope'], 'max_unique_molecules_within_one_sample')
                    fields, _ = read_rows(folder/f'TEST.k{k}{suffix}.support.tsv')
                    self.assertEqual(fields, ['object','id','molecule_support','sample_support',
                        'platform_support','molecule_ids','origins','max_sample_support'])
                self.assertEqual(json.loads((folder/'graph_semantics.json').read_text())['min_support'], 2)
                self.assertEqual(json.loads((folder/'metrics.json').read_text())['min_support'], 2)
            self.assertTrue(all(r['min_support'] == 2 for r in
                                json.loads((root/'out/k_comparison.json').read_text())))


class ReportThresholdTests(unittest.TestCase):
    def test_metadata_legacy_text_and_missing_parameter_are_distinguished(self):
        self.assertEqual(vntrdbg_report.support_threshold({'min_support': 2}, {})[0], 2)
        self.assertEqual(vntrdbg_report.support_threshold({}, {'min_support': 1})[0], 1)
        value, source = vntrdbg_report.support_threshold({}, {'supported_view':
            'at least 2 molecules within one sample per raw node/edge; reference retained explicitly'})
        self.assertEqual(value, 2)
        self.assertIn('recorded legacy', source)
        value, source = vntrdbg_report.support_threshold({}, {})
        self.assertEqual(value, 3)
        self.assertIn('inferred legacy default', source)
        with self.assertRaises(ValueError):
            vntrdbg_report.support_threshold({'min_support': 0}, {})

    def test_report_uses_each_k_threshold_and_marks_legacy_default(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            folder = run/'results/dbg/D1S80'
            folder.mkdir(parents=True)
            acceptance = run/'results/acceptance'
            acceptance.mkdir()
            (acceptance/'str-regression.json').write_text('{"status":"PASS"}')
            (run/'results/preparation').mkdir()
            (run/'source').mkdir()
            (run/'source/SHA256SUMS').write_text('synthetic test\n')
            (run/'str-commit.txt').write_text('synthetic test\n')
            (folder/'validation.json').write_text('{"status":"PASS"}')
            (folder/'input_summary.json').write_text(json.dumps(dict(bams=0,input_records=0,
                independent_molecules=0,observation_scope={})))
            comparison=[]
            for k, threshold in [(3, 1), (4, 2), (5, None)]:
                d=folder/f'k{k}'
                d.mkdir()
                graph={'paths': []}
                if threshold is not None:
                    graph['min_support']=threshold
                (d/f'D1S80.k{k}.supported.json').write_text(json.dumps(graph))
                comparison.append(dict(k=k,compact_nodes=0,compact_edges=0,supported_compact_nodes=0,
                    supported_compact_edges=0,supported_branch_nodes=0,candidate_edge_retention=0,
                    candidate_repeated_kmer_fraction=0))
            (folder/'k_comparison.json').write_text(json.dumps(comparison))
            with patch.object(vntrdbg_report, 'preview'), redirect_stdout(io.StringIO()):
                vntrdbg_report.main(run)
            report=(run/'RESULTS.md').read_text()
            self.assertIn('k=3: 1 (recorded graph metadata)', report)
            self.assertIn('k=4: 2 (recorded graph metadata)', report)
            self.assertIn('k=5: 3 (inferred legacy default; parameter not recorded)', report)
            self.assertNotIn('require three molecules', report)


class SourceNamespaceTests(unittest.TestCase):
    def prepare_fixture(self, root, mapping):
        sequence = 'CCTGGTGGTCCTCCTT'
        evidence, out = root/'evidence', root/'out'
        evidence.mkdir()
        out.mkdir()
        reference, catalog, models, bam = root/'ref.fa', root/'loci.bed', root/'models.tsv', root/'a.bam'
        reference.write_bytes(('>chr1\n'+sequence+'\n').encode())
        index = Path(str(reference)+'.fai')
        index.write_text('chr1\t16\t6\t16\t17\n', encoding='utf-8')
        catalog.write_text('chr1\t0\t16\tID=TEST\n', encoding='utf-8')
        models.write_text('locus\tmodel\nTEST\tcompact_dbg\n', encoding='utf-8')
        bam.write_bytes(b'synthetic SAM transport fixture')
        sources = root/'sources.tsv'
        if mapping is not None:
            sources.write_text('bam\tnamespace\n'+mapping, encoding='utf-8')
        parameters = dict(catalog_sha256=sha256(catalog.read_bytes()).hexdigest(),
            reference_fai_sha256=sha256(index.read_bytes()).hexdigest(),
            reference_size=reference.stat().st_size, reference_mtime_ns=reference.stat().st_mtime_ns,
            sources=str(sources) if mapping is not None else 'None')
        vntrdbg.tsv(evidence/'parameters.tsv', [dict(parameter=k,value=v) for k,v in parameters.items()])
        vntrdbg.tsv(evidence/'bam_audit.tsv', [dict(bam=str(bam))])
        for name, fields in [('molecules.tsv',['molecule_key']), ('audit.tsv',['allele']),
                             ('tokens.tsv',['allele']), ('motifs.tsv',['motif','usable','level'])]:
            vntrdbg.tsv(evidence/name, [], fields)
        (evidence/'anchor_decisions.jsonl').write_text('', encoding='utf-8')
        args = SimpleNamespace(models=models,locus='TEST',catalog=catalog,evidence_dir=evidence,
                               reference=reference,bam_dir=root,samtools='samtools',out=out)
        sam = '@SQ\tSN:chr1\tLN:16\n@RG\tID:rg1\tLB:library\n'
        sam += 'read\t0\tchr1\t1\t60\t16M\t*\t0\t0\t'+sequence+'\t*\tRG:Z:rg1\n'
        return args, bam, sam

    def test_prepare_relative_sources_and_empty_namespace_override_library(self):
        for namespace in (None, '', 'chosen'):
            with self.subTest(namespace=namespace), tempfile.TemporaryDirectory() as temporary:
                mapping = None if namespace is None else './a.bam\t'+namespace+'\n'
                args, bam, sam = self.prepare_fixture(Path(temporary), mapping)
                with patch.object(vntrdbg_inputs, 'discover_bams', return_value=[('S','HiFi',bam)]), \
                        patch.object(vntrdbg_inputs.subprocess, 'run'), \
                        patch.object(vntrdbg_inputs.subprocess, 'check_output', return_value=sam), \
                        redirect_stdout(io.StringIO()):
                    records, _ = vntrdbg_inputs.prepare(args)
                observed = [r for r in records if r['kind']=='observed']
                self.assertEqual(len(observed), 1)
                self.assertEqual(json.loads(observed[0]['molecule_key']),
                                 ['S', 'library' if namespace is None else namespace, 'read'])

    def test_prepare_rejects_duplicate_normalized_bam_before_writing(self):
        with tempfile.TemporaryDirectory() as temporary:
            args, bam, _ = self.prepare_fixture(Path(temporary), 'a.bam\tfirst\n./a.bam\tother\n')
            with patch.object(vntrdbg_inputs, 'discover_bams', return_value=[('S','HiFi',bam)]):
                with self.assertRaisesRegex(ValueError, 'Duplicate BAM'):
                    vntrdbg_inputs.prepare(args)
            self.assertEqual(list(args.out.iterdir()), [])


if __name__ == '__main__':
    unittest.main()
