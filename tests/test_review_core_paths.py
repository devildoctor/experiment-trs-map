"""审查问题回归：主/备选歧义对称，以及位点输出目录不能越过指定根目录。"""
from pathlib import Path
from types import SimpleNamespace
import os
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from vntrcore import decompose
import vntrgraph as app


class CoreAmbiguityTests(unittest.TestCase):
    def test_tied_alternative_marks_its_different_token(self):
        primary, _, audit = decompose('ACGT', ['ACAT', 'ACCT'], return_details=True)
        self.assertTrue(audit['annotation_ambiguity'])
        self.assertEqual(audit['score_gap'], 0)
        self.assertNotEqual(primary[0].family, audit['alternative'][0]['family'])
        self.assertTrue(primary[0].annotation_ambiguity)
        self.assertTrue(audit['alternative'][0]['annotation_ambiguity'])
        self.assertEqual(''.join(t.exact_bases for t in primary), 'ACGT')
        self.assertEqual(''.join(t['exact_bases'] for t in audit['alternative']), 'ACGT')

    def test_shared_token_stays_clear_in_both_tied_paths(self):
        sequence = 'TTGCACGT'
        primary, _, audit = decompose(sequence, ['TTGC', 'ACAT', 'ACCT'], return_details=True)
        self.assertTrue(audit['annotation_ambiguity'])
        self.assertEqual([t.annotation_ambiguity for t in primary], [False, True])
        self.assertEqual([t['annotation_ambiguity'] for t in audit['alternative']], [False, True])
        self.assertEqual(primary[0].exact_bases, 'TTGC')
        self.assertEqual(audit['alternative'][0]['exact_bases'], 'TTGC')

    def test_non_tied_alternative_does_not_mark_ambiguity(self):
        primary, _, audit = decompose('ACAT', ['ACAT', 'ACCT'], return_details=True)
        self.assertFalse(audit['annotation_ambiguity'])
        self.assertGreater(audit['score_gap'], 0)
        self.assertTrue(audit['alternative'])
        self.assertFalse(any(t.annotation_ambiguity for t in primary))
        self.assertFalse(any(t['annotation_ambiguity'] for t in audit['alternative']))


class LocusPathTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='review-core-', dir=ROOT)
        self.folder = Path(self.temporary.name).resolve()
        self.folder.relative_to(ROOT.resolve())
        self.addCleanup(self.temporary.cleanup)

    def catalog(self, name):
        path = self.folder / 'loci.bed'
        path.write_text(f'chr1\t100\t112\tID={name};MOTIFS=ACGT\n', encoding='utf-8')
        return path

    def symlink(self, target, link):
        try:
            link.symlink_to(target, target_is_directory=True)
        except OSError as error:
            if os.name == 'nt':
                self.skipTest(f'Windows symlink creation unavailable: {error}')
            raise

    def test_catalog_rejects_special_path_components(self):
        for name in ('.', '..'):
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, 'locus ID'):
                    app.read_catalog(self.catalog(name))

    def test_catalog_accepts_normal_id_with_punctuation(self):
        locus = app.read_catalog(self.catalog('D1S80.v2_sample-1'))[0]
        self.assertEqual(locus.name, 'D1S80.v2_sample-1')

    def test_helper_checks_id_even_without_catalog(self):
        output = self.folder / 'output'
        for name in ('.', '..', '../escape', 'a/b', 'a\\b', ''):
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, 'locus ID'):
                    app.locus_output_dir(output, name)
        self.assertFalse(output.exists())

    def test_helper_accepts_normal_child_without_creating_it(self):
        output = self.folder / 'output'
        destination = app.locus_output_dir(output, 'D1S80')
        self.assertEqual(destination, output.resolve() / 'D1S80')
        self.assertFalse(destination.exists())

    def test_helper_rejects_symlink_outside_output_root(self):
        output = self.folder / 'output'
        outside = self.folder / 'outside'
        output.mkdir()
        outside.mkdir()
        self.symlink(outside, output / 'D1S80')
        with self.assertRaisesRegex(ValueError, 'outside output root'):
            app.locus_output_dir(output, 'D1S80')
        self.assertEqual(list(outside.iterdir()), [])

    def test_helper_rejects_symlink_to_output_root_itself(self):
        output = self.folder / 'output'
        output.mkdir()
        self.symlink(output, output / 'D1S80')
        with self.assertRaisesRegex(ValueError, 'outside output root'):
            app.locus_output_dir(output, 'D1S80')

    def test_helper_allows_symlink_that_stays_inside_output_root(self):
        output = self.folder / 'output'
        target = output / 'target'
        target.mkdir(parents=True)
        self.symlink(target, output / 'D1S80')
        self.assertEqual(app.locus_output_dir(output, 'D1S80'), target.resolve())

    def test_direct_process_rejects_special_id_before_writing(self):
        output = self.folder / 'output'
        locus = app.Locus('..', 'chr1', 100, 112, ('ACGT',))
        with self.assertRaisesRegex(ValueError, 'locus ID'):
            app.process_locus(locus, SimpleNamespace(out=output), None, {})
        self.assertFalse(output.exists())

    def test_direct_process_rejects_escaped_symlink_before_writing(self):
        output = self.folder / 'output'
        outside = self.folder / 'outside'
        output.mkdir()
        outside.mkdir()
        self.symlink(outside, output / 'D1S80')
        locus = app.Locus('D1S80', 'chr1', 100, 112, ('ACGT',))
        with self.assertRaisesRegex(ValueError, 'outside output root'):
            app.process_locus(locus, SimpleNamespace(out=output), None, {})
        self.assertEqual(list(outside.iterdir()), [])


if __name__ == '__main__':
    unittest.main()
