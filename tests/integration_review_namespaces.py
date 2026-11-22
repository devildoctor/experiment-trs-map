"""Exercise namespace identity through actual BAMs and both view CLIs."""
import csv
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def main():
    samtools = sys.argv[1]
    source = Path(__file__).resolve().parents[1]

    def run(*args):
        subprocess.run([str(arg) for arg in args], check=True)

    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        rng = random.Random(193)
        motif = 'GAGGACCACCAGGAAG'
        sequence = (''.join(rng.choices('ACGT', k=200)) + motif * 4
                    + ''.join(rng.choices('ACGT', k=200)))
        reference = root / 'reference.fa'
        reference.write_text('>chrT\n' + sequence + '\n', encoding='utf-8')
        run(samtools, 'faidx', reference)
        catalog = root / 'loci.bed'
        catalog.write_text(f'chrT\t200\t264\tID=TESTVNTR;MOTIFS={motif}\n', encoding='utf-8')
        bam_dir = root / 'bam' / 'S1' / 'HiFi'
        bam_dir.mkdir(parents=True)
        bams = []
        # 同名读段来自两个独立文库，另一个 BAM 重复保存 A 文库分子。
        for label, read_groups in [('original', ['A', 'B']), ('copy', ['A'])]:
            lines = ['@HD\tVN:1.6\tSO:unsorted', f'@SQ\tSN:chrT\tLN:{len(sequence)}',
                     '@RG\tID:A\tLB:library-A\tSM:S1', '@RG\tID:B\tLB:library-B\tSM:S1']
            for group in read_groups:
                bases = sequence[100:364]
                lines.append(f'same-name\t0\tchrT\t101\t60\t264M\t*\t0\t0\t{bases}\t'
                             + 'I' * len(bases) + f'\tRG:Z:{group}')
            sam = root / f'{label}.sam'
            sam.write_text('\n'.join(lines) + '\n', encoding='utf-8')
            bam = bam_dir / f'{label}.TESTVNTR.bam'
            run(samtools, 'sort', '-o', bam, sam)
            run(samtools, 'index', bam)
            bams.append(bam)
        sources = root / 'sources.tsv'
        sources.write_text('bam\tnamespace\n' + ''.join(
            f'{bam.relative_to(root).as_posix()}\tshared\n' for bam in bams), encoding='utf-8')
        for mode, suffix in [('read', 'reads'), ('structure', 'structures')]:
            for override in [False, True]:
                output = root / f'{mode}-{override}'
                command = [sys.executable, source / 'scripts' / 'vntrgraph.py', f'--{mode}-view',
                           '--bam-dir', root / 'bam', '--reference', reference, '--catalog', catalog,
                           '--out', output, '--locus', 'TESTVNTR', '--samtools', samtools, '--workers', '1']
                if override:
                    command.extend(['--sources', sources])
                run(*command)
                locus_dir = output / 'TESTVNTR'
                data = json.loads((locus_dir / f'TESTVNTR.{suffix}.json').read_text(encoding='utf-8'))
                rows = data['rows']
                expected_molecules = 1 if override else 2
                assert len(rows) == (2 if override else 3), rows
                assert len({row['id'] for row in rows}) == len(rows)
                keys = {json.dumps(row['molecule_key'], sort_keys=True) for row in rows}
                assert len(keys) == expected_molecules, keys
                assert {row['namespace'] for row in rows} == ({'shared'} if override else {'library-A', 'library-B'})
                assert data['audit']['alignment_records'] == 3
                assert sum(len(row['fragments']) for row in rows) == 3
                for row in rows:
                    for fragment in row['fragments']:
                        assert fragment['sequence'] == motif * 4
                        assert ''.join(token['exact_bases'] for token in fragment['tokens']) == fragment['sequence']
                with (locus_dir / 'reads.tsv').open(encoding='utf-8', newline='') as handle:
                    table = list(csv.DictReader(handle, delimiter='\t'))
                assert len(table) == len(rows) and all('molecule_key' in row and 'namespace' in row for row in table)
                if mode == 'structure':
                    node = next(node for node in data['graph']['structure_nodes'] if node['id'] == 'M1')
                    assert node['support'] == expected_molecules, node
                    assert node['occurrences'] == 12, node
                    run(sys.executable, source / 'tests' / 'verify_structure_results.py', locus_dir)
        print('NAMESPACE_BAM_INTEGRATION_PASS: read/structure, independent LB, BAM copies, relative sources override and exact slices')


if __name__ == '__main__':
    main()
