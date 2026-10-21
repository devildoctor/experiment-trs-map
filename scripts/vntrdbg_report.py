"""Package completed server-side results and render bounded inspection previews."""
from pathlib import Path
import csv
import hashlib
import json
import subprocess
import struct
import sys
import tarfile
import xml.etree.ElementTree as ET


def preview(source, target, width, height):
    # librsvg 2.42 can exceed the requested box when both dimensions use -a.
    # Constrain only the limiting dimension, preserving the SVG aspect ratio.
    box=[float(x) for x in ET.parse(source).getroot().attrib['viewBox'].split()]
    limit=['-w',str(width)] if box[2]/box[3]>=width/height else ['-h',str(height)]
    subprocess.run(['rsvg-convert',*limit,str(source),'-o',str(target)],check=True)
    pixels=struct.unpack('>II',target.read_bytes()[16:24])
    assert 0<pixels[0]<=width and 0<pixels[1]<=height,(target,pixels)


def main(run):
    folder=run/'results/dbg/D1S80'
    validation=json.loads((folder/'validation.json').read_text())
    str_check=json.loads((run/'results/acceptance/str-regression.json').read_text())
    assert validation['status']==str_check['status']=='PASS'
    comparison=json.loads((folder/'k_comparison.json').read_text())
    inputs=json.loads((folder/'input_summary.json').read_text())
    review=run/'review'
    review.mkdir()
    for row in comparison:
        k=row['k']
        d=folder/f'k{k}'
        for suffix in ('','.supported'):
            stem=f'D1S80.k{k}{suffix}'
            preview(d/(stem+'.svg'),review/(stem+'.png'),1800,1400)
        graph=json.loads((d/f'D1S80.k{k}.supported.json').read_text())
        candidate=next((p for p in graph['paths'] if p['kind']=='candidate' and p.get('filter_complete',True)),None)
        if candidate:
            source=d/'paths-svg.supported'/f'D1S80.k{k}.supported.{candidate["id"]}.svg'
            preview(source,review/f'D1S80.k{k}.example-path.png',1300,4000)
    lines=['# D1S80 six-k compact dBG results','',
           f'Run: `{run.name}`',f'STR commit: `{(run/"str-commit.txt").read_text().strip()}`','',
           f'Input: {inputs["bams"]} BAMs, {inputs["input_records"]} alignment records, {inputs["independent_molecules"]} included molecules.',
           'Observation scope: '+json.dumps(inputs['observation_scope']),
           'All paths, support, overlaps and SVG inventories passed independent validation. STR output regression passed.','',
           '| k | all S | all L | supported S | supported L | supported branch kmers | candidate edge retention | repeated candidate kmers |',
           '|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in comparison:
        lines.append(f'|{r["k"]}|{r["compact_nodes"]}|{r["compact_edges"]}|{r["supported_compact_nodes"]}|{r["supported_compact_edges"]}|{r["supported_branch_nodes"]}|{r["candidate_edge_retention"]:.3%}|{r["candidate_repeated_kmer_fraction"]:.3%}|')
    lines += ['', '## Interpretation limits','',
              'These metrics compare local sequence resolution and evidence retention. They do not establish genotyping accuracy.',
              'Low repeated-kmer fraction indicates less within-candidate sequence collapse, not a confirmed haplotype.',
              'Supported-view edges require three molecules within a sample; reference context is retained explicitly.',
              'Inspect both all-evidence and supported graphs; a simpler graph can result from lost read support.',
              'The 31/33/49 comparison is the primary candidate range; select the final viewing k after visual review.',
              '', '## Files','',
              '`results/dbg/D1S80/k<K>/D1S80.k<K>.gfa`: complete evidence graph.',
              'Paired `.svg`: all topology; `.paths.svg`: every P; `.summary.paths.svg`: reference and candidate paths.',
              '`paths-svg/`: one SVG per P. `.supported.*` and `paths-svg.supported/` provide the filtered view.',
              '`review/`: bounded PNG previews. Full sequence and molecule provenance remain in GFA/JSON/TSV.',
              '`results/dbg/D1S80/k_comparison.tsv`: all numerical comparison metrics.']
    (run/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    all_files=sorted(p for p in folder.rglob('*') if p.is_file())
    (run/'results-artifacts.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(run))+'\n' for p in all_files))
    with tarfile.open(run/f'{run.name}_results.tar.gz','w:gz') as archive:
        for path in [folder,review,run/'RESULTS.md',run/'results-artifacts.sha256',run/'results/acceptance',run/'results/preparation',run/'str-commit.txt',run/'source/SHA256SUMS']:
            archive.add(path,arcname=str(path.relative_to(run)))
    print('REPORT_PACKAGE_PASS')


if __name__=='__main__':
    main(Path(sys.argv[1]))
