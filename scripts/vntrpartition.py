"""Frozen locus models and a compact, lossless view of observed partitions."""
from collections import defaultdict
from pathlib import Path
import json
import re

from vntrcore import compact, revcomp
from vntrmotifs import definition
from vntrreads import row_molecule_key, write_table
from vntrsummary import is_primary


def load_model(path, locus):
    model = json.loads(Path(path).read_text(encoding='utf-8'))
    expected = dict(schema_version=1, locus=locus.name, assembly='GRCh38',
                    chrom=locus.chrom, start=locus.start, end=locus.end,
                    boundary_status=locus.boundary_status,
                    coordinate_system='0-based half-open', orientation='reference_forward',
                    phase_policy='fixed_template_phase', algorithm='seed_dp_v2')
    for key, value in expected.items():
        if model.get(key) != value:
            raise ValueError(f'Locus model mismatch: {key}; expected {value!r}')
    if not isinstance(model.get('version'), str) or not model['version']:
        raise ValueError('Locus model needs a nonempty version')
    dictionary = []
    ids, sequences = set(), set()
    for motif in model.get('motifs', []):
        name, sequence = motif.get('id', ''), motif.get('sequence', '')
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', name) or name in ids:
            raise ValueError('Invalid or duplicate model motif ID')
        if not re.fullmatch(r'[ACGT]+', sequence) or sequence in sequences:
            raise ValueError('Invalid or duplicate model motif sequence')
        if motif.get('length') != len(sequence):
            raise ValueError(f'Motif length mismatch: {name}')
        ids.add(name); sequences.add(sequence)
        # Explicit units are frozen; no sample-driven renaming or demotion.
        record = dict(id=name, base=name, sequence=sequence, source='frozen_model',
                      support=0, **definition(sequence))
        record.update(level='base', match_orientations=['+'])
        dictionary.append(record)
    if not dictionary:
        raise ValueError('Locus model needs at least one motif')
    return model, dictionary


def unit_rows(tokens):
    return [dict(order=i+1, start=t['start'], end=t['end'], motif=t['motif'],
                 template=t['template'], template_length=len(t['template']),
                 observed_length=len(t['exact_bases']), bases=t['exact_bases'],
                 status=t['status'], edits=t['edit_distance'],
                 ambiguity=t.get('annotation_ambiguity', False))
            for i, t in enumerate(tokens)]


def check_partition(sequence, units):
    position = 0
    for u in units:
        if u['start'] != position or u['end']-u['start'] != u['observed_length']:
            raise ValueError('Partition has a gap, overlap or length mismatch')
        if sequence[u['start']:u['end']] != u['bases']:
            raise ValueError('Partition changed observed bases')
        position = u['end']
    if position != len(sequence) or ''.join(u['bases'] for u in units) != sequence:
        raise ValueError('Partition does not reconstruct the input')


def expression(units, reverse=False):
    blocks = []
    for u in reversed(units) if reverse else units:
        bases = revcomp(u['bases']) if reverse else u['bases']
        suffix = '?' if u['status']=='unknown' else '~' if u['status']=='partial' else '*' if u['edits'] else ''
        key = (bases, suffix)
        if blocks and blocks[-1][0] == key and u['status']=='unit':
            blocks[-1][1] += 1
        else: blocks.append([key, 1])
    return ' – '.join(f'({b}){s}×{n}' if n>1 else b+s for (b,s),n in blocks)


def export_core(destination, locus, rows, dictionary, reference, audit):
    motifs = [dict(id=m['id'], length=len(m['sequence']), bases=m['sequence'],
                   reverse_complement=revcomp(m['sequence'])) for m in dictionary]
    observations = []
    ref = dict(id='REFERENCE', sample='参考', platform='参考', coverage='catalog_interval',
               anchor='reference', alignment='reference', mapq=None,
               ambiguity=reference.get('annotation',{}).get('annotation_ambiguity',False),
               score_gap=reference.get('annotation',{}).get('score_gap'),
               units=unit_rows(reference['tokens']),
               alternatives=unit_rows(reference.get('annotation',{}).get('alternative',[])))
    check_partition(reference['sequence'], ref['units'])
    if ref['alternatives']:check_partition(reference['sequence'], ref['alternatives'])
    observations.append(ref)
    variants = defaultdict(lambda: dict(molecules=set(), hifi=set()))
    for row in sorted(rows, key=lambda r:(r['sample'],r['platform'],r['id'])):
        for f in row['fragments']:
            units = unit_rows(f['tokens'])
            check_partition(f['sequence'], units)
            alternative = unit_rows(f.get('annotation',{}).get('alternative',[]))
            if alternative: check_partition(f['sequence'], alternative)
            observations.append(dict(id=f['id'], sample=row['sample'], platform=row['platform'],
                                     coverage=f['coverage'], anchor=f['anchor'],
                                     alignment=f['alignment_type'], mapq=f['mapq'],
                                     flags=f.get('flags',[]), qstart=f['qstart'], qend=f['qend'],
                                     ambiguity=f.get('annotation',{}).get('annotation_ambiguity',False),
                                     score_gap=f.get('annotation',{}).get('score_gap'),
                                     units=units, alternatives=alternative))
            for u in units:
                key=(u['motif'],u['template'],u['status'],u['bases'])
                variants[key]['molecules'].add(row_molecule_key(row))
                if row['platform']=='HiFi' and is_primary(f):
                    variants[key]['hifi'].add(row_molecule_key(row))
    actual = [dict(motif=k[0],template=k[1],status=k[2],template_length=len(k[1]),
                   observed_length=len(k[3]),bases=k[3],molecules=len(v['molecules']),
                   hifi_primary_molecules=len(v['hifi']))
              for k,v in sorted(variants.items(),key=lambda kv:(-len(kv[1]['molecules']),kv[0]))]
    metadata = dict(locus=locus.name,chrom=locus.chrom,start=locus.start,end=locus.end,
                    boundary_status=locus.boundary_status,coordinate_system='0-based half-open',
                    orientation='reference_forward', model_version=audit.get('model_version'),
                    model_sha256=audit.get('model_sha256'),dictionary_mode=audit['dictionary_mode'],
                    bam_count=audit['bam_count'],read_rows=len(rows),fragments=len(observations)-1,
                    molecules=len({row_molecule_key(r) for r in rows if r['fragments']}),
                    reconstruction_checks=len(observations),
                    reference_expression=expression(ref['units']),reference_ambiguity=ref['ambiguity'],
                    reference_reverse_complement_expression=expression(ref['units'],True),
                    mapping_context=audit['mapping_context'],motifs=motifs)
    write_table(destination/'motifs.tsv',motifs,['id','length','bases','reverse_complement'])
    write_table(destination/'observed_units.tsv',actual,
                ['motif','template','template_length','observed_length','bases','status','molecules','hifi_primary_molecules'])
    write_table(destination/'partitions.tsv',
                [dict(observation=o['id'],sample=o['sample'],platform=o['platform'],**u)
                 for o in observations for u in o['units']],
                ['observation','sample','platform','order','start','end','motif','template',
                 'template_length','observed_length','bases','status','edits','ambiguity'])
    (destination/'summary.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
    payload=dict(summary=metadata, observations=observations, observed_units=actual)
    serialized=compact(payload).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
    html=Path(__file__).with_name('partition_view.html').read_text(encoding='utf-8')
    (destination/'index.html').write_text(html.replace('__DATA__',serialized).replace('__LOCUS__',locus.name),encoding='utf-8')
    return metadata
