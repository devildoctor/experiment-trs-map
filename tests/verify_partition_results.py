"""Reconstruct tokens against the original stored SAM SEQ and frozen model."""
from collections import Counter,defaultdict
from hashlib import sha256
from pathlib import Path
import json
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrpartition import load_model,check_partition,unit_rows,expression
from vntrgraph import Locus
from vntrreads import row_molecule_key
from vntrsummary import is_primary


def verify(root):
    audit=root/'audit'
    core=json.loads((root/'summary.json').read_text(encoding='utf-8'))
    payload=json.loads((audit/f'{core["locus"]}.structures.json').read_text(encoding='utf-8'))
    locus=Locus(**payload['locus'])
    model,dictionary=load_model(audit/'locus_model.json',locus)
    assert core['model_sha256']==sha256((audit/'locus_model.json').read_bytes()).hexdigest()
    assert [(m['id'],m['sequence']) for m in payload['motifs']]==[(m['id'],m['sequence']) for m in dictionary]
    assert payload['audit']['discovery_rounds_run']==0
    assert payload['audit']['dictionary_mode']=='frozen_apply'
    original={}
    with (audit/'records.jsonl').open(encoding='utf-8') as handle:
        for line in handle:
            record=json.loads(line);key=(record['row_id'],record['record_id'])
            assert key not in original
            original[key]=record
    stats=defaultdict(Counter);patterns=defaultdict(set);totals=Counter()
    for row in payload['rows']:
        for f in row['fragments']:
            raw=original[(row['id'],f['record_id'])]
            assert f['sequence']==raw['query_sequence'].upper()[f['qstart']:f['qend']]
            check_partition(f['sequence'],unit_rows(f['tokens']))
            if f['annotation']['alternative']:
                check_partition(f['sequence'],unit_rows(f['annotation']['alternative']))
            assert len(f['sequence'])==f['qend']-f['qstart']
            p=stats[row['platform']]
            p['fragments']+=1;p['bases']+=len(f['sequence'])
            p['ambiguous_fragments']+=bool(f['annotation']['annotation_ambiguity'])
            p['unknown_bases']+=sum(len(t['exact_bases']) for t in f['tokens'] if t['status']=='unknown')
            p['partial_bases']+=sum(len(t['exact_bases']) for t in f['tokens'] if t['status']=='partial')
            p['edited_units']+=sum(t['status']=='unit' and t['edit_distance']>0 for t in f['tokens'])
            totals.update((t['motif'],len(t['template']),len(t['exact_bases']),t['exact_bases'],t['status']) for t in f['tokens'])
            if row['platform']=='HiFi' and is_primary(f) and f['anchor']=='both':
                patterns[expression(unit_rows(f['tokens']),core['locus']=='ATXN2')].add(row_molecule_key(row))
    reference=payload['reference_structure']
    check_partition(reference['sequence'],unit_rows(reference['tokens']))
    if reference['annotation']['alternative']:
        check_partition(reference['sequence'],unit_rows(reference['annotation']['alternative']))
    result=dict(locus=core['locus'],model_version=model['version'],raw_records=len(original),
                fragments=sum(p['fragments'] for p in stats.values()),reconstruction_passed=True,
                frozen_dictionary_passed=True,platforms=dict(stats),
                reference_expression=expression(unit_rows(reference['tokens'])),
                reference_reverse_complement_expression=expression(unit_rows(reference['tokens']),True),
                reference_units=unit_rows(reference['tokens']),reference_annotation=reference['annotation'],
                top_hifi_both_anchor_patterns=[dict(expression=k,molecules=len(v)) for k,v in sorted(patterns.items(),key=lambda kv:(-len(kv[1]),kv[0]))[:10]],
                observed_unit_lengths=sorted({k[2] for k in totals if k[4]=='unit'}),
                mapping_context=core['mapping_context'])
    assert result['fragments']==core['fragments']
    (root/'verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in {'reference_units','top_hifi_both_anchor_patterns'}},ensure_ascii=False))


if __name__=='__main__':verify(Path(sys.argv[1]))
