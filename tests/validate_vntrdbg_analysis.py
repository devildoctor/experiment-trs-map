"""Independent invariants for a completed D1S80 joint-support run."""
import csv
import json
from pathlib import Path
import sys


def rows(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main(out, input_file):
    records = json.loads(input_file.read_text(encoding='utf-8'))
    candidates = {record['id'].removeprefix('CAND_'): record for record in records
                  if record['kind'] == 'candidate'}
    queries = rows(out/'joint_queries.tsv')
    support = {row['query_id']: row for row in rows(out/'joint_support.tsv')}
    evidence = rows(out/'joint_molecules.tsv')
    units = rows(out/'primitive_units.tsv')
    assert len(queries) == len(support)
    assert len({row['query_id'] for row in queries}) == len(queries)
    assert all(row['exact_bases'] == candidates[row['candidate_id']]['sequence'][
               int(row['start0']):int(row['end0'])] for row in queries)
    by_query = {}
    for row in evidence:
        by_query.setdefault(row['query_id'], []).append(row)
        hits = json.loads(row['occurrences_json'])
        assert (row['decision'] == 'exact_support') == bool(hits)
        assert len(hits) == int(row['occurrence_count'])
        for hit in hits:
            source = next(record for record in records if record['id'] == hit['path_id'])
            query = next(q for q in queries if q['query_id'] == row['query_id'])
            assert source['kind'] == 'observed'
            assert source['molecule_key'] == row['molecule_key']
            assert source['sequence'][hit['start0']:hit['end0']] == query['exact_bases']
            if query['kind'] == 'full_candidate_exact':
                assert source['coverage'] == 'two_flank_anchored'
                assert source['sequence'] == query['exact_bases']
    for query in queries:
        s = support[query['query_id']]
        e = by_query[query['query_id']]
        assert len({row['molecule_key'] for row in e}) == len(e)
        assert int(s['tested_molecules']) == len(e)
        assert int(s['exact_support_molecules']) == sum(row['decision']=='exact_support' for row in e)
        assert int(s['uninformative_molecules']) + int(s['exact_support_molecules']) == len(e)
        assert s['conflict_molecules'] == s['support_fraction'] == 'NA'
    for name, candidate in candidates.items():
        for scheme in ('primary', 'alternative'):
            grouped = sorted((row for row in units if row['candidate_id']==name
                              and row['decomposition']==scheme), key=lambda row:int(row['unit_index']))
            if not grouped:
                continue
            assert ''.join(row['exact_bases'] for row in grouped) == candidate['sequence']
            assert [int(row['unit_index']) for row in grouped] == list(range(1,len(grouped)+1))
            assert grouped[0]['start0'] == '0'
            assert int(grouped[-1]['end0']) == len(candidate['sequence'])
    assert len(rows(out/'unit_summary.tsv')) == len(candidates)
    print(json.dumps(dict(status='PASS', queries=len(queries), molecule_query_rows=len(evidence),
                          unit_rows=len(units), candidates=len(candidates))))


if __name__ == '__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]))
