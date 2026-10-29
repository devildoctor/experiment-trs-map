"""Independent output checks, optional exact read-inventory comparison to r3."""
import json
from pathlib import Path
from collections import defaultdict,Counter
import sys


def main():
    out=Path(sys.argv[1]);locus=out.name
    d=json.loads((out/f'{locus}.structures.json').read_text(encoding='utf-8'))
    nodes={};paths={};support={};nodekeys=defaultdict(set);abstractkeys=defaultdict(set)
    for line in (out/f'{locus}.reads.gfa').read_text(encoding='utf-8').splitlines():
        f=line.split('\t')
        if f[0]=='S':nodes[f[1]]=f[2];support[f[1]]=int(next(x[5:] for x in f if x.startswith('RC:i:')))
        elif f[0]=='P':paths[f[1]]=f[2].split(',')
    checks=0;rows=d['rows'];seen=set()
    for r in rows:
        row_key=(r['sample'],r['platform'],r['bam'],r.get('namespace',''),r['read_id']);assert row_key not in seen;seen.add(row_key)
        key=(r['sample'],r.get('namespace',''),r['read_id'])
        if 'molecule_key' in r:assert tuple(json.loads(r['molecule_key']))==key
        for f in r['fragments']:
            ts=f['tokens'];assert ''.join(nodes[x[:-1]] for x in paths[f['id']])==f['sequence']
            end=0
            for t in ts:
                assert t['start']==end;end=t['end']
                assert f['sequence'][t['start']:t['end']]==t['exact_bases']
                assert t['qstart']==f['qstart']+t['start'] and t['qend']==f['qstart']+t['end']
                nodekeys[t['exact_node']].add(key);abstractkeys[t['structure_node']].add(key)
            assert end==len(f['sequence']);checks+=1
    assert checks==len(paths)==d['audit']['sequence_reconstruction_checks_passed']
    assert all(support[n]==len(keys) for n,keys in nodekeys.items())
    assert all(n['support']==len(abstractkeys[n['id']]) for n in d['graph']['structure_nodes'])
    report=dict(status='passed',reads=len(rows),fragments=checks,exact_gfa_reconstruction=True,support_dedup=True)
    if len(sys.argv)>2:
        previous=Path(sys.argv[2]);old=json.loads((previous/f'{locus}.reads.json').read_text(encoding='utf-8'))
        # 旧版可能把同 QNAME 的不同 LB 合并为一行；比较实际记录和片段多重集合。
        def inventory(rs):
            result={}
            for r in rs:
                key=(r['sample'],r['platform'],r['bam'],r['read_id'])
                item=result.setdefault(key,dict(records=0,fragments=Counter()))
                item['records']+=r['record_count']
                item['fragments'].update((f['sequence'],f['qstart'],f['qend']) for f in r['fragments'])
            return result
        before=inventory(old['rows']);after=inventory(rows);assert before.keys()==after.keys(),'Read inventory changed'
        for k,r in before.items():
            rr=after[k];assert r['records']==rr['records']
            assert r['fragments']==rr['fragments']
        report.update(r3_inventory_identical=True,r3_extracted_bases_identical=True)
    (out/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('OUTPUT_VERIFICATION_PASS',json.dumps(report),flush=True)


if __name__=='__main__':main()
