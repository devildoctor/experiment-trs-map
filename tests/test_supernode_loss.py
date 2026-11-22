from collections import Counter
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrsupernode_loss import plan_substitutions
from vntrsupernode import CERTAIN


def run_case(strings,roles=None,motifs=None,fraction=.03,protected=False,certain=CERTAIN):
    # All strings share one slot. The last string normally has the least support.
    roles=roles or ['observation']*len(strings)
    motifs=motifs or ['A']*len(strings)
    records=[];tokens={};slots={};visits={}
    for i,s in enumerate(strings):
        p=f'P{i+1:04d}'
        records.append(dict(id=p,role=roles[i],sequence=s,molecules=[f'm{i}.{k}' for k in range(len(strings)-i)]))
        tokens[p]=[dict(start=0,end=len(s),motif=motifs[i],exact_bases=s,status='unit')]
        slots[p]=[(0,1,0)]
        visits[p]=[dict(supernode_id='S',token_interval=[0,1],structure_assessment=certain,
                        boundary_assessment=CERTAIN,position_assessment='aligned')]
    display=dict(nodes={'S':dict(id='S',protected=protected)})
    return plan_substitutions(records,tokens,slots,display,visits,dict(loss_fraction=fraction,max_unit_substitutions=2))


class LossTests(unittest.TestCase):
    def test_floored_budgets_s_only_and_frozen_representative(self):
        for n,budget in [(20,0),(33,0),(34,1),(66,1),(67,2)]:
            a='A'*n
            for cost in (1,2,3):
                b='C'*cost+a[cost:]
                rep,events,ds=run_case([a,b])
                self.assertEqual(sum(map(len,events.values())),cost if cost<=budget else 0)
                for e in events.get('P0002',[]):
                    self.assertEqual((e['op'],e['representative_occurrence']),('S','P0001:u0'))
        a='A'*100
        self.assertFalse(run_case([a,'C'+a[1:]],fraction=0)[1])
        self.assertFalse(run_case(['ACGT'*25,'CGTA'*25])[1])

    def test_protection_types_reference_and_no_chain(self):
        a='A'*100;b='CC'+a[2:];c='CCCC'+a[4:]
        self.assertFalse(run_case([a,b],motifs=['A','B'])[1])
        self.assertFalse(run_case([a,b],protected=True)[1])
        for state in ('not_assessed','ambiguous'):
            self.assertFalse(run_case([a,b],certain=state)[1])
        self.assertFalse(run_case([a,b],roles=['reference','observation'])[1])
        rep,events,_=run_case([a,b,c])
        self.assertEqual(set(events),{'P0002'})
        self.assertEqual(rep['P0003'][0]['exact_bases'],c)
        # A retained REF still uses the candidate source node. No N+E gain.
        self.assertFalse(run_case([a,b,b],roles=['observation','observation','reference'])[1])
