"""Compare v2 display budgets on an existing verified B2 source_members file."""
import argparse
from pathlib import Path
from vntrgraph import Locus
from vntrdisplay_export import run_comparison


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--model',type=Path,required=True);p.add_argument('--policy',type=Path,required=True)
    p.add_argument('--dot',default='dot');args=p.parse_args()
    import json
    m=json.loads(args.model.read_text(encoding='utf-8'))
    locus=Locus(m['locus'],m['chrom'],m['start'],m['end'],(),boundary_status=m['boundary_status'])
    run_comparison(args.input,args.model,args.policy,locus,args.out,args.dot)


if __name__=='__main__':main()
