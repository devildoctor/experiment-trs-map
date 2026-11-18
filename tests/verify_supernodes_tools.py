"""Independent ODGI FASTA and Bandage import/render checks on an existing bundle."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET


def command(*args,env=None):
    result=subprocess.run(args,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
    return result.stdout.decode()


def fasta(text):
    seqs={}
    for line in text.splitlines():
        if line.startswith('>'):
            name=line[1:].split()[0]
            if name in seqs: raise ValueError('Duplicate FASTA record '+name)
            seqs[name]=''
        elif line: seqs[name]+=line
    return seqs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for flag in ('bundle','out','odgi','bandage'): parser.add_argument('--'+flag,type=Path,required=True)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    bandenv=dict(os.environ,QT_QPA_PLATFORM='offscreen')
    if os.getenv('BANDAGE_LIBSTDCXX'):bandenv['LD_PRELOAD']=os.environ['BANDAGE_LIBSTDCXX']
    report=dict(odgi=command(args.odgi,'version').strip(),
                bandage=command(args.bandage,'--version',env=bandenv).strip(),checks={})
    synthetic=args.out/'repeated_iupac.gfa'
    synthetic.write_text('H\tVN:Z:1.0\nS\t1\tANRY\nL\t1\t+\t1\t+\t0M\nP\tRepeated\t1+,1+,1+\t*\n')
    locus=json.loads((args.bundle/'summary.json').read_bytes())['locus']
    inputs=[(m,args.bundle/(locus+'.'+m+'.gfa'),fasta((args.bundle/('paths.'+m+'.fa')).read_text()))
            for m in ('original','represented')]
    inputs.append(('repeated_iupac',synthetic,{'Repeated':'ANRY'*3}))
    for name,gfa,expected in inputs:
        graph=args.out/(name+'.og')
        result={}
        report['checks'][name]=result
        try:
            command(args.odgi,'build','-g',gfa,'-o',graph,'-t','1')
            text=command(args.odgi,'paths','-i',graph,'-f')
            (args.out/(name+'.odgi.fa')).write_text(text)
            actual=fasta(text)
            result.update(odgi_paths=len(actual),odgi_sequence_match=actual==expected)
            if actual!=expected:result['mismatch_paths']=sorted(n for n in set(actual)|set(expected) if actual.get(n)!=expected.get(n))
        except subprocess.CalledProcessError as e:
            result.update(odgi_sequence_match=False,odgi_error=e.stderr.decode(errors='replace'))
        info=command(args.bandage,'info',gfa,env=bandenv)
        (args.out/(name+'.bandage.txt')).write_text(info)
        svg=args.out/(name+'.bandage.svg')
        command(args.bandage,'image',gfa,svg,'--height','900',env=bandenv)
        result.update(bandage_info=bool(info.strip()),bandage_svg=ET.parse(svg).getroot().tag.endswith('svg'),
                      svg_bytes=svg.stat().st_size,gfa_sha256=hashlib.sha256(gfa.read_bytes()).hexdigest())
    (args.out/'tool_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    for m in ('original','represented'):
        if not all(report['checks'][m][k] for k in ('odgi_sequence_match','bandage_info','bandage_svg')):
            raise ValueError('Tool compatibility failed for '+m)


if __name__=='__main__':main()
