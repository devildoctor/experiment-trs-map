#!/usr/bin/env python3
"""VNTR point graphs: exact sequences, molecule evidence, candidate uncertainty."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from difflib import SequenceMatcher
from hashlib import sha256
from html import escape
from pathlib import Path
import csv
import json
import re
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vntrcore import (VERSION, DNA, Read, Allele, Token, compact, deduplicate_reads,
    edit_distance, similarity, align_to_seed, make_consensus, cluster_reads,
    call_candidates, discover_motifs, motif_dictionary, decompose, token_counts,
    build_graph, graph_evidence, dbg_evidence)
from vntrcore import revcomp as reverse_complement
import vntranchors
from vntrmotifs import compound_annotations

COLORS = ("#4477AA", "#EE7733", "#228833", "#AA3377", "#66CCEE", "#CCBB44")

@dataclass(frozen=True)
class Locus:
    name: str
    chrom: str
    start: int
    end: int
    motifs: tuple[str, ...]
    search_start: int | None = None
    search_end: int | None = None
    boundary_status: str = "candidate_window"

class IndexedFastaReference:
    """独立 FASTA/FAI 读取器，所有坐标均为 0-based 半开区间。"""

    def __init__(self, path: Path):
        self.path = path
        fai = Path(str(path) + ".fai")
        if not path.is_file() or not fai.is_file():
            raise ValueError("参考 FASTA 或同名 FAI 不存在")
        self.entries = {}
        for line in fai.read_text(encoding="utf-8").splitlines():
            name, length, offset, line_bases, line_width, *_ = line.split("\t")
            self.entries[name] = tuple(map(int, (length, offset, line_bases, line_width)))
        self.handle = path.open("rb")

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.handle.close()

    def contig_length(self, chrom: str) -> int:
        return self.entries[chrom][0]

    def fetch(self, chrom: str, start: int, end: int) -> str:
        if chrom not in self.entries:
            raise ValueError(f"参考基因组无染色体 {chrom}")
        length, offset, line_bases, line_width = self.entries[chrom]
        if start < 0 or end > length or end < start:
            raise ValueError(f"参考坐标越界: {chrom}:{start}-{end}")
        parts = []
        while start < end:
            row, column = divmod(start, line_bases)
            size = min(end - start, line_bases - column)
            self.handle.seek(offset + row * line_width + column)
            part = self.handle.read(size)
            if len(part) != size:
                raise ValueError("FASTA/FAI 不匹配")
            parts.append(part)
            start += size
        return b"".join(parts).decode("ascii").upper()


def read_catalog(path):
    loci, seen = [], set()
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise ValueError(f"BED row {number}: expected at least four columns")
        chrom, a, b, meta = fields[:4]
        search_start, search_end = int(a), int(b)
        info = dict(x.split("=",1) for x in meta.split(";") if "=" in x)
        name = info.get("ID",meta)
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", name) or name in seen:
            raise ValueError(f"Invalid or duplicate locus ID: {name}")
        if ("REPEAT_START" in info) != ("REPEAT_END" in info):
            raise ValueError("REPEAT_START and REPEAT_END must be supplied together")
        start = int(info.get("REPEAT_START",search_start))
        end = int(info.get("REPEAT_END",search_end))
        if not 0 <= search_start <= start < end <= search_end:
            raise ValueError(f"Invalid search/repeat coordinates: {name}")
        status = info.get("BOUNDARY","candidate_window")
        if status not in {"confirmed","candidate_window"}:
            raise ValueError(f"Unsupported BOUNDARY: {status}")
        motifs = tuple(x.upper() for x in info.get("MOTIFS","").split(",") if x)
        if any(set(x)-DNA for x in motifs):
            raise ValueError(f"Non-ACGT motif: {name}")
        loci.append(Locus(name,chrom,start,end,motifs,search_start,search_end,status))
        seen.add(name)
    if not loci:
        raise ValueError("Empty locus catalog")
    return loci


def discover_bams(folder: Path, locus: str) -> list[tuple[str, str, Path]]:
    """Read locus/sample/platform or legacy sample/platform BAM layouts."""
    found = []
    root = folder/locus if (folder/locus).is_dir() else folder
    for bam in sorted(root.glob("*/*/*.bam")):
        if bam.name.rsplit(".", 2)[-2] != locus:
            continue
        if not (Path(str(bam) + ".bai").exists() or bam.with_suffix(".bai").exists()):
            raise ValueError(f"缺少 BAM 索引: {bam}")
        sample, platform = bam.parent.parent.name, bam.parent.name
        if platform not in {"HiFi", "ONT", "ONTUL"}:
            raise ValueError(f"未知测序平台: {bam}")
        found.append((sample, platform, bam))
    return found


def reference_to_query(pos0: int, cigar: str, target: int) -> int | None:
    """CIGAR 投影仅用于侧翼搜索的初始范围，不决定重复边界。"""
    ref, query = pos0, 0
    for size_text, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        size = int(size_text)
        if op in "M=X":
            if ref <= target <= ref + size:
                return query + target - ref
            ref += size
            query += size
        elif op in "DN":
            if ref <= target <= ref + size:
                return query
            ref += size
        elif op in "IS":
            query += size
    return None


def flank_hit(flank: str, sequence: str, expected: int | None,
              minimum_identity: float) -> tuple[int, int, float] | None:
    """在 read 中独立寻找一个高分且位置唯一的侧翼。"""
    length = len(flank)
    starts = set()
    begin = sequence.find(flank)
    while begin >= 0:
        starts.add(begin)
        begin = sequence.find(flank, begin + 1)
    if starts:
        if len(starts) != 1:
            return None
        start = next(iter(starts))
        return start,start+length,1.0
    if not starts:
        # 先找短精确种子，再对整段侧翼评分；允许 ONT 的局部 indel。
        for offset in (0, max(0, length // 2 - 10), max(0, length - 20)):
            seed = flank[offset:offset + 20]
            begin = sequence.find(seed)
            while begin >= 0:
                starts.update(range(max(0, begin - offset - 5), min(len(sequence), begin - offset + 6)))
                begin = sequence.find(seed, begin + 1)
    # CIGAR is only a hint: a large insertion at a boundary can put the true
    # anchor far from the projected coordinate. Do not discard unique anchors.
    if not starts:
        return None
    scored = []
    for start in starts:
        for delta in (-8, -4, 0, 4, 8):
            end = start + length + delta
            if end > len(sequence) or end <= start:
                continue
            segment = sequence[start:end]
            score = SequenceMatcher(None, flank, segment, autojunk=False).ratio()
            if score >= minimum_identity:
                scored.append((score, start, end))
    if not scored:
        return None
    scored.sort(reverse=True)
    best = scored[0]
    if any(item[0] == best[0] and item[1:] != best[1:] for item in scored[1:]):
        return None
    # 第二个独立位置得分相近时，不能唯一定位边界。
    if any(abs(item[1] - best[1]) > 10 and item[0] >= best[0] - .02 for item in scored[1:]):
        return None
    return best[1], best[2], best[0]


def extract_sam_records(lines, sample, platform, locus, reference, flank_length,
                        min_mapq, source="", namespace=None, profiles=None):
    """Compatibility entry point with the same audited extraction rules."""
    if profiles is None:
        if not hasattr(reference, 'path'):
            raise ValueError('Reference profiles required for in-memory references')
        profiles, _ = vntranchors.reference_profiles(reference.path,[locus],reference,[flank_length])
    reads,audit,_ = vntranchors.extract_records(lines,sample,platform,locus,reference,
                                               profiles,min_mapq,source,namespace)
    return reads,audit

def write_tsv(path, fields, rows):
    with path.open("w",encoding="utf-8",newline="") as handle:
        writer = csv.DictWriter(handle,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k:compact(v) if isinstance(v,(dict,list,tuple)) else v for k,v in row.items()})


def source_hashes():
    return {p.name:sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('vntr*.py')}


def write_gfa(path,nodes,edges,paths):
    evidence = graph_evidence(nodes,edges,paths)
    node_counts = {x["id"]:x["read_support"] for x in evidence["nodes"]}
    edge_counts = {(x["start"],x["end"]):x["read_support"] for x in evidence["edges"]}
    with path.open("w",encoding="utf-8") as out:
        out.write("H\tVN:Z:1.0\tPG:Z:vntrgraph-2.0\n")
        for (_,family,exact),node in nodes.items():
            out.write(f"S\t{node}\t{exact}\tFT:Z:{family}\tRC:i:{node_counts[node]}\n")
        for a,b in sorted(edges):
            out.write(f"L\t{a}\t+\t{b}\t+\t0M\tRC:i:{edge_counts[a,b]}\n")
        for allele,ids in paths:
            out.write(f"P\t{allele.name}\t{','.join(x+'+' for x in ids)}\t*"
                      f"\tRC:i:{len(allele.read_ids)}\tPS:Z:{compact(allele.support)}"
                      f"\tAC:i:{len(allele.ambiguous_ids)}\tST:Z:{allele.call_status}\n")
    path.with_suffix(".json").write_text(compact(evidence),encoding="utf-8")


def compress_tokens(tokens):
    runs = []
    for token in tokens:
        key = (token.family,token.exact_bases,token.status)
        if runs and runs[-1][0] == key:
            runs[-1][1] += 1
        else:
            runs.append([key,1])
    assert ''.join(key[1]*count for key,count in runs) == ''.join(t.exact_bases for t in tokens)
    return runs


def write_svg(path,locus,alleles,reason=""):
    variants = list(dict.fromkeys((t.family,t.exact_bases,t.status) for a in alleles for t in a.tokens))
    labels = {key:f"U{i+1}" for i,key in enumerate(variants)}
    families = sorted({key[0] for key in variants})
    rows = [(a,compress_tokens(a.tokens)) for a in alleles]
    width = max(1150,max((120+len(runs)*118 for _,runs in rows),default=0))
    legend_height = sum(44+16*((len(k[1])+79)//80) for k in variants)
    height = max(220,150+len(rows)*155+legend_height)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="28" y="34" font-family="sans-serif" font-size="24">{escape(locus)} VNTR</text>',
           '<text x="28" y="58" font-family="sans-serif" font-size="13">Saved paths only; xN counts displayed sequence units, not a definitive copy-number call. Labels are locus-local.</text>']
    if reason:
        svg.append(f'<text x="28" y="104" font-family="sans-serif" font-size="15" fill="#994444">{escape(reason)}</text>')
    for index,(allele,runs) in enumerate(rows):
        y = 94+index*155
        title = f"{allele.name} | {allele.call_status} | {len(allele.sequence)} bp"
        svg.append(f'<text x="28" y="{y}" font-family="sans-serif" font-size="15">{escape(title)}</text>')
        for i,(key,count) in enumerate(runs):
            x = 50+i*118
            family,exact,status = key
            color = "#888888" if status == "unknown" else COLORS[families.index(family)%len(COLORS)]
            label = labels[key] + (f" ×{count}" if count > 1 else "") + ("*" if status == "partial" else "")
            svg.append(f'<rect x="{x}" y="{y+12}" width="100" height="30" rx="4" fill="{color}"/>')
            svg.append(f'<text x="{x+8}" y="{y+33}" font-family="monospace" font-size="13" fill="white" data-exact="{exact}" data-count="{count}" data-path="{escape(allele.name)}">{escape(label)}</text>')
            if i+1 < len(runs):
                svg.append(f'<path d="M{x+100},{y+27} h16 l-4,-3 m4,3 l-4,3" fill="none" stroke="#555"/>')
        counts = token_counts(allele.tokens)
        support = f"molecules={len(allele.read_ids)}; platforms={compact(allele.support)}; ambiguous={len(allele.ambiguous_ids)}"
        note = "observed units=" + compact(counts) if not allele.uncertain else "UNCERTAIN; no definitive copy number. " + ','.join(allele.flags)
        svg.append(f'<text x="28" y="{y+62}" font-family="sans-serif" font-size="12">{escape(support)}</text>')
        svg.append(f'<text x="28" y="{y+80}" font-family="sans-serif" font-size="12">{escape(note[:170])}</text>')
        annotation = allele.annotation
        if annotation:
            alternative = annotation.get('alternative_metrics') or {}
            info = (f"Primary full units={annotation['metrics']['full_units']}; alternative full units={alternative.get('full_units','none')}; "
                    f"cost gap={annotation.get('score_gap')}; fixed template phase; see annotation_alternatives.json")
            svg.append(f'<text x="28" y="{y+98}" font-family="sans-serif" font-size="12">{escape(info)}</text>')
        block_note = f"Compound observations={sum(h['count'] for h in allele.compound_details if h['basis']=='observed_base_token_pattern')}; catalog block counts={compact(dict(Counter(h['catalog_motif'] for h in allele.compound_details if 'catalog_motif' in h)))} (separate layers; do not sum)"
        svg.append(f'<text x="28" y="{y+116}" font-family="sans-serif" font-size="12">{escape(block_note[:170])}</text>')
    y = 120+len(rows)*155
    for key in variants:
        family,exact,status = key
        svg.append(f'<text x="28" y="{y}" font-family="sans-serif" font-size="13">{labels[key]} | {len(exact)} bp | {status} | exact bases:</text>')
        y += 18
        for start in range(0,len(exact),80):
            svg.append(f'<text x="45" y="{y}" font-family="monospace" font-size="12">{exact[start:start+80]}</text>')
            y += 16
        y += 26
    svg.append('</svg>')
    path.write_text('\n'.join(svg),encoding="utf-8")


def process_locus(locus,args,reference,profiles):
    destination = args.out/locus.name
    destination.mkdir(parents=True,exist_ok=True)
    parameters = dict(version=VERSION,locus=asdict(locus),coordinate_system="0-based_half-open",
                      boundary_shift=0,score_definition="descriptive_edit_cost_not_probability",
                      support_definition="unique_molecule_keys; platform subtotals may overlap",
                      source_hashes=source_hashes(),
                      catalog_sha256=sha256(args.catalog.read_bytes()).hexdigest(),
                      reference_equivalence_sha256=sha256(args.reference_equivalence.read_bytes()).hexdigest() if args.reference_equivalence else None,
                      reference_fai_sha256=sha256(Path(str(args.reference)+'.fai').read_bytes()).hexdigest(),
                      reference_size=args.reference.stat().st_size,
                      reference_mtime_ns=args.reference.stat().st_mtime_ns,
                      anchor_algorithm='graded_multilength_edit_alignment_v1',
                      bam_scope='all_records_in_supplied_locus_BAM; no_external_fetch',
                      input_flank_limit=500,anchor_lengths=[n for n in vntranchors.LENGTHS if n<=args.anchor_length],
                      identity_thresholds=vntranchors.IDENTITY,observed_fraction=vntranchors.OBSERVED_FRACTION,
                      localization_fraction=vntranchors.LOCALIZED_FRACTION,
                      localization_tolerance=vntranchors.POSITION_TOLERANCE,
                      max_anchor_locations=vntranchors.MAX_CANDIDATES,
                      reference_flank_grades={s:[(t['length'],t['grade']) for t in profiles[locus.name+':'+s]] for s in ('L','R')})
    parameters.update({k:str(v) for k,v in vars(args).items() if k != 'source_namespaces'})
    write_tsv(destination/'parameters.tsv',['parameter','value'],[dict(parameter=k,value=v) for k,v in parameters.items()])
    print(f"{locus.name}: extract anchored molecules",flush=True)
    bams = discover_bams(args.bam_dir,locus.name)
    reads,per_bam,details,failures = vntranchors.extract_bams(bams,locus,reference,args.samtools,
                                                           profiles,args.min_mapq,args.source_namespaces)
    if not bams:
        failures.append('No matching BAM files')
    with (destination/'anchor_decisions.jsonl').open('w',encoding='utf-8') as decisions:
        for detail in details:
            decisions.write(compact(detail)+'\n')
    reads,dedup = deduplicate_reads(reads)
    write_tsv(destination/'molecules.tsv',list(asdict(reads[0]))+['molecule_key'] if reads else ['molecule_key'],
              [dict(asdict(r),molecule_key=r.molecule_key) for r in reads])
    write_tsv(destination/'dedup.tsv',['molecule_key','status','observations','records','platforms','sequence_hashes'],dedup)
    print(f"{locus.name}: {len(reads)} molecules; candidate discovery and joint assignment",flush=True)
    candidates,assignments,genotypes = call_candidates(reads,args.min_support,args.ploidy,args.score_margin)
    write_tsv(destination/'assignments.tsv',['molecule_key','sample','platform','platforms','status','candidate','best','score_gap','near_candidates','scores'],assignments)
    write_tsv(destination/'genotypes.tsv',['sample','ploidy','status','supported_candidates','assignments','combinations','combination_status','phase'],genotypes)
    alleles = [a for a in candidates if len(a.read_ids) >= args.min_support]
    motifs = motif_dictionary(reference.fetch(locus.chrom,locus.start,locus.end),alleles,reads,locus.motifs,args.min_support)
    write_tsv(destination/'motifs.tsv',['motif','family','source','level','primitive_sequence','repeat_factor','phase_policy','usable','validated','evidence_status','read_support','substring_read_support','molecule_ids','support_positions'],motifs)
    base = tuple(m['motif'] for m in motifs if m['usable'] and m['level'] == 'base')
    compound = tuple(m['motif'] for m in motifs if m['source']=='catalog' and m['level']=='compound')
    supported_families = {m['family'] for m in motifs if m['validated']}
    print(f"{locus.name}: {len(alleles)} supported sequences; motif annotation and graphs",flush=True)
    for allele in alleles:
        allele.tokens,explained,allele.annotation = decompose(allele.sequence,base,return_details=True)
        allele.compound_details = compound_annotations(allele.tokens,compound)
        allele.dbg_k,allele.graph_support,allele.unitigs = dbg_evidence(allele,reads,destination/f'{allele.name}.dbg.json')
        if locus.boundary_status != 'confirmed':
            allele.flags.append('candidate_boundaries')
        if any(r.molecule_key in allele.read_ids and r.anchor_scores and min(r.anchor_scores)<1 for r in reads):
            allele.flags.append('approximate_anchor')
        if any(t.status == 'unknown' for t in allele.tokens):
            allele.flags.append('unknown_bases')
        if any(t.status != 'unknown' and t.family not in supported_families for t in allele.tokens):
            allele.flags.append('unsupported_motif_seed')
        if any(t.status == 'partial' for t in allele.tokens):
            allele.flags.append('partial_units')
        if any(t.annotation_ambiguity for t in allele.tokens):
            allele.flags.append('annotation_ambiguity')
        if allele.graph_support < .90:
            allele.flags.append('low_local_graph_support')
        allele.uncertain |= any(x not in {'stable_residual_split','residual_discovery'} for x in allele.flags)
    reason = '; '.join(failures)
    if not alleles and not reason:
        reason = ('Two nonoverlapping unique 31-mers not established in both reference flanks' if not all(any(t['grade']=='preferred' for t in profiles[locus.name+s]) for s in (':L',':R'))
                  else 'Insufficient independent spanning molecules; no full-length allele inferred')
    if alleles:
        left = reference.fetch(locus.chrom,locus.start-25,locus.start)
        right = reference.fetch(locus.chrom,locus.end,locus.end+25)
        nodes,edges,paths = build_graph(alleles,left,right)
        write_gfa(destination/f'{locus.name}.gfa',nodes,edges,paths)
    write_svg(destination/f'{locus.name}.svg',locus.name,alleles,reason)
    rows = [dict(allele=a.name,sample=a.sample,sequence=a.sequence,read_support=len(a.read_ids),
                 platform_support=a.support,molecule_ids=a.read_ids,seed_molecules=a.seed_molecules,
                 seed_source=a.seed_source,call_status=a.call_status,identity=a.identity,
                 observed_length_range=[min((len(r.sequence) for r in reads if r.molecule_key in a.read_ids),default=0),
                                        max((len(r.sequence) for r in reads if r.molecule_key in a.read_ids),default=0)],
                 consensus_residual=a.consensus_residual,platform_consistency=a.platform_consistency,ambiguous_molecules=a.ambiguous_ids,
                 uncertain=a.uncertain,flags=sorted(set(a.flags)),unit_counts=token_counts(a.tokens),
                 base_metrics=a.annotation.get('metrics',{}),annotation_score_gap=a.annotation.get('score_gap'),
                 phase_policy='fixed_template_phase',compound_counts=dict(Counter(h['catalog_motif'] for h in a.compound_details if 'catalog_motif' in h)),
                 compound_details=a.compound_details,dbg_k=a.dbg_k,graph_support=a.graph_support) for a in candidates]
    fields = list(rows[0]) if rows else ['allele','sample','call_status']
    write_tsv(destination/'candidates.tsv',fields,rows)
    write_tsv(destination/'audit.tsv',fields,[row for row in rows if row['read_support'] >= args.min_support])
    write_tsv(destination/'tokens.tsv',['allele','family','exact_bases','start','end','phase','phase_policy','orientation','edit_operations','status','annotation_ambiguity','level'],
              [dict(asdict(t),allele=a.name) for a in alleles for t in a.tokens])
    write_tsv(destination/'compound_blocks.tsv',['allele','start_token','end_token','start','end','pattern','count','units_per_block','base_unit_count','catalog_motif','basis','exact_bases'],
              [dict(h,allele=a.name) for a in alleles for h in a.compound_details])
    (destination/'annotation_alternatives.json').write_text(json.dumps({a.name:a.annotation for a in alleles},ensure_ascii=False,indent=2),encoding='utf-8')
    write_tsv(destination/'bam_audit.tsv',['bam','sample','platform','spanning_observations','reasons','molecule_decisions'],per_bam)
    status = dict(locus=locus.name,status='input_error' if failures else 'supported_sequences' if alleles else 'unresolved',
                  errors=failures,reason=reason,molecules=len(reads),candidate_count=len(candidates),path_count=len(alleles),
                  uncertain_paths=sum(a.uncertain for a in alleles),boundary_status=locus.boundary_status,
                  input_records=sum(row['reasons'].get('input_records',0) for row in per_bam),
                  molecule_decisions=dict(Counter(row['status'] for row in details)))
    (destination/'status.json').write_text(compact(status),encoding='utf-8')
    print(compact(status),flush=True)
    return bool(failures)


def main():
    parser = argparse.ArgumentParser(description='VNTR point graph with exact sequence and molecule provenance')
    parser.add_argument('run',nargs='?',default='run')
    parser.add_argument('--bam-dir',type=Path,default=Path('bamdata'))
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--catalog',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--locus',action='append',default=[])
    parser.add_argument('--samtools',default=shutil.which('samtools') or 'samtools')
    parser.add_argument('--min-support',type=int,default=3)
    parser.add_argument('--anchor-length',type=int,default=500,choices=vntranchors.LENGTHS,
                        help='Maximum flank template length; retain all 150/250/350/500 lengths up to this cap')
    parser.add_argument('--min-mapq',type=int,default=20)
    parser.add_argument('--ploidy',type=int,default=2)
    parser.add_argument('--score-margin',type=float,default=.5)
    parser.add_argument('--sources',type=Path,help='Optional TSV: bam, namespace; relative paths resolve beside this TSV')
    parser.add_argument('--reference-equivalence',type=Path,
                        help='Optional documented, base-identical ALT-to-primary coordinate blocks')
    args = parser.parse_args()
    if args.min_support < 1 or not 0 <= args.min_mapq <= 254 or not 1 <= args.ploidy <= 4 or args.score_margin <= 0:
        parser.error('Invalid support, anchor, MAPQ, ploidy or score margin')
    if not Path(args.samtools).exists() and not shutil.which(str(args.samtools)):
        parser.error('samtools not found')
    args.source_namespaces = {}
    if args.sources:
        with args.sources.open(encoding='utf-8',newline='') as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                bam = Path(row['bam'])
                key = str((bam if bam.is_absolute() else args.sources.parent/bam).resolve())
                if key in args.source_namespaces:
                    parser.error('Duplicate BAM in sources manifest')
                args.source_namespaces[key] = row['namespace']
    loci = read_catalog(args.catalog)
    if args.locus:
        loci = [x for x in loci if x.name in args.locus]
        if len(loci) != len(set(args.locus)):
            parser.error('Locus missing from BED')
    for locus in loci:
        if (args.out/locus.name).exists() and any((args.out/locus.name).iterdir()):
            parser.error(f'Output already exists: {args.out/locus.name}; choose a new run directory')
    args.out.mkdir(parents=True,exist_ok=True)
    with IndexedFastaReference(args.reference) as reference:
        for locus in loci:
            if locus.chrom not in reference.entries or locus.start < 150 or locus.end+150 > reference.contig_length(locus.chrom):
                raise ValueError(f'Locus or flanks outside reference: {locus.name}')
        print('Scanning every flank 31-mer across the supplied full reference',flush=True)
        equivalences=vntranchors.load_equivalences(args.reference_equivalence,reference)
        profiles,probes = vntranchors.reference_profiles(args.reference,loci,reference,
                                                       [n for n in vntranchors.LENGTHS if n<=args.anchor_length],equivalences)
        (args.out/'anchor_reference.json').write_text(compact(dict(profiles=profiles,probes=probes)),encoding='utf-8')
        failures = sum(process_locus(locus,args,reference,profiles) for locus in loci)
    return int(bool(failures))


if __name__ == '__main__':
    if '--structure-view' in sys.argv:
        sys.argv.remove('--structure-view')
        from vntrstructure import main as structure_main
        raise SystemExit(structure_main())
    if '--read-view' in sys.argv:
        sys.argv.remove('--read-view')
        if len(sys.argv) > 1 and sys.argv[1] == 'run':
            sys.argv.pop(1)
        from vntrreads import main as read_view_main
        raise SystemExit(read_view_main())
    raise SystemExit(main())
