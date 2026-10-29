"""Evidence-backed structure classes, local modules and reference exploration.

Classes describe ordered motif runs, not inferred genotypes. All memberships
retain a fragment ID and molecular support is deduplicated by sample/namespace/QNAME.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from html import escape
from pathlib import Path
import json

from vntrcore import canonical_motif, compact, revcomp
from vntrmotifs import exact_runs
from vntrreads import row_molecule_key


# 由结构内容生成稳定短 ID；ID 表达内容分组，不表示生物学等位基因编号。
def identity(prefix, value):
    return prefix+sha256(compact(value).encode()).hexdigest()[:12]


# 只合并相邻同标签的完整单元；未知和部分单元保持可见边界。
def blocks(tokens):
    result=[]
    for i,t in enumerate(tokens):
        label=t['motif']+('r' if t['orientation']=='-' else '')
        if t['status']=='unknown':label='U'
        if t['status']=='partial':label+='p'
        if result and t['status']=='unit' and result[-1]['status']=='unit' and result[-1]['label']==label:
            result[-1].update(count=result[-1]['count']+1,end=i+1,bp_end=t['end'])
        else:result.append(dict(label=label,count=1,start=i,end=i+1,bp_start=t['start'],bp_end=t['end'],status=t['status']))
    return result


# 骨架用 ×n 忽略次数差异，详细表达保留实际次数。
def structure_text(bs, skeleton=False):
    return ' → '.join(b['label']+('×n' if skeleton else '×'+str(b['count'])) if b['status']=='unit'
                     else f"{b['label']}({b['bp_end']-b['bp_start']}bp)" for b in bs)


# 记录参考序列、窗口精确周期和来源；目录边界标签不等于自动确认边界。
def reference_report(locus, sequence, window, window_start, dictionary, segment, strand, fasta):
    tokens=segment((sequence,dictionary,'HiFi'));bs=blocks(tokens)
    for t in tokens:
        t['reference_start']=locus.start+t['start'];t['reference_end']=locus.start+t['end']
    repeats=exact_runs(window)
    for hit in repeats:
        hit['reference_start']=window_start+hit['start'];hit['reference_end']=window_start+hit['end']
        hit['crosses_core_boundary']=any(hit['reference_start']<p<hit['reference_end'] for p in (locus.start,locus.end))
    stat=Path(fasta).stat()
    return dict(chrom=locus.chrom,start=locus.start,end=locus.end,coordinate_system='0-based half-open',
                catalog_boundary_status=locus.boundary_status,gene_strand=strand,sequence=sequence,
                gene_sequence=revcomp(sequence) if strand=='-' else sequence if strand=='+' else None,
                sequence_sha256=sha256(sequence.encode()).hexdigest(),fasta=str(Path(fasta).resolve()),
                fasta_size=stat.st_size,fasta_mtime_ns=stat.st_mtime_ns,
                fai_sha256=sha256(Path(str(fasta)+'.fai').read_bytes()).hexdigest(),
                window_start=window_start,window_end=window_start+len(window),window_sequence=window,
                tokens=tokens,blocks=bs,structure=structure_text(bs),exact_period_runs=repeats,
                unknown_bases=sum(len(t['exact_bases']) for t in tokens if t['status']=='unknown'),
                boundary_review='crossing_exact_run' if any(h['crosses_core_boundary'] for h in repeats) else 'no_exact_crossing_detected',
                method='seed-phase interruptions + independent exact periods 1..100bp, >=3 copies and >=12bp; final dictionary DP',
                limitations='Reference is one sequence, not population structure. Boundary status is catalog metadata; no automatic confirmation. Exact scan can miss imperfect repeats. Gene strand must be supplied.')


# 用于强证据筛选：主比对、MAPQ≥20，并排除重复、QC失败和锚点次序冲突。
def is_primary(f):
    return f['alignment_type']=='primary' and f['mapq']>=20 and not set(f.get('flags',[]))&{'duplicate_flag','qc_fail','anchor_order_conflict_cigar_fallback'}


# 完整结构强证据还要求 HiFi、双锚点及两侧至少 90% 的局部锚点匹配。
def strong_boundary(row,f):
    return row['platform']=='HiFi' and is_primary(f) and f['coverage']=='spans_candidate_interval' and f['anchor']=='both' and all(x is not None and x>=.9 for x in f.get('anchor_scores',[1,1]))


# 发现实际片段中的复合串联块和相邻游程窗口，不跨越未知、残片或反向候选。
def local_candidates(tokens,max_unit_tokens=12):
    """Tandem patterns and short ordered run windows, never crossing U/P.

    Tandem patterns use their primitive token period and greedy longest span.
    Shared local blocks have 2..3 adjacent runs; boundary runs are intact.
    """
    labels=[t['motif'] if t['status']=='unit' and t['orientation']=='+' else None for t in tokens]
    result=[];i=0
    while i<len(labels):
        choices=[]
        for k in range(2,min(max_unit_tokens,(len(labels)-i)//2)+1):
            unit=labels[i:i+k]
            if None in unit or len(set(unit))<2:continue
            if any(k%p==0 and unit==unit[:p]*(k//p) for p in range(1,k)):continue
            count=1
            while labels[i+count*k:i+(count+1)*k]==unit:count+=1
            if count>=2:choices.append((k*count,-k,count,unit))
        if choices:
            width,negk,count,unit=max(choices);k=-negk
            # Phase equivalence only for the module ID; observed phase retained.
            canonical=min(tuple(unit[j:]+unit[:j]) for j in range(k))
            result.append(dict(kind='tandem',key=list(canonical),pattern=unit,count=count,start=i,end=i+width))
            i+=width
        else:i+=1
    bs=blocks(tokens)
    for i in range(len(bs)):
        for width in (2,3):
            part=bs[i:i+width]
            if len(part)!=width or any(b['status']!='unit' or b['label'].endswith('r') for b in part):continue
            key=[[b['label'],b['count']] for b in part]
            result.append(dict(kind='shared_block',key=key,pattern=key,count=1,start=part[0]['start'],end=part[-1]['end']))
    return result


# 按有序骨架分类型、按次数组合分亚型；同一分子出现冲突时保留为未分配。
def structural_summary(rows,minimum=3,min_fraction=.10):
    classes={};modules={};assignments=[];denominators=defaultdict(set);molecule_paths=defaultdict(set)
    observations=[]
    for row in rows:
        molecular=row_molecule_key(row)
        for f in row['fragments']:
            if strong_boundary(row,f):denominators[row['sample']].add(molecular)
            bs=blocks(f['tokens']);reason=[]
            if f['coverage']!='spans_candidate_interval':reason.append('partial_coverage')
            if any(t['status']!='unit' for t in f['tokens']):reason.append('unresolved_or_partial_unit')
            if any(t['orientation']=='-' for t in f['tokens']):reason.append('reverse_candidate')
            if any(t.get('annotation_ambiguity') for t in f['tokens']):reason.append('annotation_ambiguity')
            if not bs:reason.append('empty')
            cid=identity('S',[b['label'] for b in bs]) if not reason else ''
            counts=[b['count'] for b in bs];sid=identity('L',[cid,counts]) if cid else ''
            if cid and is_primary(f):molecule_paths[molecular].add((cid,sid))
            observations.append((row,f,molecular,bs,cid,sid,reason))
    conflicts={key for key,paths in molecule_paths.items() if len(paths)>1}
    module_seen=set()
    for row,f,molecular,bs,cid,sid,reason in observations:
        conflict=molecular in conflicts
        if conflict:reason.append('conflicting_molecule_paths');cid=sid=''
        high=strong_boundary(row,f) and not conflict
        evidence='strong' if cid and high else 'supporting' if cid else 'unassigned'
        assignment=dict(row_id=row['id'],read_id=row['read_id'],fragment_id=f['id'],sample=row['sample'],
                        namespace=row.get('namespace',''),platform=row['platform'],molecule=molecular,class_id=cid,subtype_id=sid,
                        evidence=evidence,reason=','.join(reason) or ('boundary_or_platform_limited' if not high else ''),length=len(f['sequence']),module_ids=[])
        f['class_id']=cid;f['subtype_id']=sid;f['summary_evidence']=evidence
        assignments.append(assignment)
        if cid:
            if cid not in classes:classes[cid]=dict(id=cid,skeleton=structure_text(bs,True),labels=[b['label'] for b in bs],members=[],subtypes={})
            group=classes[cid];group['members'].append(assignment)
            sub=group['subtypes'].setdefault(sid,dict(id=sid,counts=[b['count'] for b in bs],structure=structure_text(bs),members=[],examples=[]))
            sub['members'].append(assignment)
            # Choose an observed read; prioritize strong boundary then edits.
            rank=(not high,sum(t['edit_distance'] for t in f['tokens']),row['read_id'],f['id'])
            sub['examples'].append((rank,dict(row_id=row['id'],fragment_id=f['id'],sequence=f['sequence'],structure=structure_text(bs),blocks=bs)))
        for hit in local_candidates(f['tokens']):
            ts=f['tokens'][hit['start']:hit['end']]
            if any(t.get('annotation_ambiguity') for t in ts):continue
            mid=identity('B',[hit['kind'],hit['key']])
            exact=all(t['edit_distance']==0 and t['exact_bases']==t['template'] for t in ts)
            instance=dict(assignment,start=ts[0]['start'],end=ts[-1]['end'],count=hit['count'],pattern=hit['pattern'],exact=exact,
                          sequence=''.join(t['exact_bases'] for t in ts),eligible=row['platform']=='HiFi' and is_primary(f) and not conflict)
            instance.pop('module_ids')
            # One molecule, local position and exact slice is one occurrence.
            occurrence=(mid,molecular,f['qstart']+instance['start'],instance['sequence'])
            instance['occurrence_key']=identity('O',occurrence)
            key=(occurrence,row['id'],f['id'])
            if key not in module_seen:
                module_seen.add(key)
                if mid not in modules:
                    text=('['+' → '.join(hit['key'])+']×k' if hit['kind']=='tandem' else ' → '.join(f'{m}×{n}' for m,n in hit['key']))
                    modules[mid]=dict(id=mid,kind=hit['kind'],structure=text,instances=[])
                modules[mid]['instances'].append(instance)
            if mid not in assignment['module_ids']:assignment['module_ids'].append(mid)
    def counts(members):
        bysample=defaultdict(set);byplatform=defaultdict(set)
        for a in members:bysample[a['sample']].add(a['molecule']);byplatform[a['platform']].add(a['molecule'])
        return dict(support=len({a['molecule'] for a in members}),samples={s:len(v) for s,v in sorted(bysample.items())},platforms={p:len(v) for p,v in sorted(byplatform.items())})
    finished=[]
    for group in classes.values():
        members=group.pop('members');group.update(counts(members))
        strong=defaultdict(set)
        for a in members:
            if a['evidence']=='strong':strong[a['sample']].add(a['molecule'])
        group['hifi_support']=len(set().union(*strong.values())) if strong else 0
        group['sample_evidence']={s:dict(support=len(strong[s]),denominator=len(denominators[s]),fraction=len(strong[s])/max(1,len(denominators[s]))) for s in group['samples']}
        # 完整类型要求至少一个样本同时满足支持条数和该样本内比例两个门槛。
        group['status']='supported' if any(v['support']>=minimum and v['fraction']>=min_fraction for v in group['sample_evidence'].values()) else 'candidate'
        subs=[]
        for sub in group['subtypes'].values():
            sub.update(counts(sub.pop('members')));sub['representative']=min(sub.pop('examples'),key=lambda v:v[0])[1];subs.append(sub)
        group['subtypes']=sorted(subs,key=lambda s:(-s['support'],s['id']))
        group['representative']=group['subtypes'][0]['representative']
        lengths=[a['length'] for a in members];group['length_min']=min(lengths);group['length_max']=max(lengths)
        group['run_ranges']=[[min(s['counts'][i] for s in subs),max(s['counts'][i] for s in subs)] for i in range(len(group['labels']))]
        group['composition_bp']={}
        for b in group['representative']['blocks']:group['composition_bp'][b['label']]=group['composition_bp'].get(b['label'],0)+b['bp_end']-b['bp_start']
        finished.append(group)
    completed_modules=[]
    for module in modules.values():
        instances=module['instances'];module.update(counts(instances))
        eligible={a['molecule'] for a in instances if a['eligible']}
        # 局部模块使用全批次合格分子总数，不能套用上面的完整类型样本内门槛。
        exact={a['molecule'] for a in instances if a['exact']}
        module.update(hifi_support=len(eligible),exact_support=len(exact),approximate_support=len({a['molecule'] for a in instances if not a['exact']}),occurrences=len({a['occurrence_key'] for a in instances}),
                      status='supported' if len(eligible)>=minimum else 'candidate')
        # Exact/approximate sets may overlap if a read has several occurrences.
        completed_modules.append(module)
    finished.sort(key=lambda g:(g['status']!='supported',-g['hifi_support'],-g['support'],g['id']))
    completed_modules.sort(key=lambda g:(g['status']!='supported',-g['hifi_support'],-g['support'],g['id']))
    audit=dict(class_count=len(finished),supported_classes=sum(g['status']=='supported' for g in finished),module_count=len(completed_modules),
               supported_modules=sum(g['status']=='supported' for g in completed_modules),assigned_fragments=sum(bool(a['class_id']) for a in assignments),
               unassigned_fragments=sum(not a['class_id'] for a in assignments),conflicting_molecules=len(conflicts),
               eligible_hifi_by_sample={s:len(v) for s,v in sorted(denominators.items())},minimum_support=minimum,minimum_sample_fraction=min_fraction,
               module_search='primitive tandem periods 2..12 tokens; shared blocks 2..3 intact runs; no U/P crossing; all candidates retained',
               support_definition='distinct sample/namespace/QNAME, not traversal count; complete classes and local modules have different coverage requirements',
               limitations='Observed structure classes, not phased alleles or population frequencies. Local exact/approximate support can overlap. No imputation across U/P or partial coverage.')
    return dict(classes=finished,modules=completed_modules,assignments=assignments,audit=audit)


# 静态摘要最多展示 20 个达到门槛的类型；完整候选仍在 HTML 和 TSV 中。
def summary_svg(summary,locus):
    groups=[g for g in summary['classes'] if g['status']=='supported'][:20]
    height=max(170,85+len(groups)*95)
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" width="1500" height="{height}" viewBox="0 0 1500 {height}"><rect width="100%" height="100%" fill="white"/>',
         f'<text x="24" y="30" font-size="22">{escape(locus)} observed structure paths</text>',
         '<text x="24" y="55" font-size="13">Separate context per path. Representative observed counts; schematic widths. Full candidates in HTML/TSV.</text>']
    if not groups:out.append('<text x="24" y="110">No class meets the reporting evidence threshold; all candidate paths remain available.</text>')
    for i,g in enumerate(groups):
        y=100+i*95
        out.append(f'<text x="24" y="{y}" font-size="13">{escape(g["id"])} | {g["support"]} reads / {len(g["samples"])} samples</text>')
        out.append(f'<text x="24" y="{y+25}" font-size="13">{escape(g["representative"]["structure"])}</text>')
        out.append(f'<text x="24" y="{y+47}" font-size="12">{g["length_min"]}..{g["length_max"]} bp | {len(g["subtypes"])} observed joint count combinations</text>')
    return ''.join(out)+'</svg>'
