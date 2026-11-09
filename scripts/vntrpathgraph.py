"""Context-positioned sequence graphs with exact paths and explicit loss rules.

This is a center-star alignment graph, not PGGB and not a genotype caller.
Node sharing is conditioned on the chosen observed backbone's alignment columns.
"""
from collections import Counter, defaultdict
from hashlib import sha256
import json
import re

VERSION = '1.0-path-sequence-graph'


def digest(sequence):
    return sha256(sequence.encode('ascii')).hexdigest()


def select_paths(rows, min_mapq=20, platforms=None):
    """完整双侧区间候选才进入 P；同分子冲突与局部记录保留排除原因。"""
    from vntrreads import row_molecule_key
    candidates, excluded = defaultdict(list), []
    for row in rows:
        molecule = row_molecule_key(row)
        if not row['fragments']:
            excluded.append(dict(molecule=molecule, fragment='', reason=row.get('status', 'no_fragment')))
        for f in row['fragments']:
            reasons = []
            if platforms and row['platform'] not in platforms: reasons.append('platform_not_selected')
            if f.get('alignment_type') != 'primary': reasons.append('not_primary')
            mq = f.get('mapq')
            if mq is None or mq == 255 or not min_mapq <= mq <= 254: reasons.append('mapq_not_qualified')
            if f.get('flag', 0) & (4 | 256 | 512 | 1024 | 2048): reasons.append('alignment_flag')
            scores = f.get('anchor_scores', [])
            from vntrflankalign import core_eligible
            if not core_eligible(f): reasons.append(f.get('extraction_status','incomplete_B2_core'))
            bad = {'reference_skip', 'anchor_order_conflict_cigar_fallback'} & set(f.get('flags', []))
            reasons.extend(sorted(bad))
            sequence = f.get('sequence', '')
            if not re.fullmatch('[ACGTRYSWKMBDHVN]+', sequence): reasons.append('invalid_sequence')
            if reasons:
                excluded.append(dict(molecule=molecule, fragment=f['id'], reason=','.join(reasons)))
            else:
                candidates[molecule].append(dict(sequence=sequence, sample=row['sample'], platform=row['platform'],
                    fragment=f['id'], row_id=row['id'], bam=row.get('bam', ''), molecule=molecule,
                    read_id=row['read_id'], namespace=row.get('namespace', ''),
                    anchor_scores=scores, mapq=mq, qstart=f['qstart'], qend=f['qend'],
                    record_id=f.get('record_id'), strand=f.get('strand'),
                    reference_interval=f.get('reference_interval'), boundary_status=f.get('boundary_status')))
                if f.get('extraction_version'):
                    from vntrflankalign import FIELDS
                    candidates[molecule][-1].update({k:f.get(k) for k in FIELDS})
                    candidates[molecule][-1]['core_candidate_pairs'] = [[p['qstart'],p['qend']] for p in f['core_anchor_evidence']['pairs']]
    grouped = defaultdict(list)
    for molecule, observations in sorted(candidates.items()):
        sequences = {o['sequence'] for o in observations}
        if len(sequences) != 1:
            excluded.extend(dict(molecule=molecule, fragment=o['fragment'], reason='conflicting_complete_sequences')
                            for o in observations)
            continue
        grouped[next(iter(sequences))].extend(observations)
    records = []
    for sequence, observations in grouped.items():
        molecules = sorted({o['molecule'] for o in observations})
        hifi = sorted({o['molecule'] for o in observations if o['platform'] == 'HiFi'})
        samples = {s: len({o['molecule'] for o in observations if o['sample'] == s})
                   for s in sorted({o['sample'] for o in observations})}
        by_platform = {s: len({o['molecule'] for o in observations if o['platform'] == s})
                       for s in sorted({o['platform'] for o in observations})}
        records.append(dict(sequence=sequence, sha256=digest(sequence), molecules=molecules,
                            hifi_molecules=hifi, samples=samples, platforms=by_platform, observations=observations))
    records.sort(key=lambda r: (-len(r['hifi_molecules']), -len(r['molecules']), len(r['sequence']), r['sha256']))
    for i, r in enumerate(records, 1): r['id'] = f'P{i:04d}'
    return records, excluded


def alignment_columns(query, backbone, cigar=None):
    """每个实际碱基拥有上下文坐标；不把远处相同 motif 全局并成一个点。"""
    if cigar is None:
        import edlib
        cigar = edlib.align(query, backbone, mode='NW', task='path')['cigar']
    if not cigar or not re.fullmatch(r'(\d+[=XID])+', cigar):
        raise ValueError('Missing or invalid global alignment CIGAR')
    q = r = 0
    columns = []
    for count, op in re.findall(r'(\d+)([=XID])', cigar):
        count = int(count)
        if count <= 0: raise ValueError('Zero-length CIGAR operation')
        if op in '=X':
            for _ in range(count):
                if q >= len(query) or r >= len(backbone): raise ValueError('Alignment exceeds sequence')
                if (query[q] == backbone[r]) != (op == '='): raise ValueError('CIGAR base mismatch')
                columns.append((r, 1, 0, 0, query[q])); q += 1; r += 1
        elif op == 'I':
            for k in range(count):
                if q >= len(query): raise ValueError('Insertion exceeds query')
                columns.append((r, 0, count, k, query[q])); q += 1
        else:
            r += count
    if q != len(query) or r != len(backbone): raise ValueError('Global alignment does not consume both sequences')
    if ''.join(c[-1] for c in columns) != query: raise ValueError('Alignment changed query bases')
    return columns, cigar


def build_graph(sequences, backbone, cigars=None):
    if not sequences: return dict(nodes={}, edges=set(), paths={}, coordinates={}, alignments={})
    columns, alignments = {}, {}
    for name, sequence in sorted(sequences.items()):
        columns[name], alignments[name] = alignment_columns(sequence, backbone, (cigars or {}).get(name))
    ordered = sorted({k for walk in columns.values() for k in walk})
    index = {k: i+1 for i, k in enumerate(ordered)}
    paths = {p: [index[k] for k in w] for p, w in columns.items()}
    graph = dict(nodes={index[k]: k[-1] for k in ordered},
                 coordinates={index[k]: [list(k[:-1])] for k in ordered},
                 edges={e for w in paths.values() for e in zip(w, w[1:])}, paths=paths, alignments=alignments)
    validate_graph(graph, sequences)
    return graph


def validate_graph(graph, sequences):
    nodes, edges, paths = (graph[k] for k in ('nodes', 'edges', 'paths'))
    if set(paths) != set(sequences): raise ValueError('Graph path identities differ from input')
    if any(not s or not re.fullmatch('[ACGTRYSWKMBDHVN]+', s) for s in nodes.values()):
        raise ValueError('Invalid node sequence')
    if any(a not in nodes or b not in nodes for a, b in edges): raise ValueError('Dangling graph edge')
    used = set()
    for p, walk in paths.items():
        if not walk or any(n not in nodes for n in walk): raise ValueError('Missing path node')
        used.update(zip(walk, walk[1:]))
        if any(e not in edges for e in zip(walk, walk[1:])): raise ValueError('Missing path edge')
        if ''.join(nodes[n] for n in walk) != sequences[p]: raise ValueError('Path sequence changed: '+p)
    if used != edges: raise ValueError('Graph contains an edge absent from all embedded paths')
    return dict(paths=len(paths), nodes=len(nodes), edges=len(edges), reconstructed_bases=sum(map(len, sequences.values())))


def compact_graph(graph):
    """合并非分支链；保留任何 path 的起止，不能跨越分叉或汇合。"""
    nodes, edges, paths = (graph[k] for k in ('nodes', 'edges', 'paths'))
    before, after = defaultdict(set), defaultdict(set)
    for a, b in edges: after[a].add(b); before[b].add(a)
    starts = {w[0] for w in paths.values() if w} | set(graph.get('barriers', ()))
    ends = {w[-1] for w in paths.values() if w} | set(graph.get('barriers', ()))
    def next_node(n):
        if n in ends or len(after[n]) != 1: return None
        target = next(iter(after[n]))
        return target if target not in starts and len(before[target]) == 1 else None
    incoming = {b: a for a in nodes if (b := next_node(a)) is not None}
    heads = [n for n in sorted(nodes) if n not in incoming]
    chains, visited = [], set()
    for head in heads:
        chain = []; node = head
        while node is not None and node not in visited:
            chain.append(node); visited.add(node); node = next_node(node)
        chains.append(chain)
    if len(visited) != len(nodes): raise ValueError('Unexpected cycle in center-star graph')
    mapping = {n: i for i, chain in enumerate(chains, 1) for n in chain}
    compact_paths = {}
    for p, walk in paths.items():
        converted = []
        for n in walk:
            if not converted or converted[-1] != mapping[n]: converted.append(mapping[n])
        compact_paths[p] = converted
    result = dict(nodes={i: ''.join(nodes[n] for n in chain) for i, chain in enumerate(chains, 1)},
                  coordinates={i: [c for n in chain for c in graph.get('coordinates', {}).get(n, [])]
                               for i, chain in enumerate(chains, 1)},
                  edges={(mapping[a], mapping[b]) for a, b in edges if mapping[a] != mapping[b]},
                  paths=compact_paths, alignments=graph.get('alignments', {}),
                  source_nodes={i: chain for i, chain in enumerate(chains, 1)})
    validate_graph(result, {p: ''.join(nodes[n] for n in w) for p, w in paths.items()})
    return result


def apply_merge_rules(records, policy, budget, locus):
    """有损替换必须明确到原始序列、代表序列及允许改变的 0-based 位点。"""
    if not 0 <= budget <= 1: raise ValueError('Loss budget must be between zero and one')
    original = {r['id']: r['sequence'] for r in records}
    represented = dict(original); ledger = []
    if policy is None: return represented, ledger
    if policy.get('schema_version') != 1 or policy.get('locus') != locus:
        raise ValueError('Merge policy schema or locus mismatch')
    by_hash = {r['sha256']: r for r in records}
    rules = policy.get('rules', []); sources = [r['source_sha256'] for r in rules]
    if len(sources) != len(set(sources)): raise ValueError('Duplicate merge source')
    if set(sources) & {r['target_sha256'] for r in rules}: raise ValueError('Chained or cyclic replacement is forbidden')
    for rule in rules:
        if not rule.get('reason'): raise ValueError('Each replacement needs a reason')
        if rule['source_sha256'] not in by_hash or rule['target_sha256'] not in by_hash:
            raise ValueError('Replacement source and target must be observed paths')
        source, target = (by_hash[rule[k]] for k in ('source_sha256', 'target_sha256'))
        a, b = source['sequence'], target['sequence']
        if len(a) != len(b): raise ValueError('Lossy length change is forbidden')
        positions = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
        allowed = rule.get('allowed_positions', [])
        if any(type(i) is not int or i < 0 or i >= len(a) for i in allowed): raise ValueError('Invalid allowed position')
        if not positions or not set(positions).issubset(allowed): raise ValueError('Unapproved base substitution')
        if any(a[i] not in 'ACGT' or b[i] not in 'ACGT' for i in positions): raise ValueError('Cannot replace uncertain bases')
        if len(positions) > budget*len(a)+1e-10: raise ValueError('Per-path loss budget exceeded')
        represented[source['id']] = b
        ledger.append(dict(path=source['id'], target_path=target['id'], substitutions=len(positions),
            fraction=len(positions)/len(a), positions=positions, source_bases=''.join(a[i] for i in positions),
            target_bases=''.join(b[i] for i in positions), reason=rule['reason']))
    return represented, ledger


def graph_support(graph, records):
    members = {r['id']: set(r['molecules']) for r in records}
    node_paths, edge_paths = defaultdict(set), defaultdict(set)
    visits = Counter()
    for p, walk in graph['paths'].items():
        for n in walk: node_paths[n].add(p); visits[n] += 1
        for edge in zip(walk, walk[1:]): edge_paths[edge].add(p)
    def count(names): return len(set().union(*(members[p] for p in names))) if names else 0
    return ({n: dict(paths=sorted(ps), molecules=count(ps), visits=visits[n]) for n, ps in node_paths.items()},
            {e: dict(paths=sorted(ps), molecules=count(ps)) for e, ps in edge_paths.items()})


def write_gfa(path, graph, records, mode, represented, loss_ledger=None):
    nstats, estats = graph_support(graph, records)
    original = {r['id']: r['sequence'] for r in records}
    losses = {p: sum(a != b for a,b in zip(original[p],s)) for p,s in represented.items()}
    explicit = {r['path']:r for r in loss_ledger or []}
    if loss_ledger is not None: losses = {p:explicit.get(p,{}).get('substitutions',0) for p in represented}
    lines = ['H\tVN:Z:1.0\tgm:Z:'+mode,
             '# S contains literal bases; L is 0M adjacency; only embedded P walks are observations.',
             '# mc=distinct molecule keys; pc=distinct raw paths; vc=path visits; lb=explicit substitutions.']
    for n,s in sorted(graph['nodes'].items()):
        v=nstats[n]
        lines.append(f'S\t{n}\t{s}\tLN:i:{len(s)}\tmc:i:{v["molecules"]}\tpc:i:{len(v["paths"])}\tvc:i:{v["visits"]}')
    for (a,b),v in sorted(estats.items()):
        lines.append(f'L\t{a}\t+\t{b}\t+\t0M\tmc:i:{v["molecules"]}\tps:Z:'+','.join(v['paths']))
    by_id={r['id']:r for r in records}
    for p,w in sorted(graph['paths'].items()):
        extra=''
        if loss_ledger is not None:
            loss=explicit.get(p,{})
            extra=(f'\tle:i:{loss.get("loss_bp",0)}\tli:i:{loss.get("insertions",0)}'
                   f'\tld:i:{loss.get("deletions",0)}\tol:i:{len(original[p])}\trl:i:{len(represented[p])}')
        lines.append(f'P\t{p}\t'+','.join(f'{n}+' for n in w)+
                     f'\t*\tmc:i:{len(by_id[p]["molecules"])}\tlb:i:{losses[p]}'+extra)
    path.write_text('\n'.join(lines)+'\n',encoding='ascii')


def validate_gfa(path, expected):
    """重新解析落盘文件，核对标准字段、方向、边和每条路径的实际拼接。"""
    nodes, edges, walks = {}, set(), {}
    for line in path.read_text(encoding='ascii').splitlines():
        fields=line.split('\t')
        if fields[0]=='S':
            n=int(fields[1])
            if n in nodes: raise ValueError('Duplicate segment')
            nodes[n]=fields[2]
            if f'LN:i:{len(fields[2])}' not in fields[3:]: raise ValueError('Invalid segment length')
        elif fields[0]=='L':
            if fields[2]!='+' or fields[4]!='+' or fields[5]!='0M': raise ValueError('Unsupported GFA link')
            edges.add((int(fields[1]),int(fields[3])))
        elif fields[0]=='P':
            if fields[1] in walks or fields[3]!='*': raise ValueError('Invalid path')
            if any(not n.endswith('+') for n in fields[2].split(',')): raise ValueError('Invalid orientation')
            walks[fields[1]]=[int(n[:-1]) for n in fields[2].split(',')]
    return validate_graph(dict(nodes=nodes,edges=edges,paths=walks),expected)
