"""Sequence-only D1S80 unit annotation; independent of graph construction."""
from collections import defaultdict
from hashlib import sha256
from vntrcore import decompose, revcomp

VERSION = "d1s80_joint_units_v2"
TERMINAL_VERSION = "d1s80_joint_units_terminal_v1"
UNIT_BP = 16

def template_for(token):
    if token['status'] == 'unknown':
        return ''
    family = token['family']
    oriented = family if token['orientation'] == '+' else revcomp(family)
    phase = token['phase']
    return oriented[phase:] + oriented[:phase]


# 只选原始周期确为 16 bp 的可用基本模板，不把较短周期的倍数误作 16 bp 单元。
def primitive_templates(motifs):
    from vntrmotifs import primitive_unit
    return sorted({row['motif'] for row in motifs
                   if row['usable'] in ('True', True) and row['level'] == 'base'
                   and len(row['motif']) == UNIT_BP
                   and len(primitive_unit(row['motif'])) == UNIT_BP})


# 给所有保存路径增加单元注释层；原始碱基、边界和候选序列均不改写。
def primitive_units(records, motifs, legacy_tokens, terminal_by_sample=None):
    """Annotate every saved sequence; fragment counts never become locus CN.

    The core decomposer retains two paths, not every tied annotation. A tied
    result therefore remains ambiguous even when the two displayed counts agree.
    """
    templates = primitive_templates(motifs)
    if not templates:
        raise ValueError('No usable 16-bp D1S80 templates')
    legacy = defaultdict(list)
    for row in legacy_tokens:
        legacy[row['allele']].append(row)
    units, summaries, counts = [], [], []
    for record in records:
        path_id = record['id']
        is_candidate = record['kind'] == 'candidate'
        name = path_id.removeprefix('CAND_') if is_candidate else 'NA'
        sequence = record['sequence']
        sequence_hash = sha256(sequence.encode()).hexdigest()
        if record.get('sequence_sha256', sequence_hash) != sequence_hash:
            raise ValueError(f'{path_id}: sequence hash does not match input')
        coverage = record.get('coverage', 'unknown')
        # 局部 OBS 的单元数只描述该片段，不能报告整个基因组位点的 CN。
        fragment = record['kind'] == 'observed' and coverage != 'two_flank_anchored'
        count_scope = ('observed_fragment' if fragment else
                       'observed_full_interval' if record['kind'] == 'observed' else
                       'candidate_sequence' if is_candidate else 'reference_catalog_interval')
        metadata = dict(path_id=path_id, kind=record['kind'], candidate_id=name,
                        sample=record.get('sample', ''), coverage=coverage,
                        molecule_key=record.get('molecule_key', 'NA'),
                        platforms_json=record.get('platforms', []),
                        sequence_sha256=sequence_hash, count_scope=count_scope)
        old_tokens = legacy.get(name, []) if is_candidate else []
        terminal = {}
        if terminal_by_sample is not None and not fragment:
            terminal = ({m: identifier for group in terminal_by_sample.values() for m, identifier in group.items()}
                        if record['kind'] == 'reference' else terminal_by_sample.get(record.get('sample', ''), {}))
        primary, _, audit = decompose(sequence, tuple(templates), return_details=True,
                                     right_terminal_motifs=tuple(terminal))
        # 保存主方案和一个次优方案；次优存在不一定表示并列最优。
        schemes = {'primary': [dict(vars(item)) for item in primary]}
        if audit['alternative']:
            schemes['alternative'] = [dict(item) for item in audit['alternative']]
        def signature(token):
            return tuple(token[key] for key in ('start', 'end', 'status', 'family',
                                                'orientation', 'phase', 'exact_bases',
                                                'edit_operations'))
        signatures = {scheme: {signature(token) for token in tokens}
                      for scheme, tokens in schemes.items()}
        path_units = []
        for scheme, tokens in schemes.items():
            # 逐方案验证原始序列无损重建及连续坐标，未知片段也必须保留。
            if ''.join(token['exact_bases'] for token in tokens) != sequence:
                raise AssertionError(f'{path_id} {scheme}: units do not reconstruct sequence')
            position = 0
            for index, token in enumerate(tokens, 1):
                if token['start'] != position or token['end'] != position+len(token['exact_bases']):
                    raise AssertionError(f'{path_id} {scheme}: coordinates do not partition sequence')
                if token['exact_bases'] != sequence[token['start']:token['end']]:
                    raise AssertionError(f'{path_id} {scheme}: unit differs from source slice')
                position = token['end']
                parents = [i for i, old in enumerate(old_tokens, 1)
                           if int(old['start']) < token['end'] and int(old['end']) > token['start']]
                other = 'alternative' if scheme == 'primary' else 'primary'
                token_ambiguous = (audit['annotation_ambiguity']
                                   and signature(token) not in signatures.get(other, set()))
                path_units.append(dict(**metadata, decomposition=scheme, unit_index=index,
                                       start0=token['start'], end0=token['end'],
                                       status=token['status'], family=token['family'],
                                       template_bases=template_for(token),
                                       exact_bases=token['exact_bases'],
                                       observed_unit_bp=len(token['exact_bases']),
                                       orientation=token['orientation'], phase_offset=token['phase'],
                                       phase_policy=audit['phase_policy'],
                                       edit_operations=token['edit_operations'],
                                       annotation_ambiguity=token_ambiguous,
                                       legacy_token_indexes_json=parents if is_candidate else 'NA'))
                if terminal_by_sample is not None:
                    path_units[-1]['terminal_template_id'] = (terminal.get(template_for(token), 'NA')
                                                              if token['end'] == len(sequence) else 'NA')
            if scheme == 'primary':
                full = sum(token['status'] == 'full' for token in tokens)
                partial = sum(token['status'] == 'partial' for token in tokens)
                unknown = sum(len(token['exact_bases']) for token in tokens
                              if token['status'] == 'unknown')
                alternate = audit['alternative']
                alt_count = sum(item['status'] == 'full' for item in alternate) if alternate else None
                annotation_status = ('unresolved' if unknown or partial else
                                     'ambiguous' if audit['annotation_ambiguity'] else 'resolved_count')
                count_status = 'fragment_only' if fragment else annotation_status
                summaries.append(dict(**metadata,
                                      sequence_bp=len(sequence), nominal_period_bp=UNIT_BP,
                                      length_equiv_cn=len(sequence)/UNIT_BP,
                                      full_16bp_units=full, partial_units=partial,
                                      unknown_bp=unknown, primitive_cn=full if count_status == 'resolved_count' else 'NA',
                                      primitive_cn_status=count_status,
                                      annotation_count_status=annotation_status,
                                      annotation_ambiguity=audit['annotation_ambiguity'],
                                      alternative_full_units=alt_count if alt_count is not None else 'NA',
                                      primary_score=audit['primary_score'],
                                      alternative_score=audit['alternative_score'] if alternate else 'NA',
                                      score_gap=audit['score_gap'] if alternate else 'NA',
                                      alternative_is_tied=bool(alternate) and audit['annotation_ambiguity'],
                                      alternative_limit=audit['alternative_limit'],
                                      score_definition=audit['score_definition'],
                                      phase_policy=audit['phase_policy'],
                                      legacy_full_tokens=sum(row['status'] == 'full' for row in old_tokens)
                                      if is_candidate else 'NA',
                                      method_version=VERSION))
                # 启用终端模板模式后统一保留 CN=NA；完整注释单元数另列。
                if terminal_by_sample is not None:
                    summaries[-1].update(primitive_cn='NA', cn_policy='withheld_terminal_validation',
                        primitive_cn_status='withheld' if count_status == 'resolved_count' else count_status,
                        annotated_full_units=full, available_terminal_ids_json=sorted(terminal.values()),
                        method_version=TERMINAL_VERSION)
        units.extend(path_units)
        grouped = defaultdict(lambda: [0, 0])
        for row in path_units:
            if row['decomposition'] == 'primary':
                key = (row['family'], row['status'])
                grouped[key][0] += 1
                grouped[key][1] += len(row['exact_bases'])
        for (family, status), (segments, bp) in sorted(grouped.items()):
            counts.append(dict(**metadata,
                               decomposition='primary', family=family, status=status,
                               segment_count=segments, observed_bp=bp))
    return units, summaries, counts
