"""Opt-in, fixed-phase D1S80 right-terminal templates from independent molecules.

Candidates never nominate or vote. Reference nominates but never votes. Accepted
templates remain sample-specific and must not enter the ordinary motif dictionary.
"""
from collections import Counter, defaultdict
from hashlib import sha256
import json

from vntrcore import canonical_motif, edit_distance, motif_orientation, revcomp
from vntrmotifs import primitive_unit

VERSION = 'd1s80_right_terminal_v1'
UNIT_BP = 16


# 验证样本特异的右端模板：参考可提名，只有合格独立 OBS 可以投票。
def validate_terminals(records, molecules, decisions, templates, right_flank,
                       chrom, repeat_start, repeat_end, min_support=3):
    """Return candidate rows, per-molecule evidence and accepted sample templates.

    Every support requires the unchanged query slice, a complete anchor decision,
    the original inner-repeat dictionary and an exact external 16-bp context.
    Identical duplicated observations vote once; inconsistent observations never
    vote. Unknown/failed contexts are retained with reasons, not counted as alleles.
    """
    if min_support < 3 or len(right_flank) != UNIT_BP or set(right_flank)-set('ACGT'):
        raise ValueError('Require >=3 independent molecules and an exact 16-bp right context')
    fixed = set(templates) | {revcomp(m) for m in templates}
    if not fixed or any(len(m) != UNIT_BP for m in fixed):
        raise ValueError('Require original 16-bp inner-repeat templates')
    # 重复元数据和重复记录先按分子归组，不能增加独立支持数。
    by_molecule, by_decision, observations = defaultdict(list), defaultdict(list), defaultdict(list)
    for row in molecules:
        by_molecule[(row['sample'], row['molecule_key'])].append(row)
    for row in decisions:
        by_decision[row['molecule_key']].append(row)
    for row in records:
        if row['kind'] == 'observed':
            observations[(row['sample'], row['molecule_key'])].append(row)
    references = [r for r in records if r['kind'] == 'reference']
    if len(references) != 1 or len(references[0]['sequence']) != repeat_end-repeat_start:
        raise ValueError('Require one reference matching the unchanged catalog interval')
    # 提名保留精确方向与相位，不把循环移位后的末端合并计票。
    ref_terminal = references[0]['sequence'][-UNIT_BP:]
    nominated = {ref_terminal: {'reference'}}
    nomination_molecules = defaultdict(lambda: defaultdict(set))
    assessed = {}
    for (sample, key), paths in sorted(observations.items()):
        full = [r for r in paths if r['coverage'] == 'two_flank_anchored']
        result = dict(sample_id=sample, molecule_key=key,
                      path_ids_json=sorted(r['id'] for r in paths),
                      source_records_json=sorted({s for r in paths for s in r.get('source_records', [])}),
                      platforms_json=sorted({p for r in paths for p in r.get('platforms', [])}),
                      terminal_bases='', inner_distance='NA', outer_exact=False,
                      qL='NA', qR='NA', boundary_evidence_valid=False, reason='local_only')
        assessed[(sample, key)] = result
        if not full:
            continue
        result['source_records_json'] = sorted({s for r in full for s in r.get('source_records', [])})
        result['platforms_json'] = sorted({p for r in full for p in r.get('platforms', [])})
        if any(r.get('sequence_sha256') != sha256(r['sequence'].encode()).hexdigest() for r in full):
            result['reason'] = 'frozen_sequence_hash_mismatch'
            continue
        audit = by_decision.get(key, [])
        if not audit or any(r.get('status') != 'complete_support' or r.get('sample') != sample for r in audit):
            result['reason'] = 'missing_or_conflicting_anchor_decision'
            continue
        if any(r.get('anchor_status') != 'complete_support' for r in full):
            result['reason'] = 'noncomplete_frozen_anchor_status'
            continue
        # Multiple metadata rows are not extra molecules. Disagreement is rejected.
        # 同一分子的完整序列、坐标或来源不一致时，整分子退出支持。
        signatures = {(r['sequence'], tuple(r.get('stored_query_interval', [])),
                       tuple(r.get('original_query_interval', [])),
                       tuple(sorted(r.get('source_records', [])))) for r in full}
        stored = by_molecule.get((sample, key), [])
        if len(signatures) != 1 or not stored:
            result['reason'] = 'conflicting_or_missing_complete_observation'
            continue
        source = full[0]
        if (not source.get('source_records') or
                any(not set(source['source_records']).issubset(r.get('sources', [])) for r in audit)):
            result['reason'] = 'source_records_absent_from_anchor_decision'
            continue
        # 核对保存的完整 query、切片坐标、原始方向坐标及来源记录。
        verified = []
        for row in stored:
            a, b = int(row['qL']), int(row['qR'])
            whole, sequence = row['query_sequence'], row['sequence']
            valid = (0 <= a < b <= len(whole) and whole[a:b] == sequence == source['sequence']
                     and [a, b] == source.get('stored_query_interval')
                     and json.loads(row['original_query_interval']) == source.get('original_query_interval')
                     and sorted(json.loads(row['records'])) == sorted(source.get('source_records', [])))
            if not valid:
                break
            verified.append((a, b, whole[b:b+UNIT_BP], sequence))
        if len(verified) != len(stored) or len(set(verified)) != 1:
            result['reason'] = 'source_slice_or_boundary_mismatch'
            continue
        a, b, outer, sequence = verified[0]
        result.update(qL=a, qR=b, boundary_evidence_valid=True)
        if len(sequence) < 2*UNIT_BP:
            result['reason'] = 'insufficient_inner_context'
            continue
        terminal = sequence[-UNIT_BP:]
        result['terminal_bases'] = terminal
        nominated.setdefault(terminal, set()).add('complete_observed')
        nomination_molecules[terminal][sample].add(key)
        # 内侧 16 bp 须接近原重复模板，外侧 16 bp 须与真实参考上下文精确一致。
        result['inner_distance'] = min(edit_distance(sequence[-2*UNIT_BP:-UNIT_BP], m) for m in fixed)
        result['outer_exact'] = outer == right_flank
        result['reason'] = ('inner_repeat_not_supported' if result['inner_distance'] > 1 else
                            'right_context_not_exact' if not result['outer_exact'] else 'eligible_context')
    samples = sorted({r['sample'] for r in records if r['kind'] != 'reference'})
    candidates, evidence, accepted = [], [], defaultdict(dict)
    # Existing exact oriented templates need no new terminal-only transition.
    # 只有新的、原始周期为 16 bp 的末端才可能成为专用终端模板。
    for terminal, sources in sorted(nominated.items()):
        if len(terminal) != UNIT_BP or set(terminal)-set('ACGT') or terminal in fixed:
            continue
        identifier = 'D1S80_R_' + sha256(terminal.encode()).hexdigest()[:16]
        family = canonical_motif(terminal)
        phase, orientation = motif_orientation(terminal, family)
        for sample in samples:
            # 门槛逐样本计算；平台小计可能重叠，不能相加当作总分子数。
            support, platform_support, reasons = set(), defaultdict(set), Counter()
            for (s, key), observation in sorted(assessed.items()):
                if s != sample:
                    continue
                reason = observation['reason']
                if reason == 'eligible_context':
                    reason = ('exact_terminal_with_context' if observation['terminal_bases'] == terminal
                              else 'different_terminal_state')
                if len(primitive_unit(terminal)) != UNIT_BP:
                    reason = 'nonprimitive_terminal_template'
                supports = reason == 'exact_terminal_with_context'
                if supports:
                    support.add(key)
                    for platform in observation['platforms_json']:
                        platform_support[platform].add(key)
                reasons[reason] += 1
                evidence.append(dict(terminal_id=identifier, **observation,
                                     decision='support' if supports else 'uninformative',
                                     validation_reason=reason))
            status = ('rejected_nonprimitive' if len(primitive_unit(terminal)) != UNIT_BP else
                      'accepted' if len(support) >= min_support else 'insufficient_independent_support')
            candidates.append(dict(locus_id='D1S80', terminal_id=identifier, role='terminal', side='right',
                template_bases=terminal, template_length=UNIT_BP, family=family,
                phase=phase, orientation=orientation, chrom=chrom,
                reference_start0=repeat_end-UNIT_BP if terminal == ref_terminal else 'NA',
                reference_end0=repeat_end if terminal == ref_terminal else 'NA',
                source_json=sorted(sources), sample_id=sample,
                nominated_observed_molecules=len(nomination_molecules[terminal][sample]),
                support_molecules=len(support), platform_counts_json={p: len(v) for p, v in sorted(platform_support.items())},
                min_support=min_support, status=status, rejection_reasons_json=dict(sorted(reasons.items())),
                reference_terminal=terminal == ref_terminal, dictionary_version=VERSION))
            # 接受结果仅供该样本完整区间的右端使用，不扩充普通内部模板字典。
            if status == 'accepted':
                accepted[sample][terminal] = identifier
    return candidates, evidence, dict(accepted)
