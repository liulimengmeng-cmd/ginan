#!/usr/bin/env python3
"""Account for proven raw-observation gaps without inventing solver receipts."""
from pathlib import Path
import collections, datetime as dt, json, re
import run_users as u

R = Path(__file__).resolve().parents[1]


def collect_evidence(out, case, start, count):
    assert u.sha(case['observation']) == case['observation_sha256']
    events = collections.defaultdict(list)
    trace_paths = sorted(out.glob('Network-*.TRACE'))
    assert trace_paths
    for path in trace_paths:
        with path.open(errors='replace') as f:
            for line in f:
                assert 'Failed to get ERP' not in line and 'extrapolated ERP' not in line, 'Invalid runtime EOP'
                tag = line.split(' ', 1)[0]
                if tag in u.TAGS:
                    entry = u.event(line)
                    if 'time' in entry:
                        events[tag].append(entry)
    raw_epochs = set()
    raw_count = 0
    with Path(case['observation']).open() as f:
        for line in f:
            if line.startswith('>'):
                q = line.split()
                instant = dt.datetime(*map(int, q[1:6])) + dt.timedelta(seconds=float(q[6]))
                raw_epochs.add(str(instant))
                raw_count += 1
    log = (out / 'run.log').read_text(errors='replace')
    assert 'PEA finished' in log and f'Exiting at epoch {count} as epoch count {count} has been reached' in log
    assert not any(x in log for x in ['what():', 'INVALID_DGEMV', 'Unable to open Zhang', 'Segmentation fault', 'Failed to get ERP', 'extrapolated ERP'])
    starts = [int(x) for x in re.findall(r'^Starting epoch #(\d+)\s*$', log, re.M)]
    assert starts == list(range(1, count + 1)), ('scheduler_starts', len(starts), count)
    processed = {}
    zero_receiver_epochs = set()
    lines = log.splitlines()
    for i, line in enumerate(lines):
        match = re.match(r'Processed epoch #(\d+) - GPS time: .* - (\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\.00 ', line)
        if match:
            n, epoch = int(match[1]), match[2]
            assert epoch not in processed
            assert epoch == str(start + dt.timedelta(seconds=30 * (n - 1)))
            processed[epoch] = line
            if i > 0 and lines[i - 1].startswith('Synced 0/1 receivers after '):
                zero_receiver_epochs.add(epoch)
    grid = [str(start + dt.timedelta(seconds=30 * i)) for i in range(count)]
    # The final epoch exits before the normal Processed-epoch log line.
    assert set(processed) in (set(grid), set(grid[:-1]))
    positions = set()
    position_paths = sorted(out.glob('*.POS'))
    for path in position_paths:
        for line in path.read_text().splitlines():
            if re.match(r'^\s*2026-', line):
                epoch = str(dt.datetime.fromisoformat(line.split()[0]))
                assert epoch not in positions
                positions.add(epoch)
    sources = [*trace_paths, *position_paths, out / 'run.log', out / 'execution.yaml', Path(case['observation'])]
    return {
        'grid': grid, 'events': events, 'raw_epochs': raw_epochs,
        'zero_receiver_epochs': zero_receiver_epochs, 'positions': positions,
        'processed': processed, 'raw_observation_epochs': raw_count,
        'source_files': {str(p): u.sha(p) for p in sources},
    }


def reconcile(evidence):
    grid = evidence['grid']
    events = evidence['events']
    native = events['ZHANG_USER_PPP_EPOCH']
    times = [e['time'] for e in native]
    assert len(times) == len(set(times)), 'Duplicate native PPP receipt'
    assert set(times) <= set(grid), 'Native receipt outside requested calendar'
    missing = set(grid) - set(times)
    raw_missing = set(grid) - evidence['raw_epochs']
    assert missing == raw_missing, ('unexplained_missing_PPP_receipt', sorted(missing), sorted(raw_missing))
    assert missing <= evidence['zero_receiver_epochs'], 'Missing zero-receiver scheduler proof'
    assert missing <= set(evidence['processed']), 'Missing processed-epoch proof'
    assert not (missing & evidence['positions']), 'Position exists in alleged no-observation gap'
    for tag, entries in events.items():
        assert not (missing & {e['time'] for e in entries}), ('Solver event exists in alleged no-observation gap', tag)
    bytime = {e['time']: {**e, 'origin': 'NATIVE_PPP_TRACE', 'native_ppp_receipt': True} for e in native}
    for epoch in missing:
        bytime[epoch] = {
            'time': epoch, 'status': 'NO_OBSERVATION',
            'reason': 'RAW_RINEX_EPOCH_ABSENT_AND_SCHEDULER_SYNCED_ZERO_RECEIVERS',
            'measurement_update': '0', 'ar_authorized': '0',
            'origin': 'EXTERNAL_CALENDAR_ACCOUNTING', 'native_ppp_receipt': False,
            'has_position': False, 'counts_as_fresh': False, 'counts_as_fixed': False,
            'breaks_convergence': True,
        }
    assert sorted(bytime) == grid
    return {
        'schema': 'OBSERVATION_AWARE_CALENDAR_V1', 'passed': True,
        'calendar_epochs': len(grid), 'native_PPP_status_epochs': len(native),
        'no_observation_epochs': len(missing), 'raw_observation_epochs': evidence['raw_observation_epochs'],
        'missing_epochs': sorted(missing), 'epochs': [bytime[t] for t in grid],
        'source_files': evidence['source_files'],
        'gap_scheduler_lines': {t: evidence['processed'][t] for t in sorted(missing)},
        'native_solver_receipts_unchanged': True, 'coordinate_or_prediction_fill': False,
        'gap_breaks_convergence': True, 'fresh_and_fixed_denominators_unchanged': True,
        'does_not_assert_2880_native_PPP_receipts': True,
    }


def validate_case_calendar(out, case, start, count, pilot=False):
    evidence = collect_evidence(out, case, start, count)
    events = evidence['events']
    ledger = reconcile(evidence)
    statuses = ledger['epochs']
    assert len(statuses) == count and [e['time'] for e in statuses] == evidence['grid']
    # The following native-measurement acceptance block is copied verbatim
    # from the frozen validator; the test checks its normalized source text.
    updates={e['time'] for e in statuses if e['status']=='MEASUREMENT_UPDATED'}
    assert updates,'No fresh PPP updates'
    assert all(e['measurement_update']==e['ar_authorized']=='0' for e in statuses if e['time'] not in updates)
    assert all(e['status']=='UNINITIALIZED' for e in statuses if e['time']<'2026-10-01 02:00:00'),'Unauthorised warmup fallback'
    commits={e['time'] for e in events['ZHANG_KF_MEASUREMENT_TRANSACTION'] if e.get('boundary')=='B16' and e.get('status')=='AUTHORITATIVE_COMMIT'}
    assert commits==updates,('B16_commit_mismatch',len(commits),len(updates))
    ar=events['ZHANG_USER_GONSIS_AR_SUMMARY'];authorities=events['ZHANG_USER_FLOAT_AUTHORITY']
    if case['mode']=='AR':
     assert authorities and all(e['core_bitwise_unchanged']==e['factor_sequence_unchanged']=='1' for e in authorities)
     assert all(e['feedback']=='EPOCH_ONLY' and e['float_feedback']==e['nis_evaluations']=='0' and e['time'] in updates for e in ar)
     assert all(e['ratio_threshold']=='3' for e in events['ZHANG_USER_GONSIS_POLICY'])
    else:assert not ar and not authorities
    applied=[e for e in events['ZHANG_USER_REGIONAL_EPOCH'] if e.get('status')=='APPLIED' and e['time'] in updates]
    assert applied
    applied_by_system=collections.Counter(e.get('system','UNSPECIFIED') for e in applied)
    assert set(applied_by_system)=={'GAL'},applied_by_system
    assert all(e['full_joint_R_shared']==e['receiver_common_GF_direction_preserved']=='1' and e['user_truth_geometry_used']=='0' for e in applied)
    captures=events['ZHANG_USER_GONSIS_REGIONAL_CAPTURE']
    assert all(e['marker_matched']==e['full_CP_PHASE_IONO_R']=='1' and e['actual_rows']==e['expected_rows'] for e in captures)
    pos=[]
    for p in out.glob('*.POS'):
     for line in p.read_text().splitlines():
      if re.match(r'^\s*2026-',line):pos.append(line)
    assert pos,'No finite position output'
    result={'passed':True,'calendar_epochs':count,'fresh_updates':len(updates),'POS_rows':len(pos),
     'first_fresh_epoch':min(updates),'last_fresh_epoch':max(updates),'regional_fresh_receipts_by_system':dict(applied_by_system),
     'AR_summary_systems':sorted(set(e['system'] for e in ar)),'postfix_regional_gate_receipts':len(captures),
     'FLOAT_bitwise_preserved':True,'source_warmup_not_filled':True,'scientific_accuracy_accepted':False}
    # End of the unchanged acceptance block.
    ledger_path = out / 'calendar_status.json'
    assert not ledger_path.exists() and not (out / 'runtime_acceptance.json').exists()
    assert not (out / 'event_receipts.json').exists()
    u.save(ledger_path, ledger)
    result.update(
        native_PPP_status_epochs=ledger['native_PPP_status_epochs'],
        no_observation_epochs=ledger['no_observation_epochs'],
        observation_gaps=ledger['missing_epochs'],
        calendar_accounting_file=str(ledger_path), calendar_accounting_sha256=u.sha(ledger_path),
        calendar_accounting_method=ledger['schema'], original_native_calendar_assertion_passed=False,
        coordinate_or_prediction_fill=False, native_solver_receipts_unchanged=True,
    )
    u.save(out / 'event_receipts.json', events)
    u.save(out / 'runtime_acceptance.json', result)
    return result


def load_statuses_for_analysis(root, events, runtime, grid):
    native = events['ZHANG_USER_PPP_EPOCH']
    native_bytime = {e['time']: e for e in native}
    assert len(native_bytime) == len(native)
    path = root / 'calendar_status.json'
    if not path.exists():
        assert sorted(native_bytime) == grid
        return native_bytime
    assert u.sha(path) == runtime['calendar_accounting_sha256']
    ledger = u.load(path)
    assert ledger['passed'] and ledger['schema'] == 'OBSERVATION_AWARE_CALENDAR_V1'
    assert ledger['coordinate_or_prediction_fill'] is False
    assert ledger['native_PPP_status_epochs'] == len(native)
    assert ledger['no_observation_epochs'] == len(ledger['missing_epochs'])
    assert [e['time'] for e in ledger['epochs']] == grid
    for entry in ledger['epochs']:
        if entry['native_ppp_receipt']:
            original = {k: v for k, v in entry.items() if k not in ('origin', 'native_ppp_receipt')}
            assert original == native_bytime[entry['time']]
        else:
            assert entry['time'] not in native_bytime and entry['time'] in ledger['missing_epochs']
            assert entry['status'] == 'NO_OBSERVATION'
            assert entry['measurement_update'] == entry['ar_authorized'] == '0'
            assert not entry['has_position'] and not entry['counts_as_fresh'] and not entry['counts_as_fixed']
            assert entry['breaks_convergence']
    for path, digest in ledger['source_files'].items():
        assert u.sha(path) == digest, ('Calendar evidence changed', path)
    return {e['time']: e for e in ledger['epochs']}
