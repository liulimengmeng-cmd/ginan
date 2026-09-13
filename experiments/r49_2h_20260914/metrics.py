import pathlib, json, csv, collections, datetime
BASE = pathlib.Path(__file__).resolve().parent
def quality(path):
    epochs = collections.defaultdict(list)
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                continue  # a currently incomplete final line is not a complete record
            epochs[int(row['gpst_seconds'])].append(row)
    result = []
    for epoch, rows in sorted(epochs.items()):
        by_sat = collections.defaultdict(lambda: collections.defaultdict(list))
        product_rows = [r for r in rows if r['solution'] == 'PRODUCT_FIXED']
        for row in product_rows:
            by_sat[row['satellite']][row['observable']].append(row)
        sats, components, duplicates = [], collections.defaultdict(list), []
        for sat, signals in by_sat.items():
            if any(len(signals.get(obs, [])) > 1 for obs in ('L1C', 'L2W')):
                duplicates.append(sat)
                continue  # fail closed for ambiguous duplicate records
            if any(len(signals.get(obs, [])) != 1 for obs in ('L1C', 'L2W')):
                continue
            pair = [signals[obs][0] for obs in ('L1C', 'L2W')]
            cid = pair[0]['integer_component_id'].strip()
            if cid in ('', 'NONE') or cid != pair[1]['integer_component_id'].strip():
                continue
            if not all(r['product_state'] == 'AR_VALID' and all(r[k] == '1' for k in
                       ('integer_valid', 'pppar_usable', 'ar_valid', 'dual_frequency_ar_valid')) for r in pair):
                continue
            sats.append(sat)
            components[cid].append(sat)
        result.append(dict(epoch=epoch, time=datetime.datetime.fromtimestamp(epoch, datetime.timezone.utc).strftime('%H:%M:%S'),
                           rows=len(rows), product_fixed_rows=len(product_rows), strict_ar=len(sats),
                           satellites=sorted(sats), largest_component=max(map(len, components.values()), default=0),
                           components=dict(components), duplicate_product_satellites=duplicates,
                           product_states=dict(collections.Counter(r['product_state'] for r in product_rows)),
                           product_invalid_reasons=dict(collections.Counter(r['invalid_reason'] for r in product_rows))))
    return result

def metrics(rows):
    vals = [r['strict_ar'] for r in rows]
    lc = [r['largest_component'] for r in rows]
    first = next((i for i, v in enumerate(vals) if v), len(vals))
    return dict(epochs=len(vals), strict_ar_peak=max(vals, default=0), satellite_epoch_integral=sum(vals),
                largest_component_peak=max(lc, default=0), largest_component_integral=sum(lc),
                zero_epochs=sum(v == 0 for v in vals), zero_epochs_after_first=sum(v == 0 for v in vals[first:]),
                positive_to_zero_transitions=sum(a > 0 and b == 0 for a, b in zip(vals, vals[1:])))

if __name__ == '__main__':
    recorded=json.loads((BASE/'version_quality.json').read_text())
    for version, data in recorded.items():
        actual=metrics(quality(pathlib.Path(data['path'])/'zhang_internal_products.csv'))
        assert actual==data['metrics'], (version, actual, data['metrics'])
        print(version, actual)
    rows=json.loads((BASE/'epoch_comparison.json').read_text())
    for veto in [True, False]:
        group=[x for x in rows if x['r45_veto']==veto]
        print('writer_veto',veto,'epochs',len(group),'r37',sum(x['r37'] for x in group),'r45',sum(x['r45'] for x in group))
