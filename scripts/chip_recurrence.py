"""
Aggregates a per-variant recurrence table for CHIP candidates.
variant_id = CHROM_POS_REF_ALT

Rows used: wl_whitelist or wl_manualreview are TRUE
"""
import csv
import argparse
from statistics import median

OUT_COLS = ['variant_id', 'Gene', 'HGVSp', 'Variant.Consequence', 'n_carriers', 'n_carriers_all', 'frequency',
            'median_vaf', 'min_vaf', 'max_vaf', 'median_ad', 'n_batches', 'wl_manualreview', 'recurrent_flag']

def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def fmt(x, spec='.3f'):
    return format(x, spec) if x is not None else '.'


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', required=True, help='chip_calls.csv (chip.smk chip_gather)')
    p.add_argument('--sample-qc', required=True, help='aggregate.smk sample_qc.csv')
    p.add_argument('--min-carriers', type=int, required=True)
    p.add_argument('-o', '--output', required=True)
    a = p.parse_args(argv)

    with open(a.sample_qc, newline='') as fh:
        qc = list(csv.DictReader(fh))
    qc_pass = {r['sample_id'] for r in qc if r['sample_flags'] == 'PASS'}
    batch = {r['sample_id']: r['batch'] for r in qc} if qc and 'batch' in qc[0] else None
    if not qc_pass:
        raise SystemExit(f'{a.sample_qc}: no samples with sample_flags == PASS')

    var = {}        # variant_id -> dict(genes, hgvsp, csq, carriers{sample: (vaf, ad)}, for the entire set)
    with open(a.input, newline='') as fh:
        for r in csv.DictReader(fh):
            if r['wl_whitelist'] != 'TRUE' and r['wl_manualreview'] != 'TRUE':
                continue
            v = var.setdefault(r['ID'], {'genes': [], 'hgvsp': [], 'csq': [], 'carriers': {}, 'all': set(), 'mr': False})
            if r['Gene'] not in v['genes']:
                v['genes'].append(r['Gene'])
                v['hgvsp'].append(r['HGVSp'])
                v['csq'].append(r['Variant.Consequence'])
            v['mr'] |= r['wl_manualreview'] == 'TRUE'
            s = r['Sample.ID']
            v['all'].add(s)
            if s in qc_pass:
                v['carriers'][s] = (to_float(r['Sample.AltFrac']), to_float(r['Sample.AltDepth']))

    n_qc = len(qc_pass)
    rows = []
    for vid, v in var.items():
        vafs = [x[0] for x in v['carriers'].values() if x[0] is not None]
        ads = [x[1] for x in v['carriers'].values() if x[1] is not None]
        n = len(v['carriers'])
        rows.append({
            'variant_id': vid, 'Gene': '&'.join(v['genes']), 'HGVSp': '&'.join(v['hgvsp']),
            'Variant.Consequence': '&'.join(v['csq']),
            'n_carriers': n, 'n_carriers_all': len(v['all']), 'frequency': f'{n / n_qc:.6g}',
            'median_vaf': fmt(median(vafs) if vafs else None), 'min_vaf': fmt(min(vafs, default=None)),
            'max_vaf': fmt(max(vafs, default=None)), 'median_ad': fmt(median(ads) if ads else None, 'g'),
            'n_batches': len({batch.get(s, 'NA') for s in v['carriers']}) if batch is not None else '.',
            'wl_manualreview': 'TRUE' if v['mr'] else 'FALSE',
            'recurrent_flag': 'TRUE' if n >= a.min_carriers else 'FALSE'})
    rows.sort(key=lambda r: (-r['n_carriers'], r['variant_id']))

    with open(a.output, 'w', newline='') as fo:
        w = csv.DictWriter(fo, fieldnames=OUT_COLS)
        w.writeheader()
        w.writerows(rows)
    n_t = sum(r['recurrent_flag'] == 'TRUE' for r in rows)
    print(f'{len(rows)} candidate variants, {n_qc} QC-pass samples; {n_t} with >= {a.min_carriers} '
          f'QC-pass carriers (recurrent_flag) -> {a.output}')


if __name__ == '__main__':
    main()
