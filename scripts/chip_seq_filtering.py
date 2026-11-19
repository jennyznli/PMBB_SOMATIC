"""
CHIP sequence-level filtering on a cohort VEP csv

Steps:
  01  Gene in the whitelist gene list (latest *_chip_rheum_genes.txt), Bio.type == protein_coding,
      Sample.AltFrac / Sample.Depth / Sample.AltDepth all present (not '', '.', 'NA', '-1')
  02  MANE.Select or MANE.PlusClinical filled
  03  Sample.Depth >= --min-dp
  04  Sample.AltDepth >= --min-alt-depth
  05  Sample.AltFrac > --min-vaf (strict)
"""
import os
import sys
import csv
import argparse
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chip_whitelist import GENES_PATTERN, _latest, _gene_accessions  # noqa: E402

STEPS = ['s01_gene_biotype_values', 's02_mane', 's03_dp', 's04_alt_depth', 's05_vaf']
VALUE_COLS = ('Sample.AltFrac', 'Sample.Depth', 'Sample.AltDepth')
MISSING = ('', '.', 'NA', '-1')


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def filled(x):
    return (x or '').strip() not in ('', '.')


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', required=True, help='cohort vep csv')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('--stats', required=True)
    p.add_argument('--sample-qc-out', required=True)
    p.add_argument('--samples-from', default='', help='csv w/ Sample.ID; adds samples w/ no rows in --input')
    p.add_argument('--drop-cols', default='', help='comma-separated input cols left out of --output')
    p.add_argument('--whitelist-dir', default='whitelist')
    p.add_argument('--min-dp', type=int, default=20)
    p.add_argument('--min-alt-depth', type=int, default=3)
    p.add_argument('--min-vaf', type=float, default=0.02)
    a = p.parse_args(argv)

    gene_file = _latest(a.whitelist_dir, GENES_PATTERN)
    genes = set(_gene_accessions(gene_file))
    drop = set(filter(None, a.drop_cols.split(',')))

    n, samples = Counter(), set()
    with open(a.input, newline='') as fi, open(a.output, 'w', newline='') as fo:
        r = csv.DictReader(fi)
        w = csv.DictWriter(fo, fieldnames=[c for c in r.fieldnames if c not in drop], extrasaction='ignore')
        w.writeheader()
        for x in r:
            n['input'] += 1
            samples.add(x['Sample.ID'])
            if x['Gene'] not in genes:
                continue
            n['in_gene'] += 1
            dp, ad, vaf = (to_float(x[k]) for k in ('Sample.Depth', 'Sample.AltDepth', 'Sample.AltFrac'))
            ok = [x['Bio.type'] == 'protein_coding' and all(x[k].strip() not in MISSING for k in VALUE_COLS),
                  filled(x['MANE.Select']) or filled(x['MANE.PlusClinical']),
                  dp is not None and dp >= a.min_dp,
                  ad is not None and ad >= a.min_alt_depth,
                  vaf is not None and vaf > a.min_vaf]
            cum = True
            for step, s_ok in zip(STEPS, ok):
                cum = cum and s_ok
                n[step] += cum
            if cum:
                w.writerow(x)

    if a.samples_from:
        with open(a.samples_from, newline='') as fh:
            samples |= {x['Sample.ID'] for x in csv.DictReader(fh)}
    with open(a.sample_qc_out, 'w') as fh:
        fh.write('sample_id,sample_flags\n' + ''.join(f'{s},PASS\n' for s in sorted(samples)))

    order = ['input', 'in_gene'] + STEPS
    with open(a.stats, 'w') as fs:
        fs.write('step\tn\n' + ''.join(f'{k}\t{n[k]}\n' for k in order))
    print(f'{a.input}: ' + ', '.join(f'{k}={n[k]}' for k in order) + f' -> {a.output}')
    print(f'{len(samples)} samples; min_dp={a.min_dp} min_alt_depth={a.min_alt_depth} min_vaf>{a.min_vaf}; '
          f'{len(genes)} genes in {gene_file}')


if __name__ == '__main__':
    main()
