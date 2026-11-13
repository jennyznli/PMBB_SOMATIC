"""
CHIP whitelist annotator for VEP-parsed rows 

Whitelist:
  Full_CHIP_gene_list_*.txt
  CHIP_missense_vars_cv_*.txt 
  CHIP_nonsense_FS_vars_cv_*.txt 
  CHIP_splice_vars_cv_*.txt 
  chip_whitelist_rules_*.tsv 

Added columns:
  wl_whitelist     meets any whitelist criteria
  wl_mis           missense variant
  wl_lof           loss-of-function variant
  wl_splice        splice variant
  wl_exception     exception
  wl_manualreview  needs manual review
  wl_rule          custom written fields read in from chip_whitelist_rules_*.tsv
  wl_transcript    MANE.Select == the gene list's Accession
"""

import os
import re
import csv
import glob
import argparse

REQUIRED=['Gene','Variant.Consequence','Amino.acids','Protein.position','EXON','INTRON','MANE.Select']

WL_COLS=['wl_whitelist','wl_mis','wl_lof','wl_splice','wl_exception','wl_manualreview','wl_rule','wl_transcript']

# maps VEP consequence term -> whitelist class ('&' consequences are split)
CLASS_OF={'missense_variant':'missense',
        'stop_gained':'nonsense',
        'frameshift_variant':'frameshift',
        'splice_acceptor_variant':'splice',
        'splice_donor_variant':'splice',
        'inframe_deletion':'inframe_deletion',
        'inframe_insertion':'inframe_insertion'}
CLASS_ALIASES={'lof':{'nonsense','frameshift'},'inframe':{'inframe_deletion','inframe_insertion'}}
VALID_CLASSES=set(CLASS_OF.values())|set(CLASS_ALIASES)

# whitelist file reader, stores rows as a list of dicts
def _rows(path):
    with open(path,newline='') as fh:
        lines=[l.replace('\r','') for l in fh if l.strip() and not l.startswith('#')]
    return [{k.strip():(v or '').strip() for k,v in r.items()} for r in csv.DictReader(lines,delimiter='\t')]

# returns the latest file matching a glob pattern
def _latest(d,pattern):
    hits=sorted(glob.glob(os.path.join(d,pattern)))
    if not hits:
        raise SystemExit(f'{d}: no file matching {pattern}')
    return hits[-1]

# parses protein position ranges as a list of (lo, hi) tuples
def _ranges(s):
    out=[]
    for part in (p.strip() for p in s.split(',') if p.strip()):
        lo,_,hi=part.partition('-')
        out.append((int(lo),int(hi or lo)))
    return out

# returns the first integer in a string, or None if none is found
def _first_int(s):
    m=re.match(r'\s*(\d+)',s or '')
    return int(m.group(1)) if m else None

# returns the first exon or intron number from an EXON / INTRON field
def _region_num(s):
    return _first_int(re.split(r'[|/]',s or '')[0])

# returns the whitelist classes for a VEP consequence string
def var_classes(consequence):
    return {CLASS_OF[t] for t in (consequence or '').split('&') if t in CLASS_OF}

# turns VEP amino acid and protein position fields into a single-residue protein change (compare w whitelist)
def protein_change(amino_acids,protein_position):
    ref,sep,alt=(amino_acids or '').partition('/')
    pos=(protein_position or '').strip()
    if sep and len(ref)==1 and len(alt)==1 and ref!='-' and alt!='-' and pos.isdigit():
        return f'{ref}{pos}{alt}'
    return None

class Whitelist:
    def __init__(self,wl_dir,rules_path=None):
        self.paths={'genes':_latest(wl_dir,'Full_CHIP_gene_list_*.txt'),
                    'missense':_latest(wl_dir,'CHIP_missense_vars_cv_*.txt'),
                    'lof':_latest(wl_dir,'CHIP_nonsense_FS_vars_cv_*.txt'),
                    'splice':_latest(wl_dir,'CHIP_splice_vars_cv_*.txt'),
                    'rules':rules_path or _latest(wl_dir,'chip_whitelist_rules_*.tsv')}
        self.accession={r['Gene']:r['Accession'] for r in _rows(self.paths['genes'])}
        self.missense={(r['Gene'],r['AAChange']) for r in _rows(self.paths['missense'])}
        self.lof_genes={r['Gene'] for r in _rows(self.paths['lof'])}
        self.splice_genes={r['Gene'] for r in _rows(self.paths['splice'])}
        self.rules=[]

        # parse the custom rules file
        for i,r in enumerate(_rows(self.paths['rules']),1):
            classes=set()

            # split the comma-separated list of variant classes
            for c in (x.strip() for x in r['var_class'].split(',') if x.strip()):
                if c not in VALID_CLASSES:
                    raise SystemExit(f'{self.paths["rules"]} row {i}: unknown var_class {c!r}')
                classes|=CLASS_ALIASES.get(c,{c})

            # check the action
            if r['action'] not in ('whitelist','manualreview'):
                raise SystemExit(f'{self.paths["rules"]} row {i}: action must be whitelist / manualreview')
            
            # create a label for the rule and add a restriction if present
            label=f"{r['gene']}_{r['var_class']}"
            restriction=r['exons'] or r['aa_ranges']
            if restriction:
                label+=f"_{restriction}"

            self.rules.append(dict(gene=r['gene'],classes=classes,exons=_ranges(r['exons']),
                                   aa=_ranges(r['aa_ranges']),action=r['action'],
                                   label=label))
                                
        unknown=sorted({r['gene'] for r in self.rules}-set(self.accession))
        if unknown:
            raise SystemExit(f'{self.paths["rules"]}: genes not in the gene list: {unknown}')

        # gene + class pairs a restricted whitelist rule owns -> gene-wide lists skip them
        self.restricted={(r['gene'],c) for r in self.rules if r['action']=='whitelist' and (r['exons'] or r['aa'])
                         for c in r['classes']}

    # checks whether a variant matches a specific, custom rule
    def _rule_hits(self,rule,classes,row):
        if not rule['classes']&classes:
            return False

        if rule['aa']:
            # check if the protein position matches 
            pos=_first_int(row['Protein.position'])
            if pos is None or not any(lo<=pos<=hi for lo,hi in rule['aa']):
                return False

        if rule['exons']:
            # if it has an INTRON value
            if 'splice' in classes and row['INTRON'] not in ('','.'):
                n=_region_num(row['INTRON'])          # intron n = between exons n and n+1
                # check if the intron is adjacent to a whitelisted exon
                ok=n is not None and any(lo-1<=n<=hi for lo,hi in rule['exons'])
            else:
                # check if the exon is in the whitelist
                n=_region_num(row['EXON'])
                ok=n is not None and any(lo<=n<=hi for lo,hi in rule['exons'])
            if not ok:
                return False
            
        return True

    def annotate(self,row):
        # initialize whitelist status and matching rules
        out={'wl_whitelist':False,'wl_mis':False,'wl_lof':False,'wl_splice':False,'wl_exception':False,
             'wl_manualreview':False,'wl_rule':[],'wl_transcript':'.'}

        gene=row['Gene']

        if gene in self.accession:
            classes=var_classes(row['Variant.Consequence'])
            mane=(row['MANE.Select'] or '').split('.')[0]

            # TRUE / FALSE = MANE transcript (version dropped) matches the whitelist Accession; '.' = no MANE
            if mane:
                out['wl_transcript']=mane==self.accession[gene]
            
            # check if missense and protein change is in the missense whitelist
            pc=protein_change(row['Amino.acids'],row['Protein.position'])
            if pc and (gene,pc) in self.missense:
                out['wl_mis']=True
                out['wl_rule'].append('missense')
            
            # check if variant is in a LoF gene
            lof=classes&CLASS_ALIASES['lof']
            if lof and gene in self.lof_genes and not any((gene,c) in self.restricted for c in lof):
                out['wl_lof']=True
                out['wl_rule'].append('lof')

            # check if variant is in a splice gene
            if 'splice' in classes and gene in self.splice_genes and (gene,'splice') not in self.restricted:
                out['wl_splice']=True
                out['wl_rule'].append('splice')
            
            # check if variant matches a custom rule
            for rule in self.rules:
                # go through each rule & skip if it does not match
                if rule['gene']!=gene or not self._rule_hits(rule,classes,row):
                    continue
                
                if rule['action']=='manualreview':
                    out['wl_manualreview']=True
                    out['wl_rule'].append('manual_review')
                    continue
                
                out['wl_exception']=True

                out['wl_rule'].append(f"{rule['label']}")

                hit=rule['classes']&classes

                out['wl_mis']|='missense' in hit
                out['wl_lof']|=bool(hit&CLASS_ALIASES['lof'])
                out['wl_splice']|='splice' in hit
            
            # overall whitelist flag
            out['wl_whitelist']=out['wl_mis'] or out['wl_lof'] or out['wl_splice'] or out['wl_exception']
        
        # format output values
        out['wl_rule']='; '.join(out['wl_rule']) or '.'

        return {k:('TRUE' if v is True else 'FALSE' if v is False else v) for k,v in out.items()}

# adds whitelist columns to a csv
def main(argv=None):
    p=argparse.ArgumentParser(description='Add CHIP whitelist columns (wl_*) to a VEP-parsed csv.')
    p.add_argument('-i','--input',required=True)
    p.add_argument('-o','--output',required=True)
    p.add_argument('--whitelist-dir',default='whitelist')
    p.add_argument('--rules',default=None,help='default: latest WHITELIST_DIR/chip_whitelist_rules_*.tsv')
    args=p.parse_args(argv)
    wl=Whitelist(args.whitelist_dir,args.rules)
    n=n_wl=n_mr=0
    with open(args.input,newline='') as fh, open(args.output,'w',newline='') as fo:
        reader=csv.DictReader(fh)
        missing=[c for c in REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f'{args.input} is missing columns: {missing}')
        w=csv.DictWriter(fo,fieldnames=reader.fieldnames+WL_COLS,quoting=csv.QUOTE_NONNUMERIC)
        w.writeheader()
        for row in reader:
            row.update(wl.annotate(row))
            n+=1
            n_wl+=row['wl_whitelist']=='TRUE'
            n_mr+=row['wl_manualreview']=='TRUE'
            w.writerow(row)
    print(f'{n} rows: {n_wl} whitelisted, {n_mr} manual review -> {args.output}')
    print('whitelist files: '+', '.join(f'{k}={v}' for k,v in wl.paths.items()))


if __name__=='__main__':
    main()
