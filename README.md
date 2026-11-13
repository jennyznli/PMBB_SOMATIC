# PMBB somatic + CHIP pipeline

Proposed workflow for calling somatic variants (tumor-only Mutect2) on PMBB whole-exome CRAMs and building a curated CHIP callset from them.

Two releases:

1. **Somatic callset** Every PMBB sample with light flagged QC and VEP annotation, for general use.
2. **CHIP callset** A curated subset in CH genes, with flags for artifacts and germline variants.

| Stage | Rule file | Status |
|---|---|---|
| 0. QC | `rules/qc.smk` | proposed |
| 1. Calling | `rules/mutect2.smk` | written |
| 2. Light filter | `rules/filter.smk` | proposed |
| 3. Annotation | `rules/annotate.smk` | written |
| 4. Aggregate / release | `rules/aggregate.smk` | proposed |
| 5. CHIP curation | `rules/chip.smk` | written |


### Things to clarify: 
- Panel of normals (from Daehyun potentially) needed for Mutect2 calling
- Review CHIP whitelist with Daehyun

## Layout

```
config/      one yaml per rule
rules/       Snakemake rule files
scripts/     helper scripts
whitelist/   CHIP gene list, variant whitelists, rules tsv (the newest file by date is used)
```
---

## 0. QC (`qc.smk`)

Catch bad samples before calling, and their results go into `sample_qc`.

| Check | Tool | What it catches |
|---|---|---|
| File integrity | `samtools quickcheck`, contigs vs reference `.fai` | truncated CRAMs, wrong reference build |
| Read groups | `samtools view -H` (`@RG`) | missing or merged read groups |
| Mapping / duplicates | `samtools flagstat`, `samtools stats -t targets.bed` | low % mapped, high duplicate rate |
| Coverage | `mosdepth` | mean target coverage, fraction of targets at ≥10/20/50/100x |
| Sex check | `samtools idxstats` (chrX/chrY ratio) | sample swaps |
| Report | `multiqc` | one HTML covering all samples |

## 1. Somatic variant calling (`mutect2.smk`)

This stage follows the GATK somatic short-variant best practices in tumor-only mode.

**Panel of normals (PoN).** Each normal is run through Mutect2 (`--max-mnp-distance 0`) and FilterMutectCalls. The results are combined with GenomicsDBImport and CreateSomaticPanelOfNormals. 

**Calling.** Mutect2 reassembles reads in active regions to find candidate SNVs and indels. It scores each candidate (TLOD) using the PoN and the gnomAD germline resource as priors.

**Orientation artifact model.** Mutect2 writes F1R2 counts: for each alt allele, which read orientation supports it. LearnReadOrientationModel learns how often each base change is seen on only one orientation. One-sided support is the signature of damage artifacts such as FFPE deamination (C>T) and oxidation/OxoG (G>T).

**Contamination.** GetPileupSummaries is run at common germline SNPs. At those sites, reads that don't fit the sample's own genotype point to DNA from another person. CalculateContamination turns this into a contamination estimate plus tumor segmentation. `gather_contamination` collects the per-sample estimates into one table. 

**Filtering.** FilterMutectCalls combines the contamination table, the segmentation and the orientation priors (`--ob-priors`) and sets FILTER. Its output is one filtered VCF per sample.

## 2. Light filter (`filter.smk`)

This stage acts on each sample's VCF:

- `bcftools norm`: split multiallelics and left-align indels.
- Flag, but don't remove, calls that fail these minimums:
  - DP ≥ 20
  - alt AD ≥ 3
  - alt support on both read orientations (F1R2 ≥ 1 and F2R1 ≥ 1)
- The Mutect2 FILTER is carried forward but not required; can be helpful for later review. 

Output: one normalized, flagged VCF per sample. 

## 3. Annotation (`annotate.smk`)

- Collect the unique variant sites across all samples into a sites-only VCF (no sample columns), one per chromosome.
- Run VEP once per chromosome on that VCF. It annotates with CANONICAL and MANE_SELECT transcripts, REVEL, SpliceAI, AlphaMissense, LOFTEE, UTRAnnotator, NMD, MaveDB, ClinVar and AutoGVP.
- `vep_vcf_parser.py` flattens the VEP output to CSV.

The CHIP stage requires reads a cohort-mode CSV (`parser_mode: cohort`), which has one row per sample and variant.

## 4. Aggregate / release (`aggregate.smk`, proposed)

| Table | One row per | Columns |
|---|---|---|
| `somatic_variants` | sample × variant (one file per chromosome) | sample_id, variant_id, gene, VAF, DP, AD, FILTER, QC flags |
| `variant_catalog` | variant | variant_id, gene, consequence, VEP annotations, n_carriers, frequency |
| `sample_qc` | sample | sample_id, coverage, contamination, batch, sex check, QC status |

`sample_qc` = stage 0 metrics + stage 1 `contamination.tsv`, joined on `sample_id`.

`somatic_variants` doesn't repeat the annotations; join it to `variant_catalog` on `variant_id` to get them. `variant_id` = `chr_pos_ref_alt`.

## 5. CHIP curation (`chip.smk`)

Every step after the sequence filter adds flag columns and never removes rows.

### 5a. Sequence filter (`chip_seq_filtering.py`)

| Step | Keep if |
|---|---|
| s01 | gene is on the CHIP gene list, `Bio.type == protein_coding`, and AltFrac/Depth/AltDepth are present |
| s02 | MANE Select or MANE Plus Clinical transcript |
| s05 | `Sample.AltFrac > min_vaf` (default 0.02) |

Step counts are written to `chip_seq_stats.tsv`.

### 5b. Whitelist (`chip_whitelist.py`)

Each call is matched against the files in `whitelist/`, which adds these columns:

| Flag | Meaning |
|---|---|
| `wl_mis` / `wl_lof` / `wl_splice` | matches the missense, nonsense/frameshift or splice whitelist |
| `wl_exception` | gene-specific exception |
| `wl_rule` | which rule from `chip_whitelist_rules_*.tsv` matched |
| `wl_transcript` | whether the MANE transcript matches the transcript the whitelist refers to |
| `wl_whitelist` | CHIP candidate |
| `wl_manualreview` | possible candidate that needs manual review |

### 5c. Recurrence (`chip_recurrence.py`)

For each candidate variant (`wl_whitelist` or `wl_manualreview`), this step summarizes the carriers:

- `n_carriers` (QC-pass samples), `n_carriers_all`, `frequency`
- median, min and max VAF; median alt depth; `n_batches`
- `recurrent_flag`: set when `n_carriers ≥ min_carriers`. A variant shared by many people is more likely an artifact than a real clone.

### 5d. Artifact and germline flags (`chip_artifact.R`)

| Flag | Applies to | Rule |
|---|---|---|
| `batch_flag` | recurrent variants | Fisher test of carrier status vs batch, BH-adjusted p < `p_batch`. Real CHIP shouldn't cluster in one batch. |
| `age_flag` | recurrent variants, except `exempt_genes` | logistic regression of carrier status on age plus covariates. The flag is set unless OR > 1 and p < `p_age`. |
| `germline_binom` | every candidate call | VAF can't be told apart from 0.5 (binomial p ≥ `p_binom`) |
| `germline_gnomad` | every candidate call | gnomAD max AF > `gnomad_max_af` |
| `germline_m2` | every candidate call | Mutect2 FILTER includes `germline` |
| `germline_flag` | every candidate call | any of the three germline flags |
| `review` | variant | `recurrent_flag` or `wl_manualreview` |

Each test can be switched off in `config/chip.yaml`. The batch test needs a `batch` column in `sample_qc`. The age test needs `covariates_tsv`, which has sample_id, age and any other covariates and is kept outside the repo.

Outputs:
- `chip_artifact.csv`: one row per variant, with all flags. This is the manual-review table.
- `chip_calls.artifact.csv`: one row per call, with the flags joined on.

### 5e. Final callset (`chip_final`, proposed)

- Apply the final cutoffs (min AD, recurrence, germline, manual-review decisions) to produce the CHIP callset.
- Optional OncoKB annotation, run once on the merged CSV. 

--- 
