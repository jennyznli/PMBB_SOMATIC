"""
CHIP calling from cohort VEP csv

    snakemake -s rules/chip.smk --configfile config/chip.yaml

  chip_seq_filter  CH genes, MANE, DP / alt depth / VAF
  chip_whitelist   whitelist flags
  chip_recurrence  carriers per variant (whitelist + manual_review)
  chip_artifact    batch / age / germline tests 

"""

C = config
PY = C.get("python", "python3")
D = C["outdir"]


localrules: chip_all

rule chip_all:
    input:
        f"{D}/chip_artifact.csv",


rule chip_seq_filter:
    input:
        csv=C["input"],
        samples=C.get("samples_from") or [],
        script="scripts/chip_seq_filtering.py",
    output:
        csv=temp(f"{D}/chip_seq.csv"),
        stats=f"{D}/chip_seq_stats.tsv",
        qc=f"{D}/sample_qc.csv",
    log:
        f"{D}/logs/chip_seq_filter.log",
    params:
        samples=lambda w, input: f"--samples-from {input.samples}" if input.samples else "",
    resources:
        mem_mb=8000,
    shell:
        r"""
        {PY} {input.script} -i {input.csv} {params.samples} --sample-qc-out {output.qc} \
            --whitelist-dir {C[whitelist_dir]} --min-dp {C[min_dp]} \
            --min-alt-depth {C[min_alt_depth]} --min-vaf {C[min_vaf]} \
            -o {output.csv} --stats {output.stats} > {log} 2>&1
        """


rule chip_whitelist:
    input:
        csv=rules.chip_seq_filter.output.csv,
        script="scripts/chip_whitelist.py",
    output:
        f"{D}/chip_calls.csv",
    log:
        f"{D}/logs/chip_whitelist.log",
    shell:
        "{PY} {input.script} -i {input.csv} -o {output} --whitelist-dir {C[whitelist_dir]} > {log} 2>&1"


rule chip_recurrence:
    input:
        calls=rules.chip_whitelist.output[0],
        qc=rules.chip_seq_filter.output.qc,
        script="scripts/chip_recurrence.py",
    output:
        f"{D}/chip_recurrence.csv",
    log:
        f"{D}/logs/chip_recurrence.log",
    shell:
        r"""
        {PY} {input.script} -i {input.calls} --sample-qc {input.qc} \
            --min-carriers {C[min_carriers]} -o {output} > {log} 2>&1
        """


rule chip_artifact:
    input:
        calls=rules.chip_whitelist.output[0],
        rec=rules.chip_recurrence.output[0],
        qc=rules.chip_seq_filter.output.qc,
        covariates=C["covariates_tsv"] or [],
        script="scripts/chip_artifact.R",
    output:
        variants=f"{D}/chip_artifact.csv",
        calls=f"{D}/chip_calls.artifact.csv",
    log:
        f"{D}/logs/chip_artifact.log",
    params:
        cov=f"--covariates {C['covariates_tsv']}" if C["covariates_tsv"] else "",
        exempt=",".join(C.get("exempt_genes") or []),
        batch_test=C.get("batch_test", True),
        age_test=C.get("age_test", True),
        germline_test=C.get("germline_test", True),
    shell:
        r"""
        Rscript {input.script} --calls {input.calls} --recurrence {input.rec} --sample-qc {input.qc} \
            {params.cov} --p-age {C[p_age]} --p-batch {C[p_batch]} \
            --batch-test {params.batch_test} --age-test {params.age_test} \
            --germline-test {params.germline_test} --p-binom {C[p_binom]} --gnomad-max-af {C[gnomad_max_af]} \
            --exempt "{params.exempt}" \
            --out-variants {output.variants} --out-calls {output.calls} > {log} 2>&1
        """
