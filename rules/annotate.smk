"""
VEP annotation + parsing of input VCF

    snakemake -s rules/annotate.smk --configfile config/annotate.yaml
"""

import re

VEP = config["vep"]
OUT = config["outdir"]
VCFS = config["vcfs"]
VEP_VCF = f"{OUT}/vep/{{name}}.vep.vcf"
VEP_CSV = f"{OUT}/csv/{{name}}.vep.csv"

wildcard_constraints:
    name="|".join(re.escape(n) for n in VCFS),

localrules: annotate_all

rule annotate_all:
    input:
        expand(VEP_CSV, name=VCFS),


rule vep_annotate:
    input:
        lambda w: VCFS[w.name],
    output:
        VEP_VCF,
    log:
        f"{OUT}/logs/{{name}}.vep.log",
    threads: 8
    resources:
        mem_mb=24000,
    params:
        container=VEP["container"],
        vep_dir=VEP["vep_dir"],
        resources_dir=VEP["resources_dir"],
        fasta=VEP["fasta"],
        loftee_repo=VEP["loftee_repo"],
        loftee_data=VEP["loftee_data"],
        samtools=VEP["samtools"],
        cache_version=VEP["cache_version"],
        bin_dir=lambda w: f"/scratch/$USER/bin/annotate/{w.name}",
    shell:
        r"""
        exec > {log} 2>&1
        export SINGULARITY_TMPDIR=/scratch/$USER/sing_tmp
        export SINGULARITY_CACHEDIR=/scratch/$USER/sing_cache
        mkdir -p $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR {params.bin_dir}
        # per-input copy: concurrent jobs sharing one /scratch staging path race ("Text file busy")
        cp {params.samtools} {params.bin_dir}/samtools
        chmod +x {params.bin_dir}/samtools

        LR={params.loftee_repo}
        LD={params.loftee_data}
        singularity run --pwd "$PWD" -B "$PWD":"$PWD" -H "$PWD":"$PWD" \
            --bind {params.resources_dir}:/opt/vep/resources \
            --bind {params.vep_dir}:/opt/vep/.vep \
               {params.container} vep \
            --dir /opt/vep/.vep \
            -i $PWD/{input} \
            -o $PWD/{output} \
            --fork {threads} \
            --force_overwrite --offline --cache --cache_version {params.cache_version} \
            --format vcf --vcf --everything --canonical --mane \
            --assembly GRCh38 --species homo_sapiens \
            --fasta {params.fasta} \
            --vcf_info_field ANN \
            --plugin NMD \
            --plugin REVEL,/opt/vep/.vep/revel/revel_grch38.tsv.gz \
            --plugin SpliceAI,snv=/opt/vep/.vep/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/opt/vep/.vep/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \
            --plugin UTRAnnotator,/opt/vep/.vep/Plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt \
            --custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar.autogvp.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,AutoGVP \
            --plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
            --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz \
            --plugin LoF,loftee_path:$LR,human_ancestor_fa:$LD/human_ancestor.fa.gz,conservation_file:$LD/loftee.sql,gerp_bigwig:$LD/gerp_conservation_scores.homo_sapiens.GRCh38.bw
        """


rule parse_vep:
    input:
        vcf=VEP_VCF,
        script="scripts/vep_vcf_parser.py",
    output:
        VEP_CSV,
    log:
        f"{OUT}/logs/{{name}}.parse.log",
    params:
        env=VEP["parser_env"],
        mode=config.get("parser_mode", "no_sample"),
    resources:
        mem_mb=8000,
    shell:
        r"""
        exec > {log} 2>&1
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate {params.env}
        python {input.script} -i {input.vcf} -o {output} -m {params.mode}
        """
