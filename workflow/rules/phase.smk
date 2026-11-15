# Phasing, the final phased VCFs, and their pgen copies.


# The panel and the cohort condition each other. The SHAPEIT5 docs skip phase_rare
# below 2,000 samples. Sites absent from the panel are dropped, and the log counts them.
rule phase_shapeit5:
    input:
        bcf=f"{OUTDIR}/prepare/{{region}}.bcf",
        csi=f"{OUTDIR}/prepare/{{region}}.bcf.csi",
        panel=PANEL,
        panel_index=f"{PANEL}.csi",
        map=config["genetic_maps"]["shapeit5"],
    output:
        f"{OUTDIR}/phase/{{region}}.bcf",
    params:
        region=lambda w: REGION_SPAN[w.region],
    log:
        f"{LOGDIR}/phase_shapeit5/{{region}}.txt",
    benchmark:
        f"{BENCHDIR}/phase_shapeit5/{{region}}.tsv"
    software:
        SHAPEIT5_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        SHAPEIT5_phase_common \
            --input {input.bcf:q} \
            --reference {input.panel:q} \
            --map {input.map:q} \
            --region {params.region:q} \
            --thread {threads} \
            --output {output:q}

        sync {output:q}
        """


rule finalize_phased:
    input:
        f"{OUTDIR}/phase/{{chrom}}.bcf",
    output:
        vcf=f"{OUTDIR}/phased/{{chrom}}.vcf.gz",
        tbi=f"{OUTDIR}/phased/{{chrom}}.vcf.gz.tbi",
    log:
        f"{LOGDIR}/finalize_phased/{{chrom}}.txt",
    benchmark:
        f"{BENCHDIR}/finalize_phased/{{chrom}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bcftools +fill-tags {input:q} -Ou -- --tags AN,AC,AF,MAF \
        | bcftools annotate \
            --set-id '%CHROM:%POS:%REF:%FIRST_ALT' \
            --threads {threads} \
            -Oz \
            --output {output.vcf:q} \
            --write-index=tbi

        echo "Phased records: $(bcftools index --nrecords {output.vcf:q})"

        sync {output:q}
        """


rule export_pgen:
    input:
        f"{OUTDIR}/phased/{{chrom}}.vcf.gz",
    output:
        multiext(f"{OUTDIR}/pgen/{{chrom}}", ".pgen", ".pvar", ".psam"),
    params:
        prefix=lambda w: f"{OUTDIR}/pgen/{w.chrom}",
    log:
        f"{LOGDIR}/export_pgen/{{chrom}}.txt",
    benchmark:
        f"{BENCHDIR}/export_pgen/{{chrom}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # plink2 otherwise takes half the node's memory, past the job's limit.
        plink2 \
            --vcf {input:q} \
            --make-pgen \
            --threads {threads} \
            --memory $(( {resources.mem_mb:q} * 4 / 5 )) \
            --out {params.prefix:q}

        sync {output:q}
        """
