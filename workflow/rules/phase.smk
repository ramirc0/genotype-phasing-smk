# Final phased VCFs and their pgen copies. The phaser's rule file writes phase/.


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
