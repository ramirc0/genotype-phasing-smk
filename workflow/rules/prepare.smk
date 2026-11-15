# Per-region genotype filtering before phasing.


# Multiallelics split with the other ALT set missing, so no sample gains a REF allele it lacks.
# Genotypes below the depth cutoff go missing, and the phaser imputes them.
rule prepare_genotypes:
    input:
        vcf=CALLSET,
        tbi=f"{CALLSET}.tbi",
        exclude=config["exclude_samples"],
        fasta=FASTA,
        fasta_indexes=fasta_indexes(FASTA),
    output:
        bcf=f"{OUTDIR}/prepare/{{region}}.bcf",
        csi=f"{OUTDIR}/prepare/{{region}}.bcf.csi",
    params:
        region=lambda w: REGION_SPAN[w.region],
        depth_filter=lambda w: f"FMT/DP<{FILTERS['min_dp']}",
        site_filter=lambda w: (
            f"MAF>={FILTERS['min_maf']} && F_MISSING<={FILTERS['max_missing']}"
        ),
    log:
        f"{LOGDIR}/prepare_genotypes/{{region}}.txt",
    benchmark:
        f"{BENCHDIR}/prepare_genotypes/{{region}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bcftools view \
            --regions {params.region:q} \
            --apply-filters PASS \
            --samples-file ^<(grep -v '^#' {input.exclude:q} | cut -f1) \
            -Ou \
            {input.vcf:q} \
        | bcftools annotate --remove 'INFO,^FORMAT/GT,FORMAT/DP' -Ou \
        | bcftools norm \
            --multiallelics -any \
            --multi-overlaps . \
            --fasta-ref {input.fasta:q} \
            -Ou \
        | bcftools view --exclude 'ALT="*"' -Ou \
        | bcftools +setGT -Ou -- \
            --target-gt q \
            --new-gt . \
            --include {params.depth_filter:q} \
        | bcftools annotate --remove FORMAT/DP -Ou \
        | bcftools +fill-tags -Ou -- --tags AN,AC,AF,MAF,F_MISSING \
        | bcftools view \
            --include {params.site_filter:q} \
            --threads {threads} \
            -Ob \
            --output {output.bcf:q} \
            --write-index=csi

        echo "Prepared records: $(bcftools index --nrecords {output.bcf:q})"

        sync {output:q}
        """
