# Local ancestry with gnomix.


# gnomix trains on admixed samples it simulates from the reference. Phase correction
# stays off, so its haplotypes stay those of the phased VCF that Tractor reads.
rule local_ancestry_gnomix:
    input:
        query=f"{OUTDIR}/ancestry/{{chrom}}.query.bcf",
        query_index=f"{OUTDIR}/ancestry/{{chrom}}.query.bcf.csi",
        reference=f"{OUTDIR}/ancestry/{{chrom}}.reference.bcf",
        reference_index=f"{OUTDIR}/ancestry/{{chrom}}.reference.bcf.csi",
        sample_map=f"{OUTDIR}/lai_sample_map.tsv",
        genetic_map=f"{REFDIR}/lai_maps/{{chrom}}.tsv",
        config=[LAI["gnomix_config"]] if "gnomix_config" in LAI else [],
    output:
        msp=f"{OUTDIR}/local_ancestry/{{chrom}}.msp.tsv",
        fb=f"{OUTDIR}/local_ancestry/{{chrom}}.fb.tsv",
    params:
        outdir=lambda w: f"{OUTDIR}/gnomix/{w.chrom}",
        seed=lambda w: config["seed"],
    log:
        f"{LOGDIR}/local_ancestry_gnomix/{{chrom}}.txt",
    benchmark:
        f"{BENCHDIR}/local_ancestry_gnomix/{{chrom}}.tsv"
    software:
        GNOMIX_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        WORKDIR=$(mktemp -d)
        trap 'rm -rf "$WORKDIR"' EXIT

        # gnomix reads VCF only.
        bcftools view -Oz --output "$WORKDIR/query.vcf.gz" {input.query:q}
        bcftools view \
            --samples-file <(cut -f1 {input.sample_map:q}) \
            -Oz \
            --output "$WORKDIR/reference.vcf.gz" \
            {input.reference:q}

        base=$(printf '%s' {input.config:q})
        [[ -n $base ]] || base=$CONDA_PREFIX/share/gnomix/configs/config_array.yaml
        python {workflow.basedir:q}/scripts/gnomix_config.py \
            --base "$base" \
            --seed {params.seed:q} \
            --threads {threads} \
            --output "$WORKDIR/config.yaml"

        gnomix \
            "$WORKDIR/query.vcf.gz" \
            {params.outdir:q} \
            {wildcards.chrom:q} \
            False \
            {input.genetic_map:q} \
            "$WORKDIR/reference.vcf.gz" \
            {input.sample_map:q} \
            "$WORKDIR/config.yaml"

        cp {params.outdir:q}/query_results.msp {output.msp:q}
        cp {params.outdir:q}/query_results.fb {output.fb:q}

        sync {output:q}
        """
