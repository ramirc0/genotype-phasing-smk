# Local ancestry with RFMix2.


rule local_ancestry_rfmix2:
    input:
        query=f"{OUTDIR}/ancestry/{{chrom}}.query.bcf",
        query_index=f"{OUTDIR}/ancestry/{{chrom}}.query.bcf.csi",
        reference=f"{OUTDIR}/ancestry/{{chrom}}.reference.bcf",
        reference_index=f"{OUTDIR}/ancestry/{{chrom}}.reference.bcf.csi",
        sample_map=f"{OUTDIR}/lai_sample_map.tsv",
        genetic_map=f"{REFDIR}/lai_maps/{{chrom}}.tsv",
    output:
        msp=f"{OUTDIR}/local_ancestry/{{chrom}}.msp.tsv",
        fb=f"{OUTDIR}/local_ancestry/{{chrom}}.fb.tsv",
    params:
        prefix=lambda w: f"{OUTDIR}/local_ancestry/{w.chrom}",
        generations=lambda w: LAI["rfmix_generations"],
        seed=lambda w: config["seed"],
    log:
        f"{LOGDIR}/local_ancestry_rfmix2/{{chrom}}.txt",
    benchmark:
        f"{BENCHDIR}/local_ancestry_rfmix2/{{chrom}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        WORKDIR=$(mktemp -d)
        trap 'rm -rf "$WORKDIR"' EXIT

        bcftools view \
            --samples-file <(cut -f1 {input.sample_map:q}) \
            -Ob \
            --output "$WORKDIR/reference.bcf" \
            --write-index=csi \
            {input.reference:q}

        # RFMix2 reads a long option's value only after '='. A separate value crashes it.
        rfmix \
            --query-file={input.query:q} \
            --reference-file="$WORKDIR/reference.bcf" \
            --sample-map={input.sample_map:q} \
            --genetic-map={input.genetic_map:q} \
            --output-basename={params.prefix:q} \
            --chromosome={wildcards.chrom:q} \
            --generations={params.generations:q} \
            --random-seed={params.seed:q} \
            --n-threads={threads}

        sync {output:q}
        """
