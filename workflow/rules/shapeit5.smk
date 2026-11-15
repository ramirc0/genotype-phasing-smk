# Phasing with SHAPEIT5.


# The panel and the cohort condition each other. The SHAPEIT5 docs skip phase_rare
# below 2,000 samples. Sites absent from the panel are dropped, and the log counts them.
# Outside the PARs, SHAPEIT5 sets hets of samples with one X missing and imputes them
# as homozygous.
rule phase_shapeit5:
    input:
        bcf=f"{OUTDIR}/prepare/{{region}}.bcf",
        csi=f"{OUTDIR}/prepare/{{region}}.bcf.csi",
        panel=PANEL,
        panel_index=f"{PANEL}.csi",
        map=config["genetic_maps"]["shapeit5"],
        haploids=lambda w: (
            [f"{OUTDIR}/haploid_samples.txt"] if w.region == "chrX_nonpar" else []
        ),
    output:
        f"{OUTDIR}/phase/{{region}}.bcf",
    params:
        region=lambda w: REGION_SPAN[w.region],
        haploids=lambda w, input: (
            ["--haploids", input.haploids[0]] if input.haploids else []
        ),
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
            {params.haploids:q} \
            --thread {threads} \
            --output {output:q}

        sync {output:q}
        """
