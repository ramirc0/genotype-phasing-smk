# Phasing with Eagle2.


# Eagle's auto rule runs one iteration when the cohort is under half the panel, so each
# sample would be phased against the panel alone. Two iterations add cohort haplotypes.
rule phase_eagle2:
    input:
        bcf=f"{OUTDIR}/prepare/{{region}}.bcf",
        csi=f"{OUTDIR}/prepare/{{region}}.bcf.csi",
        panel=PANEL,
        panel_index=f"{PANEL}.csi",
        map=config["genetic_maps"]["eagle2"],
    output:
        f"{OUTDIR}/phase/{{region}}.bcf",
    params:
        prefix=lambda w: f"{OUTDIR}/phase/{w.region}",
    log:
        f"{LOGDIR}/phase_eagle2/{{region}}.txt",
    benchmark:
        f"{BENCHDIR}/phase_eagle2/{{region}}.tsv"
    software:
        EAGLE2_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        eagle \
            --vcfRef {input.panel:q} \
            --vcfTarget {input.bcf:q} \
            --geneticMapFile {input.map:q} \
            --pbwtIters 2 \
            --numThreads {threads} \
            --vcfOutFormat b \
            --outPrefix {params.prefix:q}

        sync {output:q}
        """
