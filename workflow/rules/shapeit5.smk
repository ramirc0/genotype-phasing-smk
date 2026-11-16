# Phasing with SHAPEIT5.


# The first pass phases against the panel, where panel and cohort condition each other,
# and leaves out sites the panel lacks. Those are real variants the panel's QC removed,
# so a second pass phases every site from the cohort alone, holding the first pass's
# phase fixed as a scaffold. The SHAPEIT5 docs skip phase_rare below 2,000 samples.
# Both rules below run it on the span in $REGION.
TWO_PASS = r"""
        WORKDIR=$(mktemp -d)
        trap 'rm -rf "$WORKDIR"' EXIT

        SHAPEIT5_phase_common \
            --input {input.bcf:q} \
            --reference {input.panel:q} \
            --map {input.map:q} \
            --region "$REGION" \
            {params.haploids:q} \
            --thread {threads} \
            --output "$WORKDIR/scaffold.bcf"
        SHAPEIT5_phase_common \
            --input {input.bcf:q} \
            --scaffold "$WORKDIR/scaffold.bcf" \
            --map {input.map:q} \
            --region "$REGION" \
            {params.haploids:q} \
            --thread {threads} \
            --output {output:q}

        sync {output:q}
"""


rule phase_shapeit5:
    input:
        bcf=f"{OUTDIR}/prepare/{{region}}.bcf",
        csi=f"{OUTDIR}/prepare/{{region}}.bcf.csi",
        panel=PANEL,
        panel_index=f"{PANEL}.csi",
        map=config["genetic_maps"]["shapeit5"],
    output:
        f"{OUTDIR}/phase/{{region}}.bcf",
    wildcard_constraints:
        region="|".join(re.escape(r) for r in REGION_SPAN if r != "chrX_nonpar"),
    params:
        region=lambda w: REGION_SPAN[w.region],
        haploids=[],
    log:
        f"{LOGDIR}/phase_shapeit5/{{region}}.txt",
    benchmark:
        f"{BENCHDIR}/phase_shapeit5/{{region}}.tsv"
    software:
        SHAPEIT5_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        REGION={params.region:q}
        """ + TWO_PASS


# With haploid samples, SHAPEIT5's memory on the whole non-PAR grows faster than the
# span: past 243 GB, against 16 GB for chr1. The cohort's haploid samples drive it.
# Chunks of equal site counts overlap by overlap_sites, so ligation always has
# shared sites to align phase on.
rule chrx_chunks:
    input:
        f"{OUTDIR}/prepare/chrX_nonpar.bcf",
    output:
        f"{OUTDIR}/phase/chrX_nonpar/chunks.txt",
    params:
        chunks=CHRX_CHUNKS["count"],
        overlap=CHRX_CHUNKS["overlap_sites"],
    log:
        f"{LOGDIR}/chrx_chunks.txt",
    benchmark:
        f"{BENCHDIR}/chrx_chunks.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bcftools query -f '%POS\n' {input:q} \
        | awk -v count={params.chunks:q} -v overlap={params.overlap:q} '
            {{ pos[NR - 1] = $1 }}
            END {{
                size = int((NR + count - 1) / count)
                for (i = 0; i < count; i++) {{
                    lo = i * size - overlap; if (lo < 0) lo = 0
                    hi = (i + 1) * size - 1 + overlap; if (hi > NR - 1) hi = NR - 1
                    print "chrX:" pos[lo] "-" pos[hi]
                }}
            }}' > {output:q}
        cat {output:q}

        sync {output:q}
        """


# SHAPEIT5 sets hets of samples with one X missing and imputes them as homozygous.
rule phase_chrx_chunk:
    input:
        chunks=f"{OUTDIR}/phase/chrX_nonpar/chunks.txt",
        bcf=f"{OUTDIR}/prepare/chrX_nonpar.bcf",
        csi=f"{OUTDIR}/prepare/chrX_nonpar.bcf.csi",
        panel=PANEL.format(region="chrX_nonpar"),
        panel_index=f"{PANEL.format(region='chrX_nonpar')}.csi",
        map=config["genetic_maps"]["shapeit5"].format(region="chrX_nonpar"),
        haploids=f"{OUTDIR}/haploid_samples.txt",
    output:
        f"{OUTDIR}/phase/chrX_nonpar/{{chunk}}.bcf",
    wildcard_constraints:
        chunk=r"\d+",
    params:
        line=lambda w: int(w.chunk) + 1,
        haploids=lambda w, input: ["--haploids", input.haploids],
    log:
        f"{LOGDIR}/phase_chrx_chunk/{{chunk}}.txt",
    benchmark:
        f"{BENCHDIR}/phase_chrx_chunk/{{chunk}}.tsv"
    software:
        SHAPEIT5_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        REGION=$(sed -n {params.line:q}p {input.chunks:q})
        """ + TWO_PASS


rule ligate_chrx:
    input:
        expand(
            f"{OUTDIR}/phase/chrX_nonpar/{{chunk}}.bcf",
            chunk=range(CHRX_CHUNKS["count"]),
        ),
    output:
        f"{OUTDIR}/phase/chrX_nonpar.bcf",
    log:
        f"{LOGDIR}/ligate_chrx.txt",
    benchmark:
        f"{BENCHDIR}/ligate_chrx.tsv"
    software:
        SHAPEIT5_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        WORKDIR=$(mktemp -d)
        trap 'rm -rf "$WORKDIR"' EXIT

        printf '%s\n' {input:q} > "$WORKDIR/chunks.txt"
        SHAPEIT5_ligate \
            --input "$WORKDIR/chunks.txt" \
            --thread {threads} \
            --output {output:q}

        sync {output:q}
        """
