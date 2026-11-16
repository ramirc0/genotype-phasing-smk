# Local ancestry with FLARE.


# FLARE numbers ancestries in the order the reference-panel file first names them, so
# the file lists samples in local_ancestry.populations order. It reads VCF, not BCF.
# The query holds array SNPs, hence array=true. em=true fits the generations since
# admixture, so local_ancestry.rfmix_generations does not apply.
rule local_ancestry_flare:
    input:
        query=f"{OUTDIR}/ancestry/{{chrom}}.query.bcf",
        query_index=f"{OUTDIR}/ancestry/{{chrom}}.query.bcf.csi",
        reference=f"{OUTDIR}/ancestry/{{chrom}}.reference.bcf",
        reference_index=f"{OUTDIR}/ancestry/{{chrom}}.reference.bcf.csi",
        sample_map=f"{OUTDIR}/lai_sample_map.tsv",
        genetic_map=f"{REFDIR}/lai_maps/{{chrom}}.tsv",
    output:
        f"{OUTDIR}/local_ancestry/{{chrom}}.anc.vcf.gz",
    params:
        ancestries=list(LAI["populations"]),
        seed=lambda w: config["seed"],
    log:
        f"{LOGDIR}/local_ancestry_flare/{{chrom}}.txt",
    benchmark:
        f"{BENCHDIR}/local_ancestry_flare/{{chrom}}.tsv"
    software:
        FLARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        WORKDIR=$(mktemp -d)
        trap 'rm -rf "$WORKDIR"' EXIT

        for ancestry in {params.ancestries:q}; do
            awk -v a="$ancestry" '$2 == a' {input.sample_map:q}
        done > "$WORKDIR/panel.txt"
        bcftools view \
            --samples-file <(cut -f1 "$WORKDIR/panel.txt") \
            -Oz \
            --output "$WORKDIR/reference.vcf.gz" \
            {input.reference:q}
        bcftools view -Oz --output "$WORKDIR/query.vcf.gz" {input.query:q}
        awk -v OFS='\t' '{{ print $1, ".", $3, $2 }}' {input.genetic_map:q} > "$WORKDIR/map.txt"

        # FLARE's JVM otherwise sizes its heap from the node's memory, past the job's limit.
        JAVA_OPTS="-Xmx$(( {resources.mem_mb:q} * 4 / 5 ))m" flare \
            ref="$WORKDIR/reference.vcf.gz" \
            ref-panel="$WORKDIR/panel.txt" \
            gt="$WORKDIR/query.vcf.gz" \
            map="$WORKDIR/map.txt" \
            out="$WORKDIR/flare" \
            array=true \
            probs=true \
            seed={params.seed:q} \
            nthreads={threads}
        cp "$WORKDIR/flare.anc.vcf.gz" {output:q}

        sync {output:q}
        """


rule flare_msp:
    input:
        anc=f"{OUTDIR}/local_ancestry/{{chrom}}.anc.vcf.gz",
        genetic_map=f"{REFDIR}/lai_maps/{{chrom}}.tsv",
    output:
        f"{OUTDIR}/local_ancestry/{{chrom}}.msp.tsv",
    params:
        ancestries=list(LAI["populations"]),
    log:
        f"{LOGDIR}/flare_msp/{{chrom}}.txt",
    benchmark:
        f"{BENCHDIR}/flare_msp/{{chrom}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        WORKDIR=$(mktemp -d)
        trap 'rm -rf "$WORKDIR"' EXIT

        bcftools query -l {input.anc:q} > "$WORKDIR/samples.txt"
        bcftools view -h {input.anc:q} | grep '^##ANCESTRY=' > "$WORKDIR/header.txt"
        bcftools query -f '%CHROM\t%POS[\t%AN1\t%AN2]\n' {input.anc:q} > "$WORKDIR/calls.tsv"

        python {workflow.basedir:q}/scripts/flare_msp.py \
            --calls "$WORKDIR/calls.tsv" \
            --samples "$WORKDIR/samples.txt" \
            --header "$WORKDIR/header.txt" \
            --ancestries {params.ancestries:q} \
            --map {input.genetic_map:q} \
            --output {output:q}

        sync {output:q}
        """
