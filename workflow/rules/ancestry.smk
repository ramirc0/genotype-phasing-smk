# Reference selection, ADMIXTURE, and the local ancestry review. The local
# ancestry tool's rule file writes local_ancestry/.


rule reference_candidates:
    input:
        metadata=config["panel_metadata"],
        panel=PANEL.format(region=PHASED_AUTOSOMES[0]),
    output:
        f"{OUTDIR}/ancestry/candidates.tsv",
    params:
        populations=lambda w: [
            f"{ancestry}={','.join(populations)}"
            for ancestry, populations in LAI["populations"].items()
        ],
    log:
        f"{LOGDIR}/reference_candidates.txt",
    benchmark:
        f"{BENCHDIR}/reference_candidates.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        WORKDIR=$(mktemp -d)
        trap 'rm -rf "$WORKDIR"' EXIT

        # polars memory-maps its input, which a process substitution cannot offer.
        bcftools query --list-samples {input.panel:q} > "$WORKDIR/panel_samples.txt"
        python {workflow.basedir:q}/scripts/reference_candidates.py \
            --metadata {input.metadata:q} \
            --panel-samples "$WORKDIR/panel_samples.txt" \
            --populations {params.populations:q} \
            --output {output:q}

        sync {output:q}
        """


# RFMix2 matches query and reference records by position alone, so both files keep
# only the array SNPs whose REF and ALT agree exactly.
rule match_array_sites:
    input:
        query=f"{OUTDIR}/array/{{chrom}}.bcf",
        query_index=f"{OUTDIR}/array/{{chrom}}.bcf.csi",
        panel=lambda w: PANEL.format(region=w.chrom),
        panel_index=lambda w: f"{PANEL.format(region=w.chrom)}.csi",
        candidates=f"{OUTDIR}/ancestry/candidates.tsv",
        sites=f"{REFDIR}/array_sites.tsv",
    output:
        query=f"{OUTDIR}/ancestry/{{chrom}}.query.bcf",
        query_index=f"{OUTDIR}/ancestry/{{chrom}}.query.bcf.csi",
        reference=f"{OUTDIR}/ancestry/{{chrom}}.reference.bcf",
        reference_index=f"{OUTDIR}/ancestry/{{chrom}}.reference.bcf.csi",
    log:
        f"{LOGDIR}/match_array_sites/{{chrom}}.txt",
    benchmark:
        f"{BENCHDIR}/match_array_sites/{{chrom}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        WORKDIR=$(mktemp -d)
        trap 'rm -rf "$WORKDIR"' EXIT

        bcftools view \
            --samples-file <(tail -n +2 {input.candidates:q} | cut -f1) \
            --targets-file {input.sites:q} \
            --types snps \
            -Ou \
            {input.panel:q} \
        | bcftools annotate \
            --remove '^FORMAT/GT' \
            -Ob \
            --output "$WORKDIR/panel.bcf" \
            --write-index=csi
        bcftools isec \
            -n=2 \
            --collapse none \
            -Ob \
            --prefix "$WORKDIR/isec" \
            {input.query:q} \
            "$WORKDIR/panel.bcf"
        mv "$WORKDIR/isec/0000.bcf" {output.query:q}
        mv "$WORKDIR/isec/0001.bcf" {output.reference:q}
        bcftools index --csi {output.query:q}
        bcftools index --csi {output.reference:q}

        echo "Matched array SNPs: $(bcftools index --nrecords {output.query:q})"

        sync {output:q}
        """


rule make_lai_maps:
    input:
        lambda w: config["genetic_maps"]["shapeit5"].format(region=w.chrom),
    output:
        f"{REFDIR}/lai_maps/{{chrom}}.tsv",
    log:
        f"{LOGDIR}/make_lai_maps/{{chrom}}.txt",
    benchmark:
        f"{BENCHDIR}/make_lai_maps/{{chrom}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        zcat {input:q} \
            | awk -v chrom={wildcards.chrom:q} -v OFS='\t' 'NR > 1 {{ print chrom, $1, $3 }}' \
            > {output:q}

        sync {output:q}
        """


# Unsupervised, so each candidate's purity can be read off its components. Only
# unrelated cohort samples enter, so families cannot pull a component their way.
rule admixture:
    input:
        query=expand(f"{OUTDIR}/ancestry/{{chrom}}.query.bcf", chrom=PHASED_AUTOSOMES),
        reference=expand(
            f"{OUTDIR}/ancestry/{{chrom}}.reference.bcf", chrom=PHASED_AUTOSOMES
        ),
        unrelated=f"{OUTDIR}/unrelated.txt",
        long_range_ld=config["long_range_ld"],
    output:
        plink=multiext(f"{OUTDIR}/admixture/merged", ".bed", ".bim", ".fam"),
        q=f"{OUTDIR}/admixture/merged.{len(LAI['populations'])}.Q",
        p=f"{OUTDIR}/admixture/merged.{len(LAI['populations'])}.P",
    params:
        prefix=f"{OUTDIR}/admixture/merged",
        components=len(LAI["populations"]),
        seed=lambda w: config["seed"],
    log:
        f"{LOGDIR}/admixture.txt",
    benchmark:
        f"{BENCHDIR}/admixture.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        WORKDIR=$(mktemp -d)
        trap 'rm -rf "$WORKDIR"' EXIT
        memory=$(( {resources.mem_mb:q} * 4 / 5 ))

        bcftools concat --threads {threads} -Ob --output "$WORKDIR/query.bcf" --write-index=csi {input.query:q}
        bcftools concat --threads {threads} -Ob --output "$WORKDIR/reference.bcf" --write-index=csi {input.reference:q}
        bcftools merge --threads {threads} -Ob --output "$WORKDIR/merged.bcf" "$WORKDIR/query.bcf" "$WORKDIR/reference.bcf"
        bcftools query --list-samples "$WORKDIR/reference.bcf" | cat {input.unrelated:q} - > "$WORKDIR/keep.txt"

        plink2 \
            --bcf "$WORKDIR/merged.bcf" \
            --keep "$WORKDIR/keep.txt" \
            --maf 0.05 \
            --exclude bed0 {input.long_range_ld:q} \
            --indep-pairwise 1000kb 1 0.1 \
            --threads {threads} \
            --memory "$memory" \
            --out "$WORKDIR/prune"
        plink2 \
            --bcf "$WORKDIR/merged.bcf" \
            --keep "$WORKDIR/keep.txt" \
            --extract "$WORKDIR/prune.prune.in" \
            --make-bed \
            --threads {threads} \
            --memory "$memory" \
            --out {params.prefix:q}

        # ADMIXTURE writes its .Q and .P files to the working directory.
        (cd "$(dirname {params.prefix:q})" \
            && admixture -s {params.seed:q} -j{threads} "$(basename {params.prefix:q}).bed" {params.components:q})

        sync {output:q}
        """


rule select_reference:
    input:
        fam=f"{OUTDIR}/admixture/merged.fam",
        q=f"{OUTDIR}/admixture/merged.{len(LAI['populations'])}.Q",
        candidates=f"{OUTDIR}/ancestry/candidates.tsv",
    output:
        sample_map=f"{OUTDIR}/lai_sample_map.tsv",
        global_ancestry=f"{OUTDIR}/global_ancestry.tsv",
    params:
        min_purity=lambda w: [
            f"{ancestry}={LAI['min_purity'][ancestry]}" for ancestry in LAI["populations"]
        ],
    log:
        f"{LOGDIR}/select_reference.txt",
    benchmark:
        f"{BENCHDIR}/select_reference.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/select_reference.py \
            --fam {input.fam:q} \
            --q {input.q:q} \
            --candidates {input.candidates:q} \
            --min-purity {params.min_purity:q} \
            --sample-map {output.sample_map:q} \
            --global-ancestry {output.global_ancestry:q}

        sync {output:q}
        """


rule plot_local_ancestry:
    input:
        msp=expand(f"{OUTDIR}/local_ancestry/{{chrom}}.msp.tsv", chrom=PHASED_AUTOSOMES),
        admixture=f"{OUTDIR}/global_ancestry.tsv",
    output:
        global_ancestry=f"{OUTDIR}/local_ancestry_global.tsv",
        karyograms=f"{OUTDIR}/plots/karyograms.svg",
        karyograms_png=f"{OUTDIR}/plots/karyograms.png",
        ideograms=f"{OUTDIR}/plots/ideograms.pdf",
        comparison=f"{OUTDIR}/plots/global_ancestry.svg",
        comparison_png=f"{OUTDIR}/plots/global_ancestry.png",
    params:
        seed=lambda w: config["seed"],
        pdf_samples=lambda w: LAI["pdf_samples"],
        colors=lambda w: [f"{a}={c}" for a, c in LAI.get("colors", {}).items()],
    log:
        f"{LOGDIR}/plot_local_ancestry.txt",
    benchmark:
        f"{BENCHDIR}/plot_local_ancestry.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/plot_local_ancestry.py \
            --msp {input.msp:q} \
            --admixture {input.admixture:q} \
            --seed {params.seed:q} \
            --pdf-samples {params.pdf_samples:q} \
            --colors {params.colors:q} \
            --global-ancestry {output.global_ancestry:q} \
            --karyograms {output.karyograms:q} \
            --ideograms {output.ideograms:q} \
            --comparison {output.comparison:q}

        sync {output:q}
        """


# Tractor needs the phased VCF whose haplotypes local ancestry was called on. The
# local ancestry query comes from this VCF, so samples and haplotypes line up.
rule extract_tracts:
    input:
        vcf=f"{OUTDIR}/phased/{{chrom}}.vcf.gz",
        msp=f"{OUTDIR}/local_ancestry/{{chrom}}.msp.tsv",
    output:
        [
            f"{OUTDIR}/tracts/{{chrom}}.anc{i}.{kind}.txt.gz"
            for i in range(len(LAI["populations"]))
            for kind in ("dosage", "hapcount")
        ],
    params:
        outdir=f"{OUTDIR}/tracts",
        ancestries=len(LAI["populations"]),
        common=f"MAF>={FILTERS['min_maf']}" if RARE["enabled"] else "",
    log:
        f"{LOGDIR}/extract_tracts/{{chrom}}.txt",
    benchmark:
        f"{BENCHDIR}/extract_tracts/{{chrom}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        WORKDIR=$(mktemp -d)
        trap 'rm -rf "$WORKDIR"' EXIT

        # Tractor names its outputs after the VCF, so a filtered copy keeps the name.
        VCF={input.vcf:q}
        COMMON={params.common:q}
        if [[ -n "$COMMON" ]]; then
            VCF="$WORKDIR/{wildcards.chrom}.vcf.gz"
            bcftools view --include "$COMMON" -Oz --output "$VCF" {input.vcf:q}
        fi

        python {workflow.basedir:q}/scripts/vendor/extract_tracts.py \
            --vcf "$VCF" \
            --msp {input.msp:q} \
            --num-ancs {params.ancestries:q} \
            --output-dir {params.outdir:q} \
            --compress-output

        sync {output:q}
        """
