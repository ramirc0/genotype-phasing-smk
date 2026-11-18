# Final phased VCFs and their pgen copies. The phaser's rule file writes phase/.


# plink2's X F-statistic is near 1 with one X copy (XY, XYY, X0) and near 0 with two
# or more (XX, XXY, XXX), which is what phasing needs. Samples it cannot call stay diploid.
# Split multiallelics leave half-calls such as 1/., which plink2 rejects by default.
rule infer_x_ploidy:
    input:
        bcf=f"{OUTDIR}/prepare/chrX_nonpar.bcf",
        csi=f"{OUTDIR}/prepare/chrX_nonpar.bcf.csi",
    output:
        sexcheck=f"{OUTDIR}/x_ploidy.tsv",
        haploids=f"{OUTDIR}/haploid_samples.txt",
        psam=f"{OUTDIR}/infer_x_ploidy/x.psam",
    params:
        prefix=f"{OUTDIR}/infer_x_ploidy/x",
    log:
        f"{LOGDIR}/infer_x_ploidy.txt",
    benchmark:
        f"{BENCHDIR}/infer_x_ploidy.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        mkdir -p "$(dirname {params.prefix:q})"

        plink2 \
            --bcf {input.bcf:q} \
            --vcf-half-call missing \
            --impute-sex max-female-xf=0.2 min-male-xf=0.8 \
            --make-just-psam \
            --threads {threads} \
            --memory $(( {resources.mem_mb:q} * 4 / 5 )) \
            --out {params.prefix:q}

        sed '1s/^#//' {params.prefix:q}.sexcheck > {output.sexcheck:q}
        awk -F'\t' 'NR == 1 {{ for (i = 1; i <= NF; i++) col[$i] = i; next }}
            $col["SNPSEX"] == 1 {{ print $col["IID"] }}' {output.sexcheck:q} > {output.haploids:q}
        echo "Samples with one X: $(wc -l < {output.haploids:q})"

        sync {output:q}
        """


# SHAPEIT5 --haploids means to set these hets missing, but its VAR_SET_MIS only ORs
# a bit, which turns a het into a scaffold genotype held at its input phase. Half-calls
# of split multiallelics go missing too, so SHAPEIT5 imputes both as homozygous.
rule mask_haploid_hets:
    input:
        bcf=f"{OUTDIR}/prepare/chrX_nonpar.bcf",
        csi=f"{OUTDIR}/prepare/chrX_nonpar.bcf.csi",
        haploids=f"{OUTDIR}/haploid_samples.txt",
    output:
        bcf=f"{OUTDIR}/prepare/chrX_nonpar.masked.bcf",
        csi=f"{OUTDIR}/prepare/chrX_nonpar.masked.bcf.csi",
    log:
        f"{LOGDIR}/mask_haploid_hets.txt",
    benchmark:
        f"{BENCHDIR}/mask_haploid_hets.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bcftools +setGT -Ob --output {output.bcf:q} --write-index=csi {input.bcf:q} -- \
            --target-gt q \
            --new-gt . \
            --include 'GT[@{input.haploids}]="het" | GT[@{input.haploids}]="mis"'

        sync {output:q}
        """


# A het left in a sample with one X means its non-PAR haplotypes disagree.
rule concat_chrx:
    input:
        regions=[f"{OUTDIR}/phase/{region}.bcf" for region in CHRX_REGIONS],
        nonpar=f"{OUTDIR}/phase/chrX_nonpar.bcf",
        haploids=f"{OUTDIR}/haploid_samples.txt",
    output:
        f"{OUTDIR}/phase/chrX.bcf",
    log:
        f"{LOGDIR}/concat_chrx.txt",
    benchmark:
        f"{BENCHDIR}/concat_chrx.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        if [[ -s {input.haploids:q} ]]; then
            hets=$(bcftools view --samples-file {input.haploids:q} -Ou {input.nonpar:q} \
                | bcftools query --include 'GT="het"' -f '[%SAMPLE\n]' | wc -l)
            if (( hets > 0 )); then
                echo "Non-PAR hets in samples with one X: $hets" >&2
                exit 1
            fi
        fi

        bcftools concat -Ob --output {output:q} {input.regions:q}

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


rule gather_phased:
    input:
        [f"{OUTDIR}/phased/{chrom}.vcf.gz" for chrom in GENOME_ORDER],
    output:
        vcf=f"{OUTDIR}/phased/genome.vcf.gz",
        tbi=f"{OUTDIR}/phased/genome.vcf.gz.tbi",
    log:
        f"{LOGDIR}/gather_phased.txt",
    benchmark:
        f"{BENCHDIR}/gather_phased.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bcftools concat \
            --threads {threads} \
            -Oz \
            --output {output.vcf:q} \
            --write-index=tbi \
            {input:q}

        echo "Phased records: $(bcftools index --nrecords {output.vcf:q})"

        sync {output:q}
        """


rule export_pgen:
    input:
        vcf=f"{OUTDIR}/phased/{{chrom}}.vcf.gz",
        sex=lambda w: [f"{OUTDIR}/infer_x_ploidy/x.psam"] if w.chrom == "chrX" else [],
    output:
        multiext(f"{OUTDIR}/pgen/{{chrom}}", ".pgen", ".pvar", ".psam"),
    params:
        prefix=lambda w: f"{OUTDIR}/pgen/{w.chrom}",
        # plink2 needs sex to load chrX. --split-par keeps PAR calls diploid in every sample.
        chrx=lambda w, input: (
            ["--split-par", "hg38", "--update-sex", input.sex[0]] if input.sex else []
        ),
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
            --vcf {input.vcf:q} \
            {params.chrx:q} \
            --make-pgen \
            --threads {threads} \
            --memory $(( {resources.mem_mb:q} * 4 / 5 )) \
            --out {params.prefix:q}

        sync {output:q}
        """
