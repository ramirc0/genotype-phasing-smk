# Array-site genotypes, relatedness, the unrelated set and genetic PCs.


rule make_array_sites:
    input:
        config["array_manifest"],
    output:
        f"{REFDIR}/array_sites.tsv",
    log:
        f"{LOGDIR}/make_array_sites.txt",
    benchmark:
        f"{BENCHDIR}/make_array_sites.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/array_sites.py \
            --manifest {input:q} \
            --output {output:q}

        sync {output:q}
        """


# Local ancestry runs on array SNPs thinned from the WGS genotypes. A split
# multiallelic site leaves several SNP records at one position, and RFMix2 matches
# records by position alone, so only positions with one SNP record are kept.
rule array_genotypes:
    input:
        vcf=f"{OUTDIR}/phased/{{chrom}}.vcf.gz",
        tbi=f"{OUTDIR}/phased/{{chrom}}.vcf.gz.tbi",
        sites=f"{REFDIR}/array_sites.tsv",
    output:
        bcf=f"{OUTDIR}/array/{{chrom}}.bcf",
        csi=f"{OUTDIR}/array/{{chrom}}.bcf.csi",
    log:
        f"{LOGDIR}/array_genotypes/{{chrom}}.txt",
    benchmark:
        f"{BENCHDIR}/array_genotypes/{{chrom}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bcftools view \
            --targets-file <(bcftools query \
                --targets-file {input.sites:q} \
                --include 'TYPE="snp"' \
                --format '%CHROM\t%POS\n' \
                {input.vcf:q} \
                | uniq -u) \
            --types snps \
            -Ob \
            --output {output.bcf:q} \
            --write-index=csi \
            {input.vcf:q}

        echo "Array records: $(bcftools index --nrecords {output.bcf:q})"

        sync {output:q}
        """


rule cohort_array_pgen:
    input:
        expand(f"{OUTDIR}/array/{{chrom}}.bcf", chrom=PHASED_AUTOSOMES),
    output:
        multiext(f"{OUTDIR}/array/cohort", ".pgen", ".pvar", ".psam"),
    params:
        prefix=f"{OUTDIR}/array/cohort",
    log:
        f"{LOGDIR}/cohort_array_pgen.txt",
    benchmark:
        f"{BENCHDIR}/cohort_array_pgen.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        WORKDIR=$(mktemp -d)
        trap 'rm -rf "$WORKDIR"' EXIT

        bcftools concat \
            --threads {threads} \
            -Ob \
            --output "$WORKDIR/cohort.bcf" \
            {input:q}

        plink2 \
            --bcf "$WORKDIR/cohort.bcf" \
            --make-pgen \
            --threads {threads} \
            --memory $(( {resources.mem_mb:q} * 4 / 5 )) \
            --out {params.prefix:q}

        sync {output:q}
        """


# KING-robust kinship tolerates population structure. 0.0442, the default cutoff,
# is the lower kinship bound of third-degree relatives.
rule relatedness:
    input:
        multiext(f"{OUTDIR}/array/cohort", ".pgen", ".pvar", ".psam"),
    output:
        kin0=f"{OUTDIR}/king.kin0",
        unrelated=f"{OUTDIR}/unrelated.txt",
    params:
        cohort=f"{OUTDIR}/array/cohort",
        prefix=f"{OUTDIR}/relatedness/king",
        cutoff=lambda w: config["kinship_cutoff"],
    log:
        f"{LOGDIR}/relatedness.txt",
    benchmark:
        f"{BENCHDIR}/relatedness.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        mkdir -p "$(dirname {params.prefix:q})"

        plink2 \
            --pfile {params.cohort:q} \
            --make-king-table \
            --king-table-filter {params.cutoff:q} \
            --threads {threads} \
            --memory $(( {resources.mem_mb:q} * 4 / 5 )) \
            --out {params.prefix:q}
        plink2 \
            --pfile {params.cohort:q} \
            --king-cutoff {params.cutoff:q} \
            --threads {threads} \
            --memory $(( {resources.mem_mb:q} * 4 / 5 )) \
            --out {params.prefix:q}

        mv {params.prefix:q}.kin0 {output.kin0:q}
        grep -v '^#' {params.prefix:q}.king.cutoff.in.id | cut -f1 > {output.unrelated:q}
        echo "Unrelated samples: $(wc -l < {output.unrelated:q})"

        sync {output:q}
        """


# PCs come from unrelated samples on pruned common SNPs outside long-range LD regions.
# Every sample is then projected, so relatives get PCs on the same scale.
rule genetic_pcs:
    input:
        cohort=multiext(f"{OUTDIR}/array/cohort", ".pgen", ".pvar", ".psam"),
        unrelated=f"{OUTDIR}/unrelated.txt",
        long_range_ld=config["long_range_ld"],
    output:
        pcs=f"{OUTDIR}/pcs.tsv",
        eigenval=f"{OUTDIR}/pcs.eigenval",
    params:
        cohort=f"{OUTDIR}/array/cohort",
        prefix=f"{OUTDIR}/genetic_pcs/pca",
        num_pcs=lambda w: config["num_pcs"],
    log:
        f"{LOGDIR}/genetic_pcs.txt",
    benchmark:
        f"{BENCHDIR}/genetic_pcs.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        mkdir -p "$(dirname {params.prefix:q})"
        memory=$(( {resources.mem_mb:q} * 4 / 5 ))

        plink2 \
            --pfile {params.cohort:q} \
            --keep {input.unrelated:q} \
            --maf 0.05 \
            --exclude bed0 {input.long_range_ld:q} \
            --indep-pairwise 1000kb 1 0.1 \
            --threads {threads} \
            --memory "$memory" \
            --out {params.prefix:q}
        plink2 \
            --pfile {params.cohort:q} \
            --keep {input.unrelated:q} \
            --extract {params.prefix:q}.prune.in \
            --freq counts \
            --pca {params.num_pcs:q} allele-wts \
            --threads {threads} \
            --memory "$memory" \
            --out {params.prefix:q}
        # plink2 adds a PROVISIONAL_REF? column only when some REF is provisional, so find A1.
        a1=$(head -1 {params.prefix:q}.eigenvec.allele | tr '\t' '\n' | grep -nx A1 | cut -d: -f1)
        plink2 \
            --pfile {params.cohort:q} \
            --extract {params.prefix:q}.prune.in \
            --read-freq {params.prefix:q}.acount \
            --score {params.prefix:q}.eigenvec.allele 2 "$a1" header-read no-mean-imputation variance-standardize \
            --score-col-nums "$(( a1 + 1 ))-$(( a1 + {params.num_pcs:q} ))" \
            --threads {threads} \
            --memory "$memory" \
            --out {params.prefix:q}.projected

        cut -f1,4- {params.prefix:q}.projected.sscore \
            | sed -e '1s/^#//' -e '1s/_AVG//g' > {output.pcs:q}
        cp {params.prefix:q}.eigenval {output.eigenval:q}

        sync {output:q}
        """
