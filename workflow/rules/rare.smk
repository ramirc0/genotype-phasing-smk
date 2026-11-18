# Rare-variant phasing with SHAPEIT5 phase_rare, on the phaser's scaffold of sites
# above filters.min_maf. The SHAPEIT5 docs recommend it above 2,000 samples, and
# phase_common below. See docs/PIPELINE.md.


def phaser_input(region):
    """Return the prepared BCF that SHAPEIT5 phases for a region."""
    masked = ".masked" if region == "chrX_nonpar" else ""
    return f"{OUTDIR}/prepare/{region}{masked}.bcf"


# One job per region, run in chunks of chunk_sites prepared sites, each with buffer_bp
# of scaffold on both sides. Chunks never split one position, because split
# multiallelics share it and would otherwise land in two chunks.
rule phase_rare:
    input:
        bcf=lambda w: phaser_input(w.region),
        csi=lambda w: f"{phaser_input(w.region)}.csi",
        scaffold=f"{SCAFFOLD_DIR}/{{region}}.bcf",
        map=config["genetic_maps"]["shapeit5"],
        haploids=lambda w: (
            [f"{OUTDIR}/haploid_samples.txt"] if w.region == "chrX_nonpar" else []
        ),
    output:
        f"{OUTDIR}/phase/{{region}}.bcf",
    params:
        chunk_sites=RARE["chunk_sites"],
        buffer=RARE["buffer_bp"],
        haploids=lambda w, input: (
            ["--haploids", input.haploids[0]] if input.haploids else []
        ),
    log:
        f"{LOGDIR}/phase_rare/{{region}}.txt",
    benchmark:
        f"{BENCHDIR}/phase_rare/{{region}}.tsv"
    software:
        SHAPEIT5_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        WORKDIR=$(mktemp -d)
        trap 'rm -rf "$WORKDIR"' EXIT

        bcftools query -f '%CHROM\t%POS\n' {input.bcf:q} \
        | awk -v size={params.chunk_sites:q} -v buffer={params.buffer:q} -v OFS='\t' '
            {{ chrom = $1; pos[NR] = $2 }}
            END {{
                first = 1
                while (first <= NR) {{
                    last = first + size - 1
                    if (last > NR) last = NR
                    while (last < NR && pos[last + 1] == pos[last]) last++
                    from = pos[first] - buffer; if (from < 1) from = 1
                    print chrom ":" pos[first] "-" pos[last], chrom ":" from "-" (pos[last] + buffer)
                    first = last + 1
                }}
            }}' > "$WORKDIR/chunks.tsv"
        echo "Chunks: $(wc -l < "$WORKDIR/chunks.tsv")"

        # phase_rare exits with an error on a chunk whose sites are all in the scaffold.
        chunk=0
        while IFS=$'\t' read -r input_region scaffold_region; do
            sites=$(bcftools view -H -r "$input_region" {input.bcf:q} | wc -l)
            phased=$(bcftools view -H -r "$input_region" {input.scaffold:q} | wc -l)
            if (( sites > phased )); then
                SHAPEIT5_phase_rare \
                    --input {input.bcf:q} \
                    --input-region "$input_region" \
                    --scaffold {input.scaffold:q} \
                    --scaffold-region "$scaffold_region" \
                    --map {input.map:q} \
                    {params.haploids:q} \
                    --score-singletons \
                    --thread {threads} \
                    --output "$WORKDIR/$chunk.bcf"
            else
                bcftools view -r "$input_region" -Ob --output "$WORKDIR/$chunk.bcf" {input.scaffold:q}
            fi
            echo "$WORKDIR/$chunk.bcf" >> "$WORKDIR/chunks.txt"
            chunk=$((chunk + 1))
        done < "$WORKDIR/chunks.tsv"

        bcftools concat \
            --file-list "$WORKDIR/chunks.txt" \
            -Ob \
            --output {output:q}

        sync {output:q}
        """
