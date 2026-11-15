#!/usr/bin/env bash
# Build the toy inputs in resources/ from 1000 Genomes and HGDP+1kGP data.
# The test does not run this. It records how the inputs were made.
set -euo pipefail
cd "$(dirname "$0")"

jg_test=/home/ramirezc/Projects/joint-genotyping-smk/main/.tests/integration
nygc_chrx=/data/zusers/ramirezc/static/hg38/variants/1000G-30x/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.recalibrated_variants.vcf.gz
grch38=/data/zusers/ramirezc/joint-genotyping-smk/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
res=/data/zusers/ramirezc/data/genotype-phasing-smk/resources

# The chrX windows straddle both PAR bounds, so each chrX region holds sites.
chrx_windows=chrX:2700000-2860000,chrX:155650000-155760000
declare -A span=(
    [chr20]=chr20:1-347695
    [chrX_par1]=chrX:2700000-2781479
    [chrX_nonpar]="chrX:2781480-2860000,chrX:155650000-155701382"
    [chrX_par2]=chrX:155701383-155760000
)

px() { pixi exec -c conda-forge -c bioconda -s bcftools==1.25 -s htslib==1.25 -s samtools==1.25 -s python==3.13 "$@"; }

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p resources/panel resources/maps

px bcftools query -l "$jg_test/results/default/callset.final.vcf.gz" > "$tmp/samples.txt"

# The toy chr20 contig, plus chrX at full length with only the windows unmasked.
px samtools faidx "$grch38" chrX:2700000-2860000 chrX:155650000-155760000 > "$tmp/windows.fa"
px python - "$jg_test/resources/chr20_toy.fasta" "$tmp/windows.fa" "$grch38.fai" "$tmp/reference.fa" <<'EOF'
import sys

chr20_fasta, windows_fasta, fai, out = sys.argv[1:]


def records(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    yield name, "".join(seq)


def write(f, name, seq):
    f.write(f">{name}\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i : i + 60] + "\n")


chrx_length = next(int(l.split("\t")[1]) for l in open(fai) if l.startswith("chrX\t"))
chrx = bytearray(b"N" * chrx_length)
for name, seq in records(windows_fasta):
    start = int(name.split(":")[1].split("-")[0])
    chrx[start - 1 : start - 1 + len(seq)] = seq.encode()
with open(out, "w") as f:
    for name, seq in records(chr20_fasta):
        write(f, name, seq)
    write(f, "chrX", chrx.decode())
EOF
px bgzip -@ 4 -c "$tmp/reference.fa" > resources/reference.fa.gz
px samtools faidx resources/reference.fa.gz

# 1000 Genomes codes males haploid outside the PARs. GATK joint genotyping at ploidy 2 codes them diploid.
px bcftools view -Ob -o "$tmp/chr20.bcf" "$jg_test/results/default/callset.final.vcf.gz"
px bcftools view -S "$tmp/samples.txt" -r "$chrx_windows" -Ou "$nygc_chrx" \
    | px bcftools annotate -x INFO,^FORMAT/GT,FORMAT/DP -Ou \
    | px bcftools +fixploidy -Ob -o "$tmp/chrX.bcf" -- -f 2
cmp <(px bcftools query -l "$tmp/chr20.bcf") <(px bcftools query -l "$tmp/chrX.bcf")
px bcftools concat -Ob -o "$tmp/callset.bcf" "$tmp/chr20.bcf" "$tmp/chrX.bcf"
px bcftools reheader --fai resources/reference.fa.gz.fai -o "$tmp/callset.reheader.bcf" "$tmp/callset.bcf"
px bcftools view -Oz -o resources/callset.vcf.gz -W=tbi "$tmp/callset.reheader.bcf"

for region in "${!span[@]}"; do
    px bcftools view -r "${span[$region]}" -S "^$tmp/samples.txt" --force-samples \
        -Ob -o "resources/panel/$region.bcf" -W "$res/hgdp1kgp/$region.bcf"
done

px bcftools query -l resources/panel/chr20.bcf > "$tmp/panel_samples.txt"
px python - "$res/hgdp1kgp/gnomad_meta_updated.tsv" "$tmp/panel_samples.txt" resources/panel/gnomad_meta_updated.tsv <<'EOF'
import csv
import sys

metadata, samples, out = sys.argv[1:]
keep = set(open(samples).read().split())
columns = ["s", "hgdp_tgp_meta.Population", "hgdp_tgp_meta.Genetic.region", "sample_filters.release_related"]
with open(metadata) as f, open(out, "w", newline="") as g:
    writer = csv.writer(g, delimiter="\t", lineterminator="\n")
    writer.writerow(columns)
    for row in csv.DictReader(f, delimiter="\t"):
        if row["s"] in keep:
            writer.writerow([row[c] for c in columns])
EOF

gmaps=$res/genetic-maps/shapeit5
zcat "$gmaps/chr20.b38.gmap.gz" | awk 'NR == 1 || $1 <= 500000' | gzip -n > resources/maps/chr20.b38.gmap.gz
zcat "$gmaps/chrX_par1.b38.gmap.gz" | awk 'NR == 1 || $1 >= 2600000' | gzip -n > resources/maps/chrX_par1.b38.gmap.gz
zcat "$gmaps/chrX.b38.gmap.gz" | awk 'NR == 1 || $1 <= 2960000 || $1 >= 155550000' | gzip -n > resources/maps/chrX_nonpar.b38.gmap.gz
zcat "$gmaps/chrX_par2.b38.gmap.gz" | awk 'NR == 1 || $1 <= 155860000' | gzip -n > resources/maps/chrX_par2.b38.gmap.gz
zcat "$res/genetic-maps/eagle2/genetic_map_hg38_withX.txt.gz" \
    | awk 'NR == 1 || ($1 == 20 && $2 <= 500000)' | gzip -n > resources/maps/genetic_map_hg38_withX.txt.gz

awk -F, 'NR <= 8 || /^\[Controls\]/ { print; controls = controls || /^\[Controls\]/; next }
    controls || ($10 == "20" && $11 <= 347695)' \
    "$res/mega/Multi-EthnicGlobal_D2.csv" > resources/Multi-EthnicGlobal_D2.csv
