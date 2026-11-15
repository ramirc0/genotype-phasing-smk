#!/usr/bin/env bash
# Fetch the reference panel, genetic maps and array manifest into one directory.
# Each panel file and map is also linked under its region name, which the config uses.
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 <resources_dir>" >&2
    exit 2
fi
dest=$1

gnomad=https://storage.googleapis.com/gcp-public-data--gnomad
panel_url=$gnomad/resources/hgdp_1kg/phased_haplotypes_v2
metadata_url=$gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv
shapeit_maps_url=https://raw.githubusercontent.com/odelaneau/shapeit/c34d4db3e99a2f7e23deb727671ae260901a5886/resources/maps/b38
eagle_map_url=https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
mega_url=https://webdata.illumina.com/downloads/productfiles/multiethnic-global-8/v1-0/build38/multi-ethnic-global-8-d2-csv.zip

# A download keeps its final name only once complete, so a rerun skips it and resumes the rest.
fetch() {
    local url=$1 out=$2
    [[ -e $out ]] && return
    curl -fL --retry 5 -C - -o "$out.part" "$url"
    mv "$out.part" "$out"
}

mkdir -p "$dest"/{hgdp1kgp,genetic-maps/shapeit5,genetic-maps/eagle2,mega}

declare -A panel_file
for c in {1..22}; do
    panel_file[chr$c]=hgdp1kgp_chr$c.filtered.SNV_INDEL.phased.shapeit5.bcf
done
panel_file[chrX_par1]=hgdp1kgp_chrX_par1.shapeit5_common.bcf
panel_file[chrX_nonpar]=hgdp1kgp_chrX_non_par.full.shapeit5_rare.bcf
panel_file[chrX_par2]=hgdp1kgp_chrX_par2.shapeit5_common.bcf

for region in "${!panel_file[@]}"; do
    f=${panel_file[$region]}
    fetch "$panel_url/$f" "$dest/hgdp1kgp/$f"
    fetch "$panel_url/$f.csi" "$dest/hgdp1kgp/$f.csi"
    ln -sfn "$f" "$dest/hgdp1kgp/$region.bcf"
    ln -sfn "$f.csi" "$dest/hgdp1kgp/$region.bcf.csi"
done
fetch "$metadata_url" "$dest/hgdp1kgp/gnomad_meta_updated.tsv"

for map in chr{1..22} chrX chrX_par1 chrX_par2; do
    fetch "$shapeit_maps_url/$map.b38.gmap.gz" "$dest/genetic-maps/shapeit5/$map.b38.gmap.gz"
done
ln -sfn chrX.b38.gmap.gz "$dest/genetic-maps/shapeit5/chrX_nonpar.b38.gmap.gz"

fetch "$eagle_map_url" "$dest/genetic-maps/eagle2/genetic_map_hg38_withX.txt.gz"

fetch "$mega_url" "$dest/mega/multi-ethnic-global-8-d2-csv.zip"
unzip -o -d "$dest/mega" "$dest/mega/multi-ethnic-global-8-d2-csv.zip"
