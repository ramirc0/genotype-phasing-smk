# genotype-phasing-smk

A [Snakemake](https://github.com/snakemake/snakemake) pipeline that prepares the
genotype inputs of an xQTL mapping SOP from a joint-genotyped WGS
callset.

## Description

The workflow reads one joint-called VCF, such as the output of
[joint-genotyping-smk](https://github.com/ramirc0/joint-genotyping-smk),
and writes every product of the SOP that needs only genotypes. Trait QC and QTL
mapping stay in per-omics pipelines.

### Pipeline

Each chromosome runs on its own:

1. `prepare_genotypes` keeps PASS records and drops the samples in
   `exclude_samples`. It splits multiallelics with the other ALT set missing, so
   no sample gains a REF allele it lacks, and drops `*` alleles. Genotypes with
   depth below `min_dp` go missing. Sites below `min_maf` or above
   `max_missing` are dropped. Only GT is kept.
2. The phaser set by `phaser` phases against the reference panel. Both
   phasers impute missing genotypes and drop sites absent from the panel, and
   their logs count them.
   - `phase_shapeit5` (default) runs SHAPEIT5 `phase_common`. Panel and cohort
     haplotypes condition each other. The SHAPEIT5 docs skip `phase_rare` below
     2,000 samples.
   - `phase_eagle2` runs Eagle2, the SOP's phaser. Eagle's auto rule runs one
     iteration when the cohort is under half the panel, which phases each
     sample against the panel alone. The rule sets `--pbwtIters 2`, so cohort
     haplotypes inform each other.
3. `finalize_phased` recomputes AC, AN, AF and MAF and sets each ID to
   `CHROM:POS:REF:ALT`.
4. `export_pgen` writes a plink2 pgen copy, which keeps the phase. Its
   chromosome codes drop the `chr` prefix, as tensorQTL's SOP example expects.

With `chrX` in `chromosomes`, chrX runs as three regions split at the GRCh38 PAR
bounds, because samples with one X are haploid only between them:

- `infer_x_ploidy` imputes each sample's sex from its X F-statistic with plink2
  on the prepared non-PAR genotypes. The F-statistic counts X copies, which is
  what phasing needs: XY, XYY and X0 samples come out haploid, and XX, XXY and
  XXX samples diploid. Samples plink2 cannot call stay diploid.
- `phase_shapeit5` passes the haploid samples to `--haploids` outside the PARs.
  SHAPEIT5 sets their hets missing and imputes them as homozygous.
- `concat_chrx` joins the three regions before `finalize_phased`.
- `export_pgen` splits the PARs with `--split-par hg38` and takes the imputed
  sex, so PAR calls stay diploid in every sample.

chrX needs `phaser: shapeit5`. Eagle2 takes haploid calls, but RFMix2, Tractor
and the SHAPEIT5 panel code them as diploid.

Across the genome:

5. `make_array_sites` reads the autosomal SNP positions from the array
   manifest once, into `results/refs/`. `array_genotypes` keeps the phased
   SNPs at those positions, as the SOP thins WGS genotypes to an array SNP set
   before local ancestry. A split multiallelic site leaves several SNP records
   at one position, and RFMix2 matches records by position alone, so only
   positions with one SNP record are kept. `cohort_array_pgen` joins the
   chromosomes into one pgen.
6. `relatedness` computes KING-robust kinship with plink2, which tolerates
   population structure. `--king-cutoff` drops a sample from each pair above
   `kinship_cutoff` and writes the rest to `unrelated.txt`.
7. `genetic_pcs` computes PCs from the unrelated samples, on SNPs with MAF of
   at least 5% outside long-range LD regions, pruned to r² < 0.1 in 1 Mb
   windows. It then projects every sample, so relatives get PCs on the same
   scale.

8. `reference_candidates` lists the panel samples of the populations in
   `local_ancestry.populations` that gnomAD does not flag as related.
   `match_array_sites` keeps the array SNPs whose REF and ALT agree exactly
   between the cohort and the candidates, because RFMix2 matches records by
   position alone.
9. `admixture` runs unsupervised ADMIXTURE on the unrelated cohort samples and
   the candidates, with one component per ancestry, on pruned common SNPs.
   `select_reference` matches ancestries to components one to one, by the
   matching with the largest summed mean share of each ancestry's candidates.
   A candidate joins the reference when its share of its own ancestry's
   component exceeds `local_ancestry.min_purity`, as in the SOP. The cohort's
   shares go to `global_ancestry.tsv`.
10. The tool set by `local_ancestry.tool` calls local ancestry per chromosome.
    Both tools read the same files and write msp files in RFMix2's layout, which
    Tractor reads.
    - `local_ancestry_rfmix2` (default) runs RFMix2, the SOP's tool.
    - `local_ancestry_gnomix` runs gnomix. It trains a model on admixed samples
      that it simulates from the reference, then deletes them. Phase correction
      stays off, so its haplotypes stay those of the phased VCF.
11. `extract_tracts` runs Tractor's `extract_tracts.py` on each phased VCF and
    its msp file. It writes, per ancestry, the ALT dosage and the haplotype
    count of each sample at every phased site, which Tractor's regression reads.
    The local ancestry query comes from the phased VCF, so samples and
    haplotypes line up. Sites outside every msp window are skipped.
12. `plot_local_ancestry` draws both haplotypes of every chromosome for 10
    random samples, so tract lengths can be checked before mapping, as the SOP
    asks. It sums each sample's tracts into genome-wide shares and plots them
    against ADMIXTURE.

The callset MUST code every genotype as diploid, chrX included, as GATK joint
genotyping at ploidy 2 does. Relatives stay in phasing. `unrelated.txt` filters samples at
mapping.

### Layout

```
config/
  config.yaml.template   parameters + paths; copy to config/<run_id>/config.yaml
  exclude_samples.tsv.template
                         headerless: sample ID, reason
  <run_id>/              one run's config and exclusions (gitignored)
workflow/
  Snakefile              includes the rule files, onstart hook, rule all
  rules/
    common.smk           config, path constants, helpers
    prepare.smk          per-chromosome genotype filters
    shapeit5.smk         phasing with SHAPEIT5
    eagle2.smk           phasing with Eagle2
    phase.smk            final phased VCFs, pgen copies
    relatedness.smk      array-site genotypes, kinship, unrelated set, genetic PCs
    ancestry.smk         reference candidates, ADMIXTURE, reference selection, plots
    rfmix2.smk           local ancestry with RFMix2
    gnomix.smk           local ancestry with gnomix
    provenance.smk       resolved config snapshot
  scripts/
    fetch_resources.sh   downloads the panel, genetic maps, array manifest and LD regions
    array_sites.py       autosomal SNP positions of an Illumina manifest
    reference_candidates.py
                         panel samples that may enter the local ancestry reference
    select_reference.py  ADMIXTURE component matching and purity filter
    plot_local_ancestry.py
                         karyograms and genome-wide ancestry from local calls
    gnomix_config.py     sets a gnomix config's seed, cores and cleanup per job
    vendor/extract_tracts.py
                         Tractor's script, unchanged; see THIRD_PARTY_LICENSES
  envs/pixi.toml         rule envs (pixi workspace)
  envs/gnomix/           gnomix package recipe at a pinned commit, with a bug fix
.tests/integration/      toy inputs for the end-to-end test; make_test_data.sh records how they were built
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env + test, slurm and lint tasks
pyproject.toml           ruff + numpydoc config
```

### Outputs

Everything goes under `results/<run_id>/`. Logs and benchmarks go to
`logs/<run_id>/` and `benchmarks/<run_id>/`.

| Path | Contents |
|---|---|
| `prepare/<chrom>.bcf` | Filtered genotypes, the phaser's input |
| `phase/<chrom>.bcf` | Phaser output |
| `phased/<chrom>.vcf.gz` | Phased genotypes with recomputed tags and `CHROM:POS:REF:ALT` IDs |
| `pgen/<chrom>.{pgen,pvar,psam}` | plink2 copy of the phased genotypes |
| `x_ploidy.tsv` | plink2 sex imputation: X F-statistic and imputed sex per sample (chrX runs only) |
| `haploid_samples.txt` | Samples with one X, phased as haploid outside the PARs (chrX runs only) |
| `array/<chrom>.bcf` | Phased genotypes at array SNPs |
| `array/cohort.{pgen,pvar,psam}` | Array-site genotypes of all chromosomes |
| `king.kin0` | Sample pairs with KING kinship above `kinship_cutoff` |
| `unrelated.txt` | Sample IDs of the unrelated set, one per line |
| `pcs.tsv` | `num_pcs` genetic PCs of every sample |
| `pcs.eigenval` | Eigenvalues of those PCs |
| `ancestry/candidates.tsv` | Reference candidates: sample, population, ancestry |
| `ancestry/<chrom>.{query,reference}.bcf` | Cohort and candidates at the matched array SNPs |
| `admixture/merged.*` | ADMIXTURE input and its `.Q` and `.P` files |
| `lai_sample_map.tsv` | Local ancestry reference: sample, ancestry |
| `global_ancestry.tsv` | ADMIXTURE share per ancestry of each unrelated cohort sample |
| `local_ancestry/<chrom>.{msp,fb}.tsv` | Local ancestry calls in RFMix2's layout |
| `local_ancestry_global.tsv` | Genome-wide share per ancestry of each sample, from its tracts |
| `tracts/<chrom>.anc<i>.{dosage,hapcount}.txt.gz` | Tractor input per ancestry `i`, numbered as in the msp header |
| `plots/karyograms.{svg,png}` | Local ancestry of 10 random samples |
| `plots/global_ancestry.{svg,png}` | Genome-wide shares from local calls against ADMIXTURE |
| `config.snapshot.json` | The resolved config, `--config` overrides included. Rerun it with `--configfile results/<run_id>/config.snapshot.json`. |

### Conventions

The workflow follows the
[Nextstrain Snakemake style guide](https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html)
and [Snakemake's best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html):

- No `run:` blocks. Custom logic lives in `workflow/scripts/`, called from `shell:`.
- Raw, triple-quoted `shell` blocks (`r"""`), one option per line.
- Every interpolation uses `:q` quoting.
- Config reaches `shell` through `params:` (a `lambda`), never directly.
- `input:` and `output:` are literal path strings. All paths are relative.
- Rules tee logs to file and terminal (`exec &> >(tee {log:q})`).
- Shell blocks end with `sync {output:q}`. On a shared filesystem, Snakemake's
  first stat of an output from the submit host otherwise waits for the
  writing node to flush it.
- Every rule has `log:` and `benchmark:`. Each names its pixi env in `software:`.
- Every rule declares its real outputs. There are no sentinel files and no `localrules`.

## Getting Started

### Dependencies

- Linux x86-64 and [pixi](https://pixi.sh) on `PATH`.
- Network access to GitHub. Snakemake and its plugins install from git.
- A SLURM cluster for `profiles/slurm`.

### Installing

```bash
git clone <repo> genotype-phasing-smk
cd genotype-phasing-smk
pixi install
```

The launcher env (`pixi.toml` at the root) has Snakemake main and its plugins
for pixi deployment and SLURM. The rule envs live in a second pixi workspace,
`workflow/envs/`. Its `pixi.lock` pins every package. An `onstart` hook
installs it before the first job, because parallel SLURM jobs would otherwise
race to install it.

- `default` runs bcftools, plink2, ADMIXTURE, RFMix2 and the Python scripts.
- `shapeit5` and `eagle2` each hold one phaser, because both pin htslib below
  1.25.
- `gnomix` runs gnomix, built by `workflow/envs/gnomix/` from upstream commit
  `cd15f65`, because gnomix has no releases or package. The recipe patches one
  call that crashes when gnomix deletes its simulated data.

Both workspaces exclude packages published in the last 7 days. The exceptions
are certifi, bcftools and htslib, so their security fixes land at once.

Fetch the reference resources once, about 20 GB:

```bash
pixi run fetch-resources /data/zusers/ramirezc/data/genotype-phasing-smk/resources
```

| Directory | Contents |
|---|---|
| `hgdp1kgp/` | The [HGDP+1kGP phased panel](https://doi.org/10.1101/gr.278378.123) v2 from gnomAD: 4,091 samples phased with SHAPEIT5, one BCF per autosome and per chrX region, and `gnomad_meta_updated.tsv` |
| `genetic-maps/shapeit5/` | SHAPEIT5 GRCh38 maps from `odelaneau/shapeit` at commit `c34d4db` |
| `genetic-maps/eagle2/` | Eagle's `genetic_map_hg38_withX.txt.gz` |
| `mega/` | The Illumina Multi-Ethnic Global-8 D2 manifest, GRCh38 |
| `high-LD-regions-hg38-GRCh38.bed` | Long-range LD regions from [plinkQC](https://github.com/meyer-lab-cshl/plinkQC) at commit `2e64f66` |

A rerun skips finished files and resumes partial ones. Each panel BCF and map
is also linked under its region name, `chr1` to `chr22`, `chrX_par1`,
`chrX_nonpar` or `chrX_par2`, which the config paths use.

### Executing program

Each run keeps its config in `config/<run_id>/`. A run ID is the start date,
then what the run covers, such as `2026-11-15-cohort-wgs`. There is no
`config/config.yaml`, so a run without `--configfile` fails instead of falling
back to a stale default.

```bash
mkdir -p config/<run_id>
cp config/config.yaml.template config/<run_id>/config.yaml       # then edit
cp config/exclude_samples.tsv.template config/<run_id>/exclude_samples.tsv  # then edit
pixi shell                                                       # launcher env
snakemake -n -p --profile profiles/local --configfile config/<run_id>/config.yaml  # dry run
snakemake --profile profiles/local --configfile config/<run_id>/config.yaml        # local
pixi run slurm <run_id>                                          # SLURM
pixi run slurm <run_id> "all -n"                                 # targets and flags in one quoted arg
pixi run test                                                    # end-to-end test on toy data
pixi run test-options                                            # the same with the non-default options
```

`pixi run slurm <run_id>` reads `config/<run_id>/config.yaml` and names every
SLURM job `<run_id>_<uuid>`, so `squeue` and `sacct` show which run a job
belongs to.

| Key | Meaning |
|---|---|
| `outdir` | Results root. |
| `pixi_env` | Environment of `workflow/envs/pixi.toml` for most rules. |
| `run_id` | Run name, `YYYY-MM-DD-<what>`. Outputs go under `results/<run_id>/`. |
| `callset` | Joint-called VCF with a `.tbi`. |
| `fasta` | Reference FASTA of the callset, with its `.fai` (and `.gzi` if bgzipped). |
| `exclude_samples` | Headerless TSV of sample ID and reason. Lines starting with `#` are ignored. Every ID MUST be in the callset. |
| `chromosomes` | Chromosomes to phase, `chr1` to `chr22`, optionally with `chrX`. Relatedness, PCs and ancestry use the autosomes. |
| `phaser` | `shapeit5` (default) or `eagle2`. The Snakefile includes only that phaser's rule file. chrX needs `shapeit5`. |
| `panel` | Phased reference panel BCF per region, with `{region}` in the path and a `.csi` beside it. |
| `genetic_maps.shapeit5` | SHAPEIT5 genetic map per region, with `{region}` in the path. |
| `genetic_maps.eagle2` | Eagle2 genetic map, one file for all chromosomes. |
| `array_manifest` | Illumina manifest CSV whose autosomal SNPs form the array SNP set. |
| `long_range_ld` | BED of long-range LD regions left out of the PCs. |
| `kinship_cutoff` | KING kinship above which a pair counts as related. Default `0.0442`, the third-degree bound, as in the SOP. |
| `num_pcs` | Number of genetic PCs. Default `20`. The SOP uses 10 for local ancestry and 5 for bulk eQTL. |
| `panel_metadata` | gnomAD HGDP+1kGP metadata TSV. It MUST have the columns `s`, `hgdp_tgp_meta.Population` and `sample_filters.release_related`. |
| `local_ancestry.tool` | `rfmix2` (default) or `gnomix`. The Snakefile includes only that tool's rule file. |
| `local_ancestry.gnomix_config` | Optional gnomix config to start from. Default: gnomix's `configs/config_array.yaml`. The rule sets its seed, cores and cleanup. |
| `local_ancestry.populations` | Panel populations of each ancestry, which may enter the reference. The ancestry order sets the codes in the msp files. |
| `local_ancestry.min_purity` | Smallest ADMIXTURE share of its own ancestry for a candidate to join the reference. Default `0.9` for AFR and EUR and `0.8` for NAT, as in the SOP. |
| `local_ancestry.rfmix_generations` | RFMix2 `--generations`. Default `8`. The SOP raises or lowers it when tracts look too fragmented or too uniform. |
| `seed` | Seed for ADMIXTURE, RFMix2 and the karyogram samples. |
| `filters.min_dp` | Genotypes with lower depth go missing. Default `10`, as in the SOP. |
| `filters.min_maf` | Sites with lower minor allele frequency are dropped. Default `0.01`, as in the SOP. |
| `filters.max_missing` | Sites with a larger missing genotype fraction are dropped. Default `0.05`. |

## Help

- Every rule reran. Any change to the rule env manifest or lock reruns every
  rule. So does moving the workdir, because the env hash includes the
  workspace's absolute path.
- A job failed. The profiles set `--show-failed-logs`. The full log sits under
  `logs/<run_id>/`.
- A SLURM node keeps failing jobs. List it under `slurm-exclude-failed-nodes`
  in `profiles/slurm/config.yaml`.

## Authors

Christian Ramirez ([@ramirc0](https://github.com/ramirc0))

## License

MIT. See [LICENSE](LICENSE). `workflow/scripts/vendor/extract_tracts.py` comes
from [Tractor](https://github.com/Atkinson-Lab/Tractor) under its MIT license.
See [THIRD_PARTY_LICENSES](THIRD_PARTY_LICENSES).

## Acknowledgments

- [snakemake-template](https://github.com/ramirc0/snakemake-template)
- [Nextstrain Snakemake style guide](https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html)
- [snakemake-software-deployment-plugin-pixi](https://github.com/cademirch/snakemake-software-deployment-plugin-pixi)
