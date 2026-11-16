# Getting started

Setup, resources, config and running. [PIPELINE.md](PIPELINE.md) describes
each rule and its outputs.

## Environments

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
are certifi, bcftools and htslib, so their security fixes land at once, and
pqplot, which this project's author also maintains.

## Resources

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

## Runs

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

## Config

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
| `kinship_cutoff` | KING kinship above which a pair counts as related. Default `0.0442`, the lower bound of third-degree relatives. |
| `num_pcs` | Number of genetic PCs. Default `20`. |
| `panel_metadata` | gnomAD HGDP+1kGP metadata TSV. It MUST have the columns `s`, `hgdp_tgp_meta.Population` and `sample_filters.release_related`. |
| `local_ancestry.tool` | `rfmix2` (default) or `gnomix`. The Snakefile includes only that tool's rule file. |
| `local_ancestry.gnomix_config` | Optional gnomix config to start from. Default: gnomix's `configs/config_array.yaml`. The rule sets its seed, cores and cleanup. |
| `local_ancestry.populations` | Panel populations of each ancestry, which may enter the reference. The ancestry order sets the codes in the msp files. |
| `local_ancestry.min_purity` | Smallest ADMIXTURE share of its own ancestry for a candidate to join the reference. Default `0.9` for AFR and EUR and `0.8` for NAT. |
| `local_ancestry.pdf_samples` | Samples with a page in `plots/ideograms.pdf`: a count, picked at random with `seed`, or `all`. Default `all`. |
| `local_ancestry.colors` | Optional color of each ancestry in the plots. It MUST name every ancestry of `local_ancestry.populations` and MAY name others. Default: the theme's color cycle. |
| `local_ancestry.rfmix_generations` | RFMix2 `--generations`. Default `8`. Raise or lower it when tracts look too fragmented or too uniform. |
| `seed` | Seed for ADMIXTURE, RFMix2 and the karyogram samples. |
| `filters.min_dp` | Genotypes with lower depth go missing. Default `10`. |
| `filters.min_maf` | Sites with lower minor allele frequency are dropped. Default `0.01`. |
| `filters.max_missing` | Sites with a larger missing genotype fraction are dropped. Default `0.05`. |

## Layout

```
config/
  config.yaml.template   parameters + paths; copy to config/<run_id>/config.yaml
  exclude_samples.tsv.template
                         headerless: sample ID, reason
  <run_id>/              one run's config and exclusions (gitignored)
docs/                    this guide, the pipeline reference, and rule graphs as Graphviz DOT with renders
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

## Conventions

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
