# genotype-phasing-smk

A [Snakemake](https://github.com/snakemake/snakemake) pipeline that prepares the
genotype inputs of an xQTL mapping SOP from a joint-genotyped WGS
callset.

## Description

The workflow reads one joint-called VCF, such as the output of
[joint-genotyping-smk](https://github.com/ramirc0/joint-genotyping-smk),
and writes every product of the SOP that needs only genotypes. Trait QC and QTL
mapping stay in per-omics pipelines.

### Layout

```
config/
  config.yaml.template   parameters + paths; copy to config/<run_id>/config.yaml
  <run_id>/              one run's config (gitignored)
workflow/
  Snakefile              includes the rule files, onstart hook, rule all
  rules/
    common.smk           config, path constants, helpers
    provenance.smk       resolved config snapshot
  scripts/
    fetch_resources.sh   downloads the panel, genetic maps and array manifest
  envs/pixi.toml         rule envs (pixi workspace)
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
race to install it. Both workspaces exclude packages published in the last 7
days. The exceptions are certifi, bcftools and htslib, so their security fixes
land at once.

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
pixi shell                                                       # launcher env
snakemake -n -p --profile profiles/local --configfile config/<run_id>/config.yaml  # dry run
snakemake --profile profiles/local --configfile config/<run_id>/config.yaml        # local
pixi run slurm <run_id>                                          # SLURM
pixi run slurm <run_id> "all -n"                                 # targets and flags in one quoted arg
pixi run test                                                    # end-to-end test on toy data
```

`pixi run slurm <run_id>` reads `config/<run_id>/config.yaml` and names every
SLURM job `<run_id>_<uuid>`, so `squeue` and `sacct` show which run a job
belongs to.

| Key | Meaning |
|---|---|
| `outdir` | Results root. |
| `pixi_env` | Environment of `workflow/envs/pixi.toml` for most rules. |
| `run_id` | Run name, `YYYY-MM-DD-<what>`. Outputs go under `results/<run_id>/`. |

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

MIT. See [LICENSE](LICENSE).

## Acknowledgments

- [snakemake-template](https://github.com/ramirc0/snakemake-template)
- [Nextstrain Snakemake style guide](https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html)
- [snakemake-software-deployment-plugin-pixi](https://github.com/cademirch/snakemake-software-deployment-plugin-pixi)
