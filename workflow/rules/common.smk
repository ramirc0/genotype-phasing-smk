# Config, path constants, and shared helpers.

import re

from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RUN_ID = config.get("run_id", "default")
OUTDIR = f"{config['outdir']}/{RUN_ID}"
REFDIR = f"{config['outdir']}/refs"
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)
SHAPEIT5_ENV = pixi(workspace="../envs", env="shapeit5", locked=True)
EAGLE2_ENV = pixi(workspace="../envs", env="eagle2", locked=True)

CALLSET = config["callset"]
FASTA = config["fasta"]
PANEL = config["panel"]
FILTERS = config["filters"]

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
CHROMOSOMES = config["chromosomes"]
if not CHROMOSOMES or len(set(CHROMOSOMES)) != len(CHROMOSOMES):
    raise WorkflowError("chromosomes MUST be a non-empty list without repeats.")
if not set(CHROMOSOMES) <= set(AUTOSOMES):
    raise WorkflowError("chromosomes MUST be autosomes, chr1 to chr22.")

# Each region is phased on its own. bcftools and SHAPEIT5 read the span after the name.
REGION_SPAN = {chrom: chrom for chrom in CHROMOSOMES}

PHASERS = ("shapeit5", "eagle2")
if config["phaser"] not in PHASERS:
    raise WorkflowError(f"phaser MUST be one of {', '.join(PHASERS)}.")


wildcard_constraints:
    region="|".join(map(re.escape, REGION_SPAN)),
    chrom="|".join(map(re.escape, CHROMOSOMES)),


def fasta_indexes(fasta):
    """Return the index files htslib needs beside a FASTA.

    Parameters
    ----------
    fasta : str
        FASTA path, plain or bgzipped.

    Returns
    -------
    list of str
        The `.fai`, plus the `.gzi` of a bgzipped FASTA.
    """
    indexes = [f"{fasta}.fai"]
    if fasta.endswith(".gz"):
        indexes.append(f"{fasta}.gzi")
    return indexes
