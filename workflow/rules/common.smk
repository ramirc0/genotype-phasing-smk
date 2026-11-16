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
GNOMIX_ENV = pixi(workspace="../envs", env="gnomix", locked=True)

CALLSET = config["callset"]
FASTA = config["fasta"]
PANEL = config["panel"]
FILTERS = config["filters"]

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
CHROMOSOMES = config["chromosomes"]
if not CHROMOSOMES or len(set(CHROMOSOMES)) != len(CHROMOSOMES):
    raise WorkflowError("chromosomes MUST be a non-empty list without repeats.")
if not set(CHROMOSOMES) <= set(AUTOSOMES) | {"chrX"}:
    raise WorkflowError("chromosomes MUST be chr1 to chr22, optionally with chrX.")

# Relatedness, PCs and ancestry use autosomes only.
PHASED_AUTOSOMES = [chrom for chrom in CHROMOSOMES if chrom in AUTOSOMES]

# Each region is phased on its own. bcftools and SHAPEIT5 read the span after the name.
# chrX splits at the GRCh38 PAR bounds, because samples with one X are haploid only between them.
CHRX_REGIONS = {
    "chrX_par1": "chrX:10001-2781479",
    "chrX_nonpar": "chrX:2781480-155701382",
    "chrX_par2": "chrX:155701383-156030895",
}
REGION_SPAN = {chrom: chrom for chrom in PHASED_AUTOSOMES}
if "chrX" in CHROMOSOMES:
    REGION_SPAN |= CHRX_REGIONS

PHASERS = ("shapeit5", "eagle2")
if config["phaser"] not in PHASERS:
    raise WorkflowError(f"phaser MUST be one of {', '.join(PHASERS)}.")
# Eagle2 takes haploid calls, but RFMix2, Tractor and the SHAPEIT5 panel code them as diploid.
if "chrX" in CHROMOSOMES and config["phaser"] != "shapeit5":
    raise WorkflowError("chrX needs phaser: shapeit5.")

LAI = config["local_ancestry"]
LAI_TOOLS = ("rfmix2", "gnomix")
if LAI["tool"] not in LAI_TOOLS:
    raise WorkflowError(f"local_ancestry.tool MUST be one of {', '.join(LAI_TOOLS)}.")
if list(LAI["min_purity"]) != list(LAI["populations"]):
    raise WorkflowError(
        "local_ancestry.min_purity MUST list the ancestries of "
        "local_ancestry.populations, in the same order."
    )
if "colors" in LAI and not set(LAI["populations"]) <= set(LAI["colors"]):
    raise WorkflowError(
        "local_ancestry.colors MUST name every ancestry of local_ancestry.populations."
    )


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
