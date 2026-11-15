"""List the panel samples that may enter the local ancestry reference.

A candidate belongs to one of the configured populations, is in the panel BCF,
and is not flagged as related in the gnomAD HGDP+1kGP metadata. ADMIXTURE then
decides which candidates are pure enough to keep.
"""

import argparse

import polars as pl

POPULATION = "hgdp_tgp_meta.Population"
RELATED = "sample_filters.release_related"


def build_parser():
    """Build the reference_candidates.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the reference_candidates.py flags.
    """
    p = argparse.ArgumentParser(description="List reference panel candidates.")
    p.add_argument("--metadata", required=True, help="gnomAD HGDP+1kGP metadata TSV.")
    p.add_argument(
        "--panel-samples", required=True, help="Panel sample IDs, one per line."
    )
    p.add_argument(
        "--populations",
        required=True,
        nargs="+",
        metavar="ANCESTRY=POP,POP",
        help="Panel populations of each ancestry, in ancestry order.",
    )
    p.add_argument(
        "--output", required=True, help="TSV to write: sample, population, ancestry."
    )
    return p


def parse_populations(specs):
    """Map each panel population to its ancestry.

    Parameters
    ----------
    specs : list of str
        One `ANCESTRY=POP,POP` string per ancestry.

    Returns
    -------
    polars.DataFrame
        Columns `population` and `ancestry`.
    """
    rows = []
    for spec in specs:
        ancestry, populations = spec.split("=", 1)
        rows += [(population, ancestry) for population in populations.split(",")]
    return pl.DataFrame(rows, schema=["population", "ancestry"], orient="row")


def main(argv=None):
    """Write the candidates and print their counts per population.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    ancestry_of = parse_populations(args.populations)
    panel = pl.read_csv(
        args.panel_samples, has_header=False, new_columns=["sample"]
    ).get_column("sample")

    candidates = (
        pl.read_csv(args.metadata, separator="\t", infer_schema_length=0)
        .filter(pl.col(RELATED) != "true", pl.col("s").is_in(panel.implode()))
        .select(pl.col("s").alias("sample"), pl.col(POPULATION).alias("population"))
        .join(ancestry_of, on="population")
        .sort("ancestry", "population", "sample")
    )
    candidates.write_csv(args.output, separator="\t")
    with pl.Config(tbl_rows=-1):
        print(candidates.group_by("ancestry", "population").len().sort("ancestry"))


if __name__ == "__main__":
    main()
