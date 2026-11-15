"""Pick the local ancestry reference from an unsupervised ADMIXTURE run.

ADMIXTURE numbers its components arbitrarily. Ancestries and components are
matched one to one, by the matching with the largest summed mean share of each
ancestry's candidates. A candidate joins the reference when its share of its own
ancestry's component exceeds that ancestry's purity cutoff, as in the SOP. The
other samples, the cohort's, get their ADMIXTURE shares per ancestry.
"""

import argparse
import itertools
import sys

import polars as pl


def build_parser():
    """Build the select_reference.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the select_reference.py flags.
    """
    p = argparse.ArgumentParser(description="Select the local ancestry reference.")
    p.add_argument(
        "--fam", required=True, help="ADMIXTURE input .fam file, as plink2 writes it."
    )
    p.add_argument("--q", required=True, help="ADMIXTURE .Q file.")
    p.add_argument("--candidates", required=True, help="Reference candidates TSV.")
    p.add_argument(
        "--min-purity",
        required=True,
        nargs="+",
        metavar="ANCESTRY=FRACTION",
        help="Purity cutoff of each ancestry, in ancestry order.",
    )
    p.add_argument(
        "--sample-map", required=True, help="Reference sample map TSV to write."
    )
    p.add_argument(
        "--global-ancestry", required=True, help="Cohort ancestry shares TSV to write."
    )
    return p


def match_components(shares, ancestries):
    """Match each ancestry to its ADMIXTURE component.

    Parameters
    ----------
    shares : polars.DataFrame
        One row per candidate: column `ancestry`, then one column per component.
    ancestries : list of str
        Ancestries in config order, as many as components.

    Returns
    -------
    dict of str to str
        Component column of each ancestry.
    """
    components = [c for c in shares.columns if c != "ancestry"]
    means = {
        row.pop("ancestry"): row
        for row in shares.group_by("ancestry").mean().iter_rows(named=True)
    }
    print("Mean component share of the candidates of each ancestry:")
    print(shares.group_by("ancestry").mean().sort("ancestry"))
    best = max(
        itertools.permutations(components),
        key=lambda order: sum(means[a][c] for a, c in zip(ancestries, order)),
    )
    component_of = dict(zip(ancestries, best))
    for ancestry, component in component_of.items():
        top = max(components, key=means[ancestry].get)
        if top != component:
            print(
                f"Warning: {ancestry} candidates peak on {top}, "
                f"but {component} gives the best one-to-one matching."
            )
    print(f"Components: {component_of}")
    return component_of


def main(argv=None):
    """Write the reference sample map and the cohort's ancestry shares.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    min_purity = dict(spec.split("=", 1) for spec in args.min_purity)
    ancestries = list(min_purity)

    samples = pl.read_csv(args.fam, separator="\t", has_header=False).select(
        pl.nth(1).alias("sample")
    )
    q = pl.read_csv(args.q, separator=" ", has_header=False)
    q = q.rename({c: f"component{i}" for i, c in enumerate(q.columns)})
    shares = pl.concat([samples, q], how="horizontal")
    candidates = pl.read_csv(args.candidates, separator="\t").select(
        "sample", "ancestry"
    )

    labelled = shares.join(candidates, on="sample")
    component_of = match_components(labelled.drop("sample"), ancestries)

    reference = (
        labelled.with_columns(
            own_share=pl.coalesce(
                pl.when(pl.col("ancestry") == a).then(pl.col(component_of[a]))
                for a in ancestries
            ),
            cutoff=pl.col("ancestry").replace_strict(
                {a: float(min_purity[a]) for a in ancestries}
            ),
        )
        .filter(pl.col("own_share") > pl.col("cutoff"))
        .with_columns(
            order=pl.col("ancestry").replace_strict(
                ancestries, list(range(len(ancestries)))
            )
        )
        .sort("order", "sample")
    )
    counts = reference.group_by("ancestry").len()
    print("Reference samples per ancestry:")
    print(counts.sort("ancestry"))
    missing = set(ancestries) - set(counts.get_column("ancestry"))
    if missing:
        sys.exit(
            f"No candidate passes the purity cutoff for {', '.join(sorted(missing))}."
        )
    reference.select("sample", "ancestry").write_csv(
        args.sample_map, separator="\t", include_header=False
    )

    shares.join(candidates, on="sample", how="anti").select(
        "sample", *(pl.col(component_of[a]).alias(a) for a in ancestries)
    ).write_csv(args.global_ancestry, separator="\t")


if __name__ == "__main__":
    main()
