"""Summarize local ancestry calls and plot them for review.

The SOP asks to check that ancestry tracts have plausible lengths before mapping.
This script draws both haplotypes of every chromosome for a few random samples,
and compares each sample's genome-wide ancestry, summed from its tracts, with
its ADMIXTURE estimate.
"""

import argparse
import math
import random
import re

import matplotlib.pyplot as plt
import poikilos as pk
import polars as pl


def build_parser():
    """Build the plot_local_ancestry.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_local_ancestry.py flags.
    """
    p = argparse.ArgumentParser(description="Summarize and plot local ancestry.")
    p.add_argument("--msp", required=True, nargs="+", help="msp.tsv per chromosome.")
    p.add_argument("--admixture", required=True, help="Cohort ADMIXTURE shares TSV.")
    p.add_argument("--samples", type=int, default=10, help="Samples to draw.")
    p.add_argument("--seed", type=int, default=1, help="Seed for picking samples.")
    p.add_argument(
        "--global-ancestry", required=True, help="Genome-wide ancestry TSV to write."
    )
    p.add_argument("--karyograms", required=True, help="Karyogram SVG to write.")
    p.add_argument("--comparison", required=True, help="Comparison SVG to write.")
    return p


def read_msp(path):
    """Read an msp.tsv file.

    Parameters
    ----------
    path : str
        An msp.tsv file in RFMix2's layout.

    Returns
    -------
    names : dict of int to str
        Ancestry of each code.
    tracts : polars.DataFrame
        Columns `chrom`, `start` and `end`, then one code column per haplotype,
        named `<sample>.0` and `<sample>.1`.
    """
    with open(path) as f:
        names = {
            int(code): name for name, code in re.findall(r"(\S+)=(\d+)", f.readline())
        }
    tracts = pl.read_csv(path, separator="\t", skip_rows=1)
    tracts = tracts.drop(tracts.columns[3:6]).rename(
        {"#chm": "chrom", "spos": "start", "epos": "end"}
    )
    return names, tracts.with_columns(pl.col("chrom").cast(pl.String))


def tract_lengths(tracts):
    """Sum one chromosome's tract lengths per haplotype and ancestry code.

    Parameters
    ----------
    tracts : polars.DataFrame
        Tracts from `read_msp`.

    Returns
    -------
    polars.DataFrame
        Columns `haplotype`, `code` and `length`.
    """
    return (
        tracts.with_columns(length=pl.col("end") - pl.col("start"))
        .unpivot(
            index="length",
            on=tracts.columns[3:],
            variable_name="haplotype",
            value_name="code",
        )
        .group_by("haplotype", "code")
        .agg(pl.col("length").sum())
    )


def genome_wide(lengths, names):
    """Turn tract lengths into each sample's genome-wide ancestry shares.

    Parameters
    ----------
    lengths : polars.DataFrame
        `tract_lengths` output of every chromosome.
    names : dict of int to str
        Ancestry of each code.

    Returns
    -------
    polars.DataFrame
        One row per sample, one column per ancestry, in code order.
    """
    shares = (
        lengths.with_columns(
            sample=pl.col("haplotype").str.replace(r"\.[01]$", ""),
            ancestry=pl.col("code").replace_strict(names),
        )
        .group_by("sample", "ancestry")
        .agg(pl.col("length").sum())
        .with_columns(pl.col("length") / pl.col("length").sum().over("sample"))
        .pivot(on="ancestry", index="sample", values="length")
    )
    return shares.select(
        "sample", *(pl.col(name).fill_null(0) for name in names.values())
    ).sort("sample")


def plot_karyograms(tracts, samples, colors, path):
    """Draw both haplotypes of every chromosome for each sample.

    Parameters
    ----------
    tracts : polars.DataFrame
        Columns `chrom`, `start`, `end`, `haplotype` and `ancestry`.
    samples : list of str
        Samples to draw, one panel each.
    colors : dict of str to str
        Color of each ancestry.
    path : str
        SVG to write. The PNG goes beside it.
    """
    chroms = sorted(
        tracts.get_column("chrom").unique(), key=lambda c: int(c.removeprefix("chr"))
    )
    nrows = math.ceil(len(samples) / 2)
    fig, axes = plt.subplots(nrows, 2, figsize=(10, 3.5 * nrows), squeeze=False)
    for ax, sample in zip(axes.flat, samples):
        for row, chrom in enumerate(chroms):
            for hap in (0, 1):
                spans = tracts.filter(
                    pl.col("chrom") == chrom, pl.col("haplotype") == f"{sample}.{hap}"
                )
                for ancestry, color in colors.items():
                    ax.broken_barh(
                        [
                            (start / 1e6, (end - start) / 1e6)
                            for start, end in spans.filter(
                                pl.col("ancestry") == ancestry
                            )
                            .select("start", "end")
                            .iter_rows()
                        ],
                        (row + 0.45 * hap, 0.4),
                        facecolors=color,
                    )
        ax.set_yticks([row + 0.425 for row in range(len(chroms))], chroms)
        ax.invert_yaxis()
        ax.set(title=sample, xlabel="Position (Mb)")
    for ax in axes.flat[len(samples) :]:
        ax.set_axis_off()
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colors.values()]
    fig.legend(handles, colors, loc="outside right upper")
    for ax in axes.flat[: len(samples)]:
        pk.despine(ax, categorical_y=True)
    pk.fit_legend(fig)
    pk.save_figure(fig, path)


def plot_comparison(local, admixture, colors, path):
    """Plot each sample's ancestry from its tracts against ADMIXTURE.

    Parameters
    ----------
    local : polars.DataFrame
        Output of `genome_wide`.
    admixture : polars.DataFrame
        One row per sample, one column per ancestry.
    colors : dict of str to str
        Color of each ancestry.
    path : str
        SVG to write. The PNG goes beside it.
    """
    both = local.join(admixture, on="sample", suffix="_admixture")
    fig, axes = plt.subplots(
        1, len(colors), figsize=(3 * len(colors), 3), squeeze=False
    )
    for ax, (ancestry, color) in zip(axes[0], colors.items()):
        ax.scatter(
            both.get_column(f"{ancestry}_admixture"),
            both.get_column(ancestry),
            s=6,
            color=color,
            linewidths=0.25,
        )
        ax.axline((0, 0), slope=1, color="0.6", linewidth=0.75)
        ax.set(title=ancestry, xlabel="ADMIXTURE share", ylabel="Local ancestry share")
    for ax in axes[0]:
        pk.despine(ax)
    pk.save_figure(fig, path)


def main(argv=None):
    """Write the genome-wide table and both figures.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    names, first = read_msp(args.msp[0])
    samples = sorted({c.rsplit(".", 1)[0] for c in first.columns[3:]})
    picked = sorted(
        random.Random(args.seed).sample(samples, min(args.samples, len(samples)))
    )
    picked_haplotypes = [f"{s}.{hap}" for s in picked for hap in (0, 1)]

    lengths, picked_tracts = [], []
    for path in args.msp:
        _, tracts = read_msp(path)
        lengths.append(tract_lengths(tracts))
        picked_tracts.append(
            tracts.select("chrom", "start", "end", *picked_haplotypes).unpivot(
                index=["chrom", "start", "end"],
                variable_name="haplotype",
                value_name="code",
            )
        )
    local = genome_wide(pl.concat(lengths), names)
    local.write_csv(args.global_ancestry, separator="\t")

    pk.use("plain")
    colors = dict(zip(names.values(), pk.PALETTES["plain"].cycle))
    karyogram_tracts = pl.concat(picked_tracts).with_columns(
        ancestry=pl.col("code").replace_strict(names)
    )
    plot_karyograms(karyogram_tracts, picked, colors, args.karyograms)
    admixture = pl.read_csv(args.admixture, separator="\t")
    plot_comparison(local, admixture, colors, args.comparison)


if __name__ == "__main__":
    main()
