"""Summarize local ancestry calls and plot them for review.

Ancestry tracts need plausible lengths before mapping. This script draws both
haplotypes of every chromosome for a few random samples, as plain bars and on
banded chromosomes, and compares each sample's genome-wide ancestry, summed
from its tracts, with its ADMIXTURE estimate.
"""

import argparse
import io
import math
import random
import re

import matplotlib.pyplot as plt
import poikilos as pk
import polars as pl
import pqplot as pqp
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfWriter


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
    p.add_argument(
        "--colors",
        nargs="*",
        default=[],
        metavar="ANCESTRY=COLOR",
        help="Color of each ancestry. Default: the theme's color cycle.",
    )
    p.add_argument(
        "--pdf-samples",
        default="all",
        help="Samples to give a PDF page: a count, or 'all'.",
    )
    p.add_argument(
        "--ideograms", required=True, help="PDF to write, one karyogram per page."
    )
    p.add_argument("--comparison", required=True, help="Comparison SVG to write.")
    return p


def pick(samples, count, seed):
    """Pick a sorted random subset of samples.

    Parameters
    ----------
    samples : list of str
        Samples to pick from.
    count : int
        Samples to pick, or all of them when there are fewer.
    seed : int
        Seed of the pick.

    Returns
    -------
    list of str
        The picked samples, sorted.
    """
    return sorted(random.Random(seed).sample(samples, min(count, len(samples))))


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
    return names, tracts.with_columns(
        pl.col("chrom").cast(pl.String), pl.col(tracts.columns[3:]).cast(pl.Int8)
    )


def sample_tracts(tracts, sample, names):
    """Gather one sample's tracts from every chromosome.

    Long tables are built one sample at a time, because all haplotypes at once
    would hold every segment of every haplotype in memory. An msp segment ends
    wherever any haplotype changes ancestry, so a haplotype's consecutive
    segments of one ancestry merge into one tract.

    Parameters
    ----------
    tracts : list of polars.DataFrame
        `read_msp` tracts, one per chromosome.
    sample : str
        Sample whose tracts to gather.
    names : dict of int to str
        Ancestry of each code.

    Returns
    -------
    polars.DataFrame
        Columns `chrom`, `start`, `end`, `copy` (0 or 1) and `ancestry`.
    """
    segments = pl.concat(
        t.select(
            "chrom",
            "start",
            "end",
            pl.col(f"{sample}.0").alias("0"),
            pl.col(f"{sample}.1").alias("1"),
        ).unpivot(
            index=["chrom", "start", "end"], variable_name="copy", value_name="code"
        )
        for t in tracts
    ).sort("chrom", "copy", "start")
    changed = pl.col("code") != pl.col("code").shift().over("chrom", "copy")
    return (
        segments.with_columns(tract=changed.fill_null(True).cum_sum())
        .group_by("tract", maintain_order=True)
        .agg(
            pl.col("chrom").first(),
            pl.col("start").first(),
            pl.col("end").last(),
            pl.col("copy").first().cast(pl.Int8),
            pl.col("code").first().replace_strict(names).alias("ancestry"),
        )
        .drop("tract")
    )


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


def plot_karyograms(tracts, samples, names, colors, path):
    """Draw both haplotypes of every chromosome for each sample.

    Parameters
    ----------
    tracts : list of polars.DataFrame
        `read_msp` tracts, one per chromosome.
    samples : list of str
        Samples to draw, one panel each.
    names : dict of int to str
        Ancestry of each code.
    colors : dict of str to str
        Color of each ancestry.
    path : str
        SVG to write. The PNG goes beside it.
    """
    chroms = sorted(
        (t.get_column("chrom")[0] for t in tracts),
        key=lambda c: int(c.removeprefix("chr")),
    )
    nrows = math.ceil(len(samples) / 2)
    fig, axes = plt.subplots(nrows, 2, figsize=(10, 3.5 * nrows), squeeze=False)
    for ax, sample in zip(axes.flat, samples):
        own = sample_tracts(tracts, sample, names)
        for row, chrom in enumerate(chroms):
            for hap in (0, 1):
                spans = own.filter(pl.col("chrom") == chrom, pl.col("copy") == hap)
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


def write_ideogram_pages(tracts, samples, names, local, colors, path):
    """Write a PDF with one page per sample: its tracts on banded hg38 chromosomes.

    Each page carries a bookmark and a page label named after its sample, so a
    viewer's outline or page box jumps to a sample.

    Parameters
    ----------
    tracts : list of polars.DataFrame
        `read_msp` tracts, one per chromosome.
    samples : list of str
        Samples to give a page, in page order.
    names : dict of int to str
        Ancestry of each code.
    local : polars.DataFrame
        Output of `genome_wide`, for the shares in each page title.
    colors : dict of str to str
        Color of each ancestry.
    path : str
        PDF to write.
    """
    hg38 = pqp.Genome.bundled("hg38")
    present = {t.get_column("chrom")[0] for t in tracts}
    chroms = [c for c in hg38.autosomes if c in present]
    shares = {row["sample"]: row for row in local.iter_rows(named=True)}
    pages = io.BytesIO()
    with PdfPages(pages) as pdf:
        for sample in samples:
            kary = pqp.Karyogram(hg38, chroms=chroms, groups=2, copies=2)
            kary.paint(
                sample_tracts(tracts, sample, names),
                start="start",
                end="end",
                color="ancestry",
                palette=colors,
                copy="copy",
            )
            kary.legend(palette=colors)
            summary = ", ".join(f"{a} {shares[sample][a]:.2f}" for a in colors)
            kary.figure.suptitle(f"{sample}: {summary}")
            pdf.savefig(kary.figure)
            plt.close(kary.figure)

    writer = PdfWriter(clone_from=pages)
    for page, sample in enumerate(samples):
        writer.add_outline_item(sample, page)
        writer.set_page_label(page, page, prefix=sample)
    writer.page_mode = "/UseOutlines"
    writer.write(path)


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
    picked = pick(samples, args.samples, args.seed)
    paged = (
        samples
        if args.pdf_samples == "all"
        else pick(samples, int(args.pdf_samples), args.seed)
    )
    tracts = [read_msp(path)[1] for path in args.msp]
    local = genome_wide(pl.concat(tract_lengths(t) for t in tracts), names)
    local.write_csv(args.global_ancestry, separator="\t")

    pk.use("plain")
    chosen = dict(pair.split("=", 1) for pair in args.colors)
    cycle = pk.PALETTES["plain"].cycle
    colors = {name: chosen.get(name, c) for name, c in zip(names.values(), cycle)}
    plot_karyograms(tracts, picked, names, colors, args.karyograms)
    write_ideogram_pages(tracts, paged, names, local, colors, args.ideograms)
    admixture = pl.read_csv(args.admixture, separator="\t")
    plot_comparison(local, admixture, colors, args.comparison)


if __name__ == "__main__":
    main()
