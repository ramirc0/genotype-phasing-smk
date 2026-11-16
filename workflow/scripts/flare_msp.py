"""Turn FLARE's per-marker local ancestry calls into an msp file in RFMix2's layout.

FLARE reports the ancestry of each haplotype at each marker. Consecutive markers
where no haplotype changes ancestry merge into one segment. As in RFMix2's msp
files, each segment runs from its first marker to the next segment's first
marker, and the last one ends at the last marker. Tractor reads a site into the
segment with ``spos <= pos < epos``, so every site between markers lands in one.
Ancestry codes follow the order given by ``--ancestries``.
"""

import argparse
import re

import numpy as np
import polars as pl


def build_parser():
    """Build the flare_msp.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the flare_msp.py flags.
    """
    p = argparse.ArgumentParser(description="Write FLARE calls as an msp file.")
    p.add_argument(
        "--calls",
        required=True,
        help="Headerless TSV: CHROM, POS, then AN1 and AN2 of each sample.",
    )
    p.add_argument(
        "--samples", required=True, help="Sample IDs in the calls' column order."
    )
    p.add_argument(
        "--header",
        required=True,
        help="File holding FLARE's '##ANCESTRY=<...>' line.",
    )
    p.add_argument(
        "--ancestries", required=True, nargs="+", help="Ancestries in msp code order."
    )
    p.add_argument(
        "--map", required=True, help="Genetic map TSV: chromosome, position, cM."
    )
    p.add_argument("--output", required=True, help="msp TSV to write.")
    return p


def code_lookup(header, ancestries):
    """Map FLARE's ancestry indices to the indices of `ancestries`.

    Parameters
    ----------
    header : str
        FLARE's ``##ANCESTRY=<AFR=0,EUR=1,...>`` line.
    ancestries : list of str
        Ancestries in msp code order. It MUST hold the same names as `header`.

    Returns
    -------
    numpy.ndarray
        Element ``i`` is the msp code of FLARE's ancestry ``i``.
    """
    flare = {name: int(code) for name, code in re.findall(r"(\w+)=(\d+)", header)}
    if sorted(flare) != sorted(ancestries):
        raise ValueError(f"FLARE ancestries {sorted(flare)} differ from {ancestries}")
    lookup = np.empty(len(flare), dtype=np.int8)
    for code, name in enumerate(ancestries):
        lookup[flare[name]] = code
    return lookup


def segments(codes):
    """Find runs of markers where no haplotype changes ancestry.

    Parameters
    ----------
    codes : numpy.ndarray
        Ancestry codes, one row per marker and one column per haplotype.

    Returns
    -------
    starts : numpy.ndarray
        First marker of each run.
    ends : numpy.ndarray
        The marker after each run's last.
    """
    changes = np.flatnonzero(np.any(codes[1:] != codes[:-1], axis=1)) + 1
    starts = np.concatenate(([0], changes))
    ends = np.concatenate((changes, [len(codes)]))
    return starts, ends


def main(argv=None):
    """Read FLARE's calls and write them as an msp file.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.samples) as f:
        samples = f.read().split()
    with open(args.header) as f:
        lookup = code_lookup(f.read(), args.ancestries)

    calls = pl.read_csv(args.calls, separator="\t", has_header=False)
    chrom = calls[0, 0]
    pos = calls.get_column(calls.columns[1]).to_numpy()
    codes = lookup[calls.select(calls.columns[2:]).to_numpy().astype(np.int8)]

    genetic = pl.read_csv(
        args.map, separator="\t", has_header=False, new_columns=["chrom", "pos", "cm"]
    )
    cm = np.interp(pos, genetic["pos"].to_numpy(), genetic["cm"].to_numpy())

    starts, ends = segments(codes)
    epos = np.concatenate((pos[starts[1:]], [pos[-1]]))
    egpos = np.concatenate((cm[starts[1:]], [cm[-1]]))
    haplotypes = [f"{s}.{hap}" for s in samples for hap in (0, 1)]
    with open(args.output, "w") as f:
        codes_line = "\t".join(f"{a}={i}" for i, a in enumerate(args.ancestries))
        f.write(f"#Subpopulation order/codes: {codes_line}\n")
        f.write(
            "#chm\tspos\tepos\tsgpos\tegpos\tn snps\t" + "\t".join(haplotypes) + "\n"
        )
        for i, (start, end) in enumerate(zip(starts, ends)):
            fields = [
                chrom,
                pos[start],
                epos[i],
                f"{cm[start]:.2f}",
                f"{egpos[i]:.2f}",
                end - start,
            ]
            f.write("\t".join(map(str, fields)) + "\t")
            f.write("\t".join(map(str, codes[start])) + "\n")


if __name__ == "__main__":
    main()
