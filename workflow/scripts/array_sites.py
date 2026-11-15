"""Write the positions of an Illumina array manifest's autosomal SNPs.

The manifest is Illumina's sectioned CSV: a ``[Heading]`` block, the ``[Assay]``
table and a ``[Controls]`` block. The output is a headerless TSV of ``chr``-prefixed
chromosome and position, sorted and without repeats, which bcftools reads as a
targets file.
"""

import argparse
import csv

AUTOSOMES = {str(i) for i in range(1, 23)}
INDELS = {"[D/I]", "[I/D]"}


def build_parser():
    """Build the array_sites.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the array_sites.py flags.
    """
    p = argparse.ArgumentParser(
        description="Extract autosomal SNP positions from an Illumina manifest."
    )
    p.add_argument("--manifest", required=True, help="Illumina manifest CSV.")
    p.add_argument(
        "--output", required=True, help="Sites TSV to write: chromosome, position."
    )
    return p


def assay_rows(path):
    """Yield the rows of a manifest's ``[Assay]`` table.

    Parameters
    ----------
    path : str
        Illumina manifest CSV.

    Yields
    ------
    dict of str to str
        One assay row, keyed by the table's column names.
    """
    with open(path, newline="") as f:
        for line in f:
            if line.startswith("[Assay]"):
                break
        for row in csv.DictReader(f):
            if row["IlmnID"].startswith("["):
                return
            yield row


def main(argv=None):
    """Read the manifest and write its autosomal SNP positions.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    sites = {
        (int(row["Chr"]), int(row["MapInfo"]))
        for row in assay_rows(args.manifest)
        if row["Chr"] in AUTOSOMES
        and row["SNP"] not in INDELS
        and row["MapInfo"] != "0"
    }
    with open(args.output, "w") as f:
        f.writelines(f"chr{chrom}\t{pos}\n" for chrom, pos in sorted(sites))


if __name__ == "__main__":
    main()
