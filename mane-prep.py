#!/usr/bin/env python3
"""
Handles the full preparation of data files for use with bcftools.

Takes mane transcripts and gencode transcript information to generate an annotation bed file.
Additionally generated header and column files for the `bcftools annotate` command.
"""
import click
import pandas as pd
import pysam
from natsort import natsort_keygen

BED_HEADER = ("#CHROM", "BEG", "END", "AnnotMANE")

INFO_LINE = '##INFO=<ID=AnnotMANE,Number=.,Type=String,Description="End|Transcript|ENSG|NCBI|AltName|ManeType|Strand|ElementType|ElementNum">'

COLS = ("CHROM", "BEG", "END", "AnnotMANE")


def compression(filepath):
    """Identify if file is compressed and return suitable value for pandas compression attribute."""
    with open(filepath, "rb") as test_f:
        if test_f.read(2) == b"\x1f\x8b":
            return "gzip"
    return None


def sorted_df(data):
    """Handle the conversion of array to data frame and sorting."""
    df = pd.DataFrame(data, columns=BED_HEADER)
    df["BEG"] = df["BEG"].astype(int)
    df["END"] = df["END"].astype(int)
    df.sort_values(["#CHROM", "BEG", "END"], ascending=[True, True, True], inplace=True, key=natsort_keygen())
    return df


def to_tabix(outfile, frames):
    """Write multiple dataframes to file appending as appropriate then compresses and indexes with tabix."""
    mode = "w"
    for df in frames:
        df.to_csv(outfile, sep="\t", index=False, header=True, mode=mode)
        mode = "a"
    pysam.tabix_index(outfile, preset="bed", force=True)


def annot_files(outfile):
    """Generate the header and column input files for bcftools."""
    with open(f"{outfile}.hdr", "w") as ofh:
        print(INFO_LINE, file=ofh)
    with open(f"{outfile}.cols", "w") as ofh:
        for l in COLS:
            print(l, file=ofh)


@click.command()
@click.option("--mane", required=True, help="Mane transcript list from UCSC Table Browser (opt gz)")
@click.option("--gencode", required=True, help="Gencode comprehensive from UCSC Table Browser (opt gz)")
@click.option("--outstub", required=True, help="Path prefix for output, no file extensions, folders must exist")
def mane(mane, gencode, outstub):
    """Driver function to generate annotation ready files from UCSC table outputs."""
    df_mane = pd.read_csv(mane, sep="\t", compression=compression(mane))
    cols_mane = ["name", "geneName", "geneName2", "ncbiGene", "maneStat"]
    df_mane = df_mane[cols_mane]
    df_gencode = pd.read_csv(gencode, sep="\t", compression=compression(gencode))
    cols_gencode = ["name", "chrom", "strand", "exonCount", "exonStarts", "exonEnds"]
    df_gencode = df_gencode[cols_gencode]

    intersected_df = pd.merge(df_mane, df_gencode, on="name", how="inner")

    """
    Need to split each record by the exonStarts/exonEnds, Calculate the introns, handle the -ve strand numbering
    """
    result_core = []
    result_alts = []
    for row in intersected_df.itertuples():
        exonStarts = [int(x) for x in row.exonStarts.split(",") if x != ""]
        exonEnds = [int(x) for x in row.exonEnds.split(",") if x != ""]

        blocks = []
        for i in range(int(row.exonCount)):
            exon = i + 1
            if row.strand == "-":
                exon = row.exonCount - i
            blocks.append(
                {
                    "type": "exon",
                    "num": exon,
                    "start": exonStarts[i],
                    "end": exonEnds[i],
                }
            )

            if i == row.exonCount - 1:
                break

            intron = i + 1
            if row.strand == "-":
                intron = row.exonCount - (1 + i)
            blocks.append(
                {
                    "type": "intron",
                    "num": intron,
                    "start": exonEnds[i] + 1,
                    "end": exonStarts[i + 1] - 1,
                }
            )

        attributes = [row.name, row.geneName, row.ncbiGene, row.geneName2, row.maneStat]

        for ent in blocks:
            annotation = "|".join([*attributes, row.strand, str(ent["type"]), str(ent["num"])])
            item = [row.chrom, str(ent["start"] - 1), str(ent["end"]), annotation]  # half-open bed style output
            if "_" in row.chrom:
                result_alts.append(item)
            else:
                result_core.append(item)

    to_tabix(f"{outstub}.bed", [sorted_df(result_core), sorted_df(result_alts)])
    annot_files(outstub)


if __name__ == "__main__":
    mane()
