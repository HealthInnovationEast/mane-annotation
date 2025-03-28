#!/usr/bin/env python3
"""
Custom annotation script to only annotate ends of events.

Additionally annotatess each element with low, high or both depending on status.
Only if low and high lists are identical is "both" used for any annotation.
"""
import re
import sys
from typing import List

import click
import pysam

bnd_re = re.compile("([\[\]])([^:]+):(\d+)[\[\]]")


def _breakend_coord(rec: pysam.VariantRecord):
    """
    Calculate the position of the related breakend.

    Breakend (BND) notation:
     - T[chr:pos[ → piece extending to the right of chr:pos is joined after T
     - T]chr:pos] → reverse complemented piece extending to the left of chr:pos is joined after T
     - ]chr:pos]T → piece extending to the left of chr:pos is joined before T
     - [chr:pos[T → reverse complemented piece extending to the right of chr:pos is joined before T

    Although we only technically care about the position we still need `][` to indicate if we should +/- 1
    """
    alt = rec.alts[0]
    (direct, chrom, pos) = bnd_re.search(alt).group(1, 2, 3)
    pos = int(pos)
    start, end = None, None
    if direct == "[":
        start = pos
        end = pos + 1
    else:
        start = pos - 1
        end = pos
    return chrom, start, end


def _end_hits(tbx: pysam.TabixFile, mane_type: str, chrom: str, start: int, end: int) -> List[str]:
    """
    Return list of hits for a coordinate set.

    Handles checking the contig is known to the tabix index.
    """
    hits = []
    if chrom in tbx.contigs:
        for row in tbx.fetch(chrom, start, end):
            if mane_type != "All":
                if mane_type not in row[3]:
                    continue
            hits.append(row[3])
    return hits


def _format_hits(tbx: pysam.TabixFile, mane_type: str, rec: pysam.VariantRecord):
    # setup bp2 coords
    (end_chrom, end_start) = (rec.chrom, rec.stop)
    end_stop = end_start + 1
    if "SVTYPE" in rec.info.keys() and rec.info["SVTYPE"] == "BND":
        (end_chrom, end_start, end_stop) = _breakend_coord(rec)
    # get hits
    bp1 = _end_hits(tbx, mane_type, rec.chrom, rec.pos, rec.pos + 1)
    bp2 = _end_hits(tbx, mane_type, end_chrom, end_start, end_stop)
    rec.info["AnnotMANEbp1"] = ",".join(bp1)
    rec.info["AnnotMANEbp2"] = ",".join(bp2)


def _expand_mane_type(mane_type: str) -> str:
    value = None
    if mane_type == "Clinical":
        value = "Plus_Clinical"
    else:
        value = mane_type
    return value


def _expand_header(annotations: str, vcf_head: pysam.VariantHeader):
    hdr_path = annotations.replace("bed.gz", "") + "hdr"
    with open(hdr_path, "rt") as hdr_in:
        info_lines = hdr_in.readlines()
    for l in info_lines:
        vcf_head.add_line(l)
    return vcf_head


@click.command()
@click.option("--annotations", required=True, help="Mane annotation bed file (output of mane-prep.py)")
@click.option("--input", required=True, help="Input vcf (optionally compressed)")
@click.option("--output", required=True, help="Output vcf (index generated if vcf.gz)")
@click.option(
    "--mode",
    default="All",
    help="Type of MANE annotation to include",
    type=click.Choice(["All", "Select", "Clinical"], case_sensitive=False),
)
def annotate(annotations, input: str, output: str, mode: str):
    """Annotates a VCF file with end specific events."""
    tbx = pysam.TabixFile(annotations, parser=pysam.asBed())
    contigs = tbx.contigs
    vcf_in = pysam.VariantFile(input)
    print("\tIndex not needed on input, ignore above warning", file=sys.stderr)

    ### This needs to be read from the name-annotations.hdr file
    vcf_header = _expand_header(annotations, vcf_in.header)
    vcf_out = pysam.VariantFile(output, "w", header=vcf_header)

    mane_type = _expand_mane_type(mode)
    for rec in vcf_in.fetch():
        if rec.chrom in contigs:
            _format_hits(tbx, mane_type, rec)
        vcf_out.write(rec)
    vcf_out.close()
    if output.endswith(".gz"):
        pysam.tabix_index(output, preset="vcf", force=True)


if __name__ == "__main__":
    annotate()
