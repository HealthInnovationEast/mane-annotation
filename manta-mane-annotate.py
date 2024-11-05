#!/usr/bin/env python3
"""
Custom annotation script to only annotate ends of events.

Additionally annotatess each element with low, high or both depending on status.
Only if low and high lists are identical is "both" used for any annotation.
"""
import sys

import click
import pysam


def _format_hits(tbx: pysam.TabixFile, rec: pysam.VariantRecord):
    low, high, hits = [], [], []
    for row in tbx.fetch(rec.chrom, rec.pos, rec.pos + 1):
        low.append(row[3])
    for row in tbx.fetch(rec.chrom, rec.stop, rec.stop + 1):
        high.append(row[3])
    if "".join(low) == "".join(high):
        for h in low:
            hits.append(f".|{h}")
    else:
        for h in low:
            hits.append(f"low|{h}")
        for h in high:
            hits.append(f"high|{h}")
    if hits:
        rec.info["AnnotMANE"] = ",".join(hits)


@click.command()
@click.option("--annotations", required=True, help="Mane annotation bed file (output of mane-prep.py)")
@click.option("--input", required=True, help="Input vcf (optionally compressed)")
@click.option("--output", required=True, help="Output vcf (index generated if vcf.gz)")
def annotate(annotations, input: str, output: str):
    """Annotates a VCF file with end specific events."""
    tbx = pysam.TabixFile(annotations, parser=pysam.asBed())
    contigs = tbx.contigs
    vcf_in = pysam.VariantFile(input)
    print("\tIndex not needed on input, ignore above warning", file=sys.stderr)
    vcf_in.header.add_line(
        '##INFO=<ID=AnnotMANE,Number=.,Type=String,Description="End|Transcript|ENSG|NCBI|AltName|Strand|ElementType|ElementNum">'
    )
    vcf_out = pysam.VariantFile(output, "w", header=vcf_in.header)

    for rec in vcf_in.fetch():
        if rec.chrom in contigs:
            _format_hits(tbx, rec)
        vcf_out.write(rec)
    vcf_out.close()
    if output.endswith(".gz"):
        pysam.tabix_index(output, preset="vcf", force=True)


if __name__ == "__main__":
    annotate()
