# mane-annotation

Preparation scripts for mane annotation using bcftools

## Setup

```
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Input files

Generate the following files via the USCS Table Browser: https://genome.ucsc.edu/cgi-bin/hgTables

Fields not specified should be cleared:

- mane-hg38.tsv.gz
  - Clade = Mammal
  - Genome = Human
  - Assembly = Dec. 2013 (GRCh38/hg38)
  - Group = Genes and Gene Predictions
  - Track = MANE
  - Table = mane
  - Region = Genome
  - Output format = All fields from selected table
  - Output filename = mane-hg38.tsv.gz
  - Output field separator = tsv
  - File type returned = Gzip compressed
- gencode-hg38.tsv.gz
  - Clade = Mammal
  - Genome = Human
  - Assembly = Dec. 2013 (GRCh38/hg38)
  - Group = Genes and Gene Predictions
  - Track = All GENCODE V46
  - Table = Comprehensive (wgEncodeGencodeCompV46)
  - Region = Genome
  - Output format = All fields from selected table
  - Output filename = gencode-hg38.tsv.gz
  - Output field separator = tsv
  - File type returned = Gzip compressed

## Prepare annotations

```
./mane-prep.py --mane mane-hg38.tsv.gz --gencode gencode-hg38.tsv.gz --outstub mane-annotations
```

## Annotate VCF

The included Manta specific annotation script enures only annotation spanning the ends of an event are included.

### `manta-mane-annotate.py`

```
./manta-mane-annotate.py \
 --annotations mane-annotations.bed.gz \
 --input input.vcf.gz \
 --output annotated.vcf.gz
```

### `bcftools annotate`

```
bcftools annotate \
 --annotations mane-annotations.bed.gz \
 --header-lines mane-annotations.hdr \
 --columns-file mane-annotations.cols \
 --output-type z \
 --output annotated.vcf.gz \
 input.vcf.gz
```

Depending on type of variants being processed you may also want to include:

```
 --merge-logic AnnotMANE:unique
```

This is not recommended for long variants such as tandem duplications.  You would need to split those out of the input
to prevent very long records/bloat.

NOTE: `--merge-logic` is experimental, default behaviour is the "first hit".

## Annotation format

The annotation is appended to the `INFO` field with the tag `AnnotMANE`:

Following the format is described in the header INFO line:

```
##INFO=<ID=AnnotMANE,Number=.,Type=String,Description="End|Transcript|ENSG|NCBI|AltName|Strand|ElementType|ElementNum">
```

Where:

| Field            | Description                        | Examples             |
| ---------------- | ---------------------------------- | -------------------- |
| End <sup>+</sup> | Which end of event this applies to | `low`, `high`, `.`   |
| Transcript       | ENST value                         | `ENST00000303635.12` |
| ENSG             | ENSG gene name                     | `ENSG00000171735.21` |
| NCBI             | NCBI gene name/ID                  | `GeneID:23261`       |
| AltName          | Colloquial name                    | `CAMTA1`             |
| Strand           | Coding strand                      | `+` or `-`           |
| ElementType      | Type of feature                    | `intron` or `exon`   |
| ElementNum       | Feature number                     | `1`..`N`             |

<sup>+</sup> Only indicates low/high when annotating events with `manta-mane-annotate.py`, where same annotation for both
`.` is expected.
