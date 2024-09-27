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

```
bcftools annotate \
 --annotations mane-annotations.bed.gz \
 --header-lines mane-annotations.hdr \
 --columns-file mane-annotations.cols \
 --output-type z \
 --output Manta_SAMPLE.candidateSV.vcf.gz \
 Manta_SAMPLE.candidateSV.vcf.gz
```
