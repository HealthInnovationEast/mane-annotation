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

See `--help` for full command details.

## Annotate VCF

The included Manta specific annotation script enures only annotation spanning the ends of an event are included.

### `manta-mane-annotate.py`

```
./manta-mane-annotate.py \
 --annotations mane-annotations.bed.gz \
 --input input.vcf.gz \
 --output annotated.vcf.gz
 --mode all
```

See `--help` for full command details.

#### `--mode`

Allows for limiting the types of MANE annotation that will be applied, selected via `maneStatus` of the table browser output.

## Annotation format

The annotation is appended to the `INFO` field with the tags `AnnotMANEbp1`, `AnnotMANEbp2`.  Both are required to handle
long events.  In the case of `BND` calls the `bp2` entry of record `A` will correspond to `bp1` of record `B`
and vice-versa.

Following the format is described in the header INFO line:

```
##INFO=<ID=AnnotMANEbp1,Number=.,Type=String,Description="Transcript|ENSG|NCBI|AltName|Strand|ElementType|ElementNum">
##INFO=<ID=AnnotMANEbp2,Number=.,Type=String,Description="Transcript|ENSG|NCBI|AltName|Strand|ElementType|ElementNum">
```

Where:

| Field       | Description       | Examples                 |
| ----------- | ----------------- | ------------------------ |
| Transcript  | ENST value        | `ENST00000303635.12`     |
| ENSG        | ENSG gene name    | `ENSG00000171735.21`     |
| NCBI        | NCBI gene name/ID | `GeneID:23261`           |
| AltName     | Colloquial name   | `CAMTA1`                 |
| ManeType    | MANE status       | `Select`, `Plus_Clinical | | Strand      | Coding strand     | `+`or`-`              | | ElementType | Type of feature   |`intron`or`exon`      | | ElementNum  | Feature number    |`1`..`N\`                 |
