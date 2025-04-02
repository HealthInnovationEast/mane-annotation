# CHANGES

## 2.0.2

- Only output info fields if an annotation exists, better conformance with VCF spec

## 2.0.1

- Corrected header file generation (missing bp1/bp2)

## 2.0.0

- Dropped support for bedtools method due to restrictions for long events
- Preparation of files using `mane-prep.py` changed to support selection of MANE status.
  - Regeneration required: `*.bed.gz`, `*.bed.gz.tbi`, `*.hdr`
  - `*.cols` not needed or generated anymore.
- `manta-mane-annotate.py`
  - Split annotation of ends to separate INFO fields `AnnotMANEbp1` & `AnnotMANEbp2` to better support bedpe conversion.

## 1.0.0

Breakend behaviour is expanded.

## 0.2.0

- Adds `manta-mane-annotate.py` script to allow better control of how annotation is applied.
- Minor change to annotation string, see [README.md](README.md#annotation-format)

## 0.1.2

Documentation regarding tandem duplications and experimental bcftools option.

## 0.1.1

Documentation fix

## 0.1.0

Initial release
