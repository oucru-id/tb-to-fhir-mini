# TB to FHIR Genomics VCF/Deeplex Excel Pipeline

A Nextflow pipeline for converting Mycobacterium tuberculosis VCF and Deeplex Excel genomic data to FHIR Genomics format with drug resistance classification and lineage analysis.

## Overview

This pipeline processes TB genomic data from two sources:
- **VCF files**: Variant call format files from sequencing platforms
- **Deeplex Excel files**: Post-analysis TB mutation data from the Deeplex platform

### Key Features

- VCF normalization and variant filtering
- TB lineage classification (barcode SNP analysis)
- Drug resistance prediction and classification
- FHIR Bundle generation
- FHIR validation

## Requirements

### System Dependencies

- **Nextflow** ≥ 20.10.0
- **Java** 11+
- **Python** 3.8+
- `bcftools` ≥ 1.14
- `samtools` ≥ 1.14
- `tabix` ≥ 1.11
- `bgzip` ≥ 1.11

## Installation

From the repo
```bash
git clone https://github.com/oucru-id/tb-to-fhir-mini.git
cd tb-to-fhir-mini
```
Docker installation
```bash
docker pull ghcr.io/javiadividya/tb-to-fhir-mini:1.0
```

## Directory Structure

```
tb-to-fhir-vcfdeeplex/
├── main.nf                          # Main workflow
├── nextflow.config                  # Configuration
├── workflows/
│   ├── vcf.nf                       # VCF processing
│   ├── lineage.nf                   # Lineage classification
│   ├── fhir.nf                      # FHIR conversion
│   ├── deeplex.nf                   # Deeplex processing
│   ├── merge_clinical_data.nf       # Clinical metadata merge
│   ├── validate_fhir.nf             # FHIR validation
│   ├── report.nf                    # Report generation
│   └── utils.nf                     # Utility functions
├── scripts/
│   ├── annotated_to_fhir.py         # VCF to FHIR converter
│   ├── xlsx_json_converter.py       # Deeplex to FHIR converter
│   ├── merge_clinical_fhir.py       # Clinical data merge (VCF)
│   ├── merge_clinical_deeplex.py    # Clinical data merge (Deeplex)
│   ├── lineage_classifier.py        # TB lineage classification
│   ├── generate_sample_report.py    # Sample report generation
│   ├── clinical_metadata_parser.py  # Metadata parser
│   └── get_versions.py              # Version collection
└── data/
    ├── H37Rv.fasta                  # M. tuberculosis reference genome
    ├── WHO-UCN-TB-2023.csv          # TB mutation database
    ├── lineage_fix.bed              # Lineage barcode SNP positions
    ├── repetitive_regions.bed       # Repetitive region sequences
    ├── clinical_metadata.csv        # VCF sample metadata
    └── clinical_metadata_deeplex.csv # Deeplex sample metadata
```

## Input Data

### VCF Files

Place VCF files in the `data/VCF/` directory. Supported formats:
- `.vcf` (uncompressed)
- `.vcf.gz` (gzip compressed)

```bash
data/VCF/
├── sample_001.vcf.gz
├── sample_002.vcf.gz
└── sample_003.vcf
```

### Deeplex Files

Place Excel files in `data/Deeplex/` directory:

```bash
data/Deeplex/
├── deeplex_batch_001.xlsx
└── deeplex_batch_002.xlsx
```

## Usage

### Basic Run

```bash
nextflow run main.nf
```

### With Custom Paths

```bash
nextflow run main.nf \
    --vcf_dir /path/to/vcf \
    --deeplex_dir /path/to/deeplex \
    --results_dir /path/to/results \
```

### Resume Failed Runs

```bash
nextflow run main.nf -resume
```

## Workflow Steps

### VCF Processing

1. **VCF_NORMALIZE**: Renaming, sorting, and normalization
2. **filter_variants**: Exclude repetitive regions, filter by variant type, depth, and quality
3. **annotate**: Add functional annotations from mutation database
4. **LINEAGE**: Extract barcode SNPs and classify TB lineage
5. **FHIR**: Convert to FHIR Observation bundles
6. **MERGE_CLINICAL_DATA**: Integrate clinical metadata
7. **VALIDATE**: Validate FHIR against profiles

### Deeplex Processing

1. **CONVERT_DEEPLEX**: Extract data from Excel files
2. **MERGE_CLINICAL_DEEPLEX**: Add clinical metadata
3. **VALIDATE**: Validate FHIR output

## Output Structure

```
results/
├── fhir/                              # VCF-derived FHIR
│   ├── sample_001.fhir.json
│   └── sample_002.fhir.json
├── fhir_merged/                       # Clinical data merged VCF FHIR Bundles
│   ├── sample_001.merged.fhir.json
│   └── sample_002.merged.fhir.json
├── fhir_deeplex/                      # Deeplex-derived FHIR
│   └── deeplex_batch_001.json
├── fhir_deeplex_merged/               # Clinical data merged Deeplex FHIR Bundles
│   └── deeplex_batch_001.merged.fhir.json
├── lineage/                           # Lineage classification results
│   ├── sample_001.lineage.json
│   └── sample_002.lineage.json
├── reports/                           # Sample summary reports
│   ├── sample_001.summary_report.txt
│   └── sample_002.summary_report.txt
├── validated_fhir/                    # Validation results
│   ├── sample_001.validation_report.json
│   └── sample_002.validation_report.json
└── runningstat/                       # Nextflow execution reports
    ├── execution.html
    ├── timeline.html
    └── dag.html
```

## Support
[GitHub Issues](https://github.com/oucru-id/tb-to-fhir-mini/issues)
