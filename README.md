# Metag2MicrobiomeAnalyst

A Python script for processing metagenomic `calc.VIS.txt` files to generate taxonomy and read count tables for downstream analysis.

## Overview

This workflow processes taxonomic classification results from metagenomic sequencing data. It extracts taxonomy information from `calc.VIS.txt` files, normalizes taxonomic names, creates Operational Taxonomic Unit (OTU) mappings, and generates both a taxonomy reference table and a sample-by-OTU read count matrix.

## Workflow Diagram

```
Input Files
├── data/
│   ├── metag1/
│   │   └── calc.VIS.txt
│   ├── metag2/
│   │   └── calc.VIS.txt
│   └── ...
│
└── metadata_bsg3.txt

    ↓

metag2microbiomeanalyst.py
├── Step 1: Load metadata (extract sample IDs)
├── Step 2: Scan for calc.VIS.txt files
├── Step 3: Parse and normalize taxonomy
├── Step 4: Create OTU mappings (grouped at genus level)
├── Step 5: Build read count table
└── Step 6: Save output files

    ↓

Output Files
├── output_taxonomy_table.txt
└── output_read_count_table.csv
```

## Requirements

- Python 3.7+
- pandas
- Required packages: `pip install pandas`

## Usage

### Basic Command

```bash
python metag2microbiomeanalyst.py \
    --input-dir DATA_DIRECTORY \
    --metadata METADATA_FILE \
    --output-prefix OUTPUT_PREFIX
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--input-dir` | Directory containing sample subdirectories with `calc.VIS.txt` files |
| `--metadata` | Path to metadata file (tab-separated, must contain `#NAME` column) |
| `--output-prefix` | Prefix for output files. Two files will be created: `<prefix>_taxonomy_table.txt` and `<prefix>_read_count_table.csv` |

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--sample-pattern` | `metag*` | Pattern to match sample subdirectories. Change if your files use a different naming convention (e.g., `sample*`) |

### Example

```bash
# Process files in ./data directory with default metag* pattern
python metag2microbiomeanalyst.py \
    --input-dir ./data \
    --metadata metadata_bsg3.txt \
    --output-prefix results/analysis_2024

# Output files:
# - results/analysis_2024_taxonomy_table.txt
# - results/analysis_2024_read_count_table.csv

# Using a custom sample pattern
python metag2microbiomeanalyst.py \
    --input-dir ./samples \
    --metadata metadata.txt \
    --output-prefix output/results \
    --sample-pattern sample_*
```

## Input File Formats

### calc.VIS.txt Files

Individual sample classification results from metagenomic analysis. Each file is tab-separated with:

- **Column 0**: Read count (integer)
- **Column 1-10**: Taxonomy levels (domain, phylum, class, subclass, order, suborder, family, genus, species, strain)

Example:
```
100	Bacteria	Firmicutes	Bacilli	Bacillales	Bacillaceae	Bacillus	Bacillus_subtilis	strain1
50	Bacteria	Proteobacteria	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Escherichia	Escherichia_coli	strain2
```

### Metadata File

Tab-separated file containing sample information. **Must contain a `#NAME` column** with sample IDs that match the subdirectory names containing the `calc.VIS.txt` files.

Example (`metadata_bsg3.txt`):
```
#NAME	patient_ID	isolation_code	sex	time_point	...
metag1_bsg3	K10	602	2	2.5y	...
metag2_bsg3	K11	603	1	3.0y	...
```

In this example, the script will look for `calc.VIS.txt` files in directories named `metag1_bsg3`, `metag2_bsg3`, etc.

## Output Files

### Taxonomy Table (`*_taxonomy_table.txt`)

Tab-separated file with one row per unique OTU (grouped at genus level).

| Column | Description |
|--------|-------------|
| `#TAXONOMY` | OTU ID (format: OTUXXXXXX) |
| `domain` | Normalized domain name |
| `phylum` | Normalized phylum name |
| `class` | Normalized class name |
| `subclass` | Normalized subclass name |
| `order` | Normalized order name |
| `suborder` | Normalized suborder name |
| `family` | Normalized family name |
| `genus` | Normalized genus name |

### Read Count Table (`*_read_count_table.csv`)

Comma-separated file in wide format with samples as columns and OTUs as rows.

| Column | Description |
|--------|-------------|
| `#NAME` | OTU ID |
| `sample1` | Read count for sample1 |
| `sample2` | Read count for sample2 |
| ... | ... |

## Data Processing Details

### Taxonomy Normalization

The script performs several normalization steps on taxonomy names:

1. **Missing values**: Converted to `'unclassified'`
2. **HTML entities**: `&quot;` → `"`, `&amp;` → `&`
3. **Surrounding quotes**: Stripped from start and end
4. **Known unknown terms**: The following are standardized to `'unclassified'`:
   - `unclassified`
   - `unknown`
   - `uncultured`
   - `unidentified`
   - `incertae sedis`
   - Empty strings
5. **Punctuation**: Non-alphanumeric characters stripped from start and end

### OTU Grouping

- OTUs are **grouped at the genus level**
- Species and strain differences are **ignored**
- All reads matching the same taxonomy up to genus are aggregated into a single OTU
- Each unique genus-level taxonomy gets a unique OTU ID (e.g., `OTU000001`, `OTU000002`)

### Sample Matching

- The script extracts sample IDs from the **directory names** containing `calc.VIS.txt` files
- Only samples that exist in the metadata `#NAME` column are included in the output
- Sample order in the output follows the order in the metadata file.
