#!/usr/bin/env python3
"""
Process calc.VIS.txt files from metagenomic data to generate taxonomy and read count tables.

This script:
1. Scans a directory for calc.VIS.txt files (typically in metag* subdirectories)
2. Parses taxonomy information from these files
3. Normalizes taxonomy names (handling unclassified, punctuation, etc.)
4. Creates OTU (Operational Taxonomic Unit) mappings, for easier parsing
5. Generates a wide read count table using sample IDs from metadata
6. Outputs taxonomy table and read count table to specified files

Usage:
    python process_vis_files.py --input-dir DATA_DIR --metadata METADATA_FILE --output-prefix OUTPUT_PREFIX

Example:
    python process_vis_files.py --input-dir ./data --metadata metadata.txt --output-prefix results/my_output
    # Generates: results/my_output_taxonomy_table.txt and results/my_output_read_count_table.csv
"""

import argparse
import os
import glob
import pandas as pd
import re
from collections import defaultdict

# Taxonomy levels used for grouping and output
ALL_TAXONOMY_LEVELS = ['domain', 'phylum', 'class', 'subclass', 'order', 'suborder', 'family', 'genus', 'species', 'strain']

# Levels to group by (up to genus level)
GROUPBY_LEVELS = ['domain', 'phylum', 'class', 'subclass', 'order', 'suborder', 'family', 'genus']

# Output columns for taxonomy table
OUTPUT_COLS = ['#TAXONOMY', 'domain', 'phylum', 'class', 'subclass', 'order', 'suborder', 'family', 'genus']


def normalize_taxon_name(name):
    """
    Normalize a taxonomy level name for consistency.
    
    Handles:
    - Missing/NaN values -> 'unclassified'
    - HTML entities (&quot;, &amp;) -> proper characters
    - Strips surrounding quotes and apostrophes
    - Known unknown terms -> 'unclassified'
    - Strips non-alphanumeric characters from start/end
    
    Args:
        name: Raw taxonomy name string (may be NaN)
        
    Returns:
        Normalized taxonomy name string
    """
    if pd.isna(name):
        return 'unclassified'
    
    name = str(name).strip()
    
    # Replace HTML entities
    name = name.replace('&quot;', '"').replace('&amp;', '&')
    
    # Strip surrounding quotes and apostrophes
    name = name.strip('"\'')
    
    # Map known unknown terms to standardized 'unclassified'
    unknown_terms = ['unclassified', 'unknown', 'uncultured', 'unidentified', 'incertae sedis', '']
    if name.lower() in unknown_terms:
        return 'unclassified'

    # Strip non-alphanumeric characters from start and end
    name = re.sub(r'^[^a-zA-Z0-9\-]+', '', name)
    name = re.sub(r'[^a-zA-Z0-9\-]+$', '', name)
    
    return name


def parse_vis_file(filepath):
    """
    Parse a calc.VIS.txt file and extract count and taxonomy entries.
    
    File format: Each line is tab-separated with:
    - Column 0: Read count
    - Columns 1-10: Taxonomy levels (domain through strain)
    
    Args:
        filepath: Path to calc.VIS.txt file
        
    Returns:
        List of tuples: (count, taxonomy_level_1, ..., taxonomy_level_10)
    """
    entries = []
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            # Pad with empty strings if line has fewer than 11 columns
            while len(parts) < 11:
                parts.append('')
            count = parts[0]
            taxonomy = parts[1:11]
            try:
                count_val = int(count)
            except ValueError:
                count_val = 1
            entries.append((count_val,) + tuple(taxonomy))
    return entries


def extract_taxonomy_and_otus(input_dir, sample_pattern='metag*'):
    """
    Extract taxonomy information and create OTU mappings from calc.VIS.txt files.
    
    Scans the input directory for subdirectories matching the sample pattern,
    then finds calc.VIS.txt files within those subdirectories.
    
    Args:
        input_dir: Base directory to search for calc.VIS.txt files
        sample_pattern: Pattern to match sample subdirectories (default: 'metag*')
        
    Returns:
        tuple: (taxonomy_df, taxonomy_cols, taxon_to_otu)
        - taxonomy_df: DataFrame with taxonomy table
        - taxonomy_cols: List of taxonomy level names
        - taxon_to_otu: Dict mapping taxonomy tuples to OTU IDs
    """
    print("  Extracting taxonomy and OTU IDs from calc.VIS.txt files...")
    
    # Build pattern to find calc.VIS.txt files in sample subdirectories
    vis_pattern = os.path.join(input_dir, sample_pattern, 'calc.VIS.txt')
    all_files = glob.glob(vis_pattern)
    
    if not all_files:
        print(f"  ERROR: No calc.VIS.txt files found matching pattern: {vis_pattern}")
        return None, None, None
    
    print(f"  Found {len(all_files)} calc.VIS.txt files")

    # Parse all entries from all files
    all_entries = []
    for f in all_files:
        entries = parse_vis_file(f)
        all_entries.extend(entries)
    
    print(f"  Total taxonomy entries: {len(all_entries)}")

    # Normalize taxonomy and aggregate counts
    normalized_to_count = defaultdict(int)
    
    for entry in all_entries:
        count = entry[0]
        taxonomy = entry[1:]
        
        # Normalize each taxonomy level
        normalized_taxonomy = []
        for i, t in enumerate(taxonomy):
            if i < len(GROUPBY_LEVELS):
                normalized_taxonomy.append(normalize_taxon_name(t))
        
        # Ensure we have exactly len(GROUPBY_LEVELS) levels
        normalized_taxonomy = normalized_taxonomy[:len(GROUPBY_LEVELS)]
        while len(normalized_taxonomy) < len(GROUPBY_LEVELS):
            normalized_taxonomy.append('unclassified')
        
        taxon_key = tuple(normalized_taxonomy)
        normalized_to_count[taxon_key] += count

    unique_taxonomies = sorted(normalized_to_count.keys())
    print(f"  Unique taxonomies (up to genus): {len(unique_taxonomies)}")
    
    # Create OTU mappings and taxonomy table
    taxon_to_otu = {}
    taxonomy_rows = []
    
    for idx, taxon_key in enumerate(unique_taxonomies, 1):
        otu_id = f'OTU{idx:06d}'
        taxon_to_otu[taxon_key] = otu_id
        
        row_dict = {'#TAXONOMY': otu_id}
        for i, level in enumerate(OUTPUT_COLS[1:]):
            if i < len(taxon_key):
                val = taxon_key[i]
                row_dict[level] = val if val != '' else ''
            else:
                row_dict[level] = ''
        taxonomy_rows.append(row_dict)
    
    taxonomy_df = pd.DataFrame(taxonomy_rows)
    
    return taxonomy_df, GROUPBY_LEVELS, taxon_to_otu


def create_wide_read_count_table(taxon_to_otu, metadata_sample_ids, input_dir, sample_pattern='metag*'):
    """
    Create a wide-format read count table with samples as columns and OTUs as rows.
    
    Args:
        taxon_to_otu: Dict mapping taxonomy tuples to OTU IDs
        metadata_sample_ids: List of sample IDs from metadata file
        input_dir: Base directory containing calc.VIS.txt files
        sample_pattern: Pattern to match sample subdirectories
        
    Returns:
        DataFrame with read counts (samples as columns, OTUs as rows)
    """
    print("  Creating wide read count table...")
    
    # Find all calc.VIS.txt files
    vis_pattern = os.path.join(input_dir, sample_pattern, 'calc.VIS.txt')
    all_files = glob.glob(vis_pattern)
    
    # Dictionary to accumulate read counts: {otu_id: {sample_id: count}}
    read_count_dict = defaultdict(lambda: defaultdict(int))

    for file_path in all_files:
        # Extract sample ID from directory name
        sample_id = os.path.basename(os.path.dirname(file_path))
        
        # Only process samples that exist in metadata
        if sample_id not in metadata_sample_ids:
            continue
        
        entries = parse_vis_file(file_path)
        for entry in entries:
            count = entry[0]
            taxonomy = entry[1:]
            
            # Normalize taxonomy
            normalized_taxonomy = []
            for i, t in enumerate(taxonomy):
                if i < len(GROUPBY_LEVELS):
                    normalized_taxonomy.append(normalize_taxon_name(t))
            
            normalized_taxonomy = normalized_taxonomy[:len(GROUPBY_LEVELS)]
            while len(normalized_taxonomy) < len(GROUPBY_LEVELS):
                normalized_taxonomy.append('unclassified')
            
            taxon_key = tuple(normalized_taxonomy)
            
            # If this taxonomy exists in our OTU mapping, add to count
            if taxon_key in taxon_to_otu:
                otu_id = taxon_to_otu[taxon_key]
                read_count_dict[otu_id][sample_id] += count

    # Build DataFrame rows
    rows = []
    for otu_id in sorted(read_count_dict.keys()):
        row = {'#NAME': otu_id}
        for sample in metadata_sample_ids:
            row[sample] = read_count_dict[otu_id].get(sample, 0)
        rows.append(row)
    
    read_count_df = pd.DataFrame(rows)
    
    # Ensure column order: #NAME first, then samples in metadata order
    cols = ['#NAME'] + metadata_sample_ids
    read_count_df = read_count_df[cols]
    
    return read_count_df


def main():
    """Main workflow function."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Process calc.VIS.txt files to generate taxonomy and read count tables.'
    )
    parser.add_argument(
        '--input-dir',
        required=True,
        help='Directory containing metag* subdirectories with calc.VIS.txt files'
    )
    parser.add_argument(
        '--metadata',
        required=True,
        help='Path to metadata file (tab-separated, must contain #NAME column)'
    )
    parser.add_argument(
        '--output-prefix',
        required=True,
        help='Output file prefix (e.g., "results/output" generates output_taxonomy_table.txt and output_read_count_table.csv)'
    )
    parser.add_argument(
        '--sample-pattern',
        default='metag*',
        help='Pattern to match sample subdirectories (default: metag*)'
    )
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.isdir(args.input_dir):
        print(f"ERROR: Input directory does not exist: {args.input_dir}")
        return
    
    # Validate metadata file
    if not os.path.exists(args.metadata):
        print(f"ERROR: Metadata file does not exist: {args.metadata}")
        return
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print("=" * 70)
    print("VIS File Processor - Processing calc.VIS.txt files")
    print("=" * 70)

    # Step 1: Load metadata
    print("\n[Step 1] Loading metadata...")
    try:
        metadata_df = pd.read_csv(args.metadata, sep='\t', engine='python')
        metadata_sample_ids = metadata_df['#NAME'].tolist()
        print(f"  Loaded {len(metadata_sample_ids)} samples from metadata")
    except Exception as e:
        print(f"  ERROR: Failed to load metadata: {e}")
        return
    
    # Step 2: Extract taxonomy
    print("\n[Step 2] Extracting taxonomy...")
    taxonomy_df, taxonomy_cols, taxon_to_otu = extract_taxonomy_and_otus(
        args.input_dir, args.sample_pattern
    )
    
    if taxonomy_df is None:
        print("  ERROR: Failed to extract taxonomy!")
        return
    
    # Step 3: Create read count table
    print("\n[Step 3] Creating read count table...")
    read_count_df = create_wide_read_count_table(
        taxon_to_otu, metadata_sample_ids, args.input_dir, args.sample_pattern
    )
    
    # Step 4: Save outputs
    print("\n[Step 4] Saving outputs...")
    
    taxonomy_output = f"{args.output_prefix}_taxonomy_table.txt"
    read_count_output = f"{args.output_prefix}_read_count_table.csv"
    
    taxonomy_df.to_csv(taxonomy_output, sep='\t', index=False)
    read_count_df.to_csv(read_count_output, index=False)
    
    print(f"\n" + "=" * 70)
    print("Workflow complete!")
    print("=" * 70)
    print(f"Taxonomy table: {taxonomy_output}")
    print(f"  Shape: {taxonomy_df.shape}")
    print(f"  Columns: {list(taxonomy_df.columns)}")
    print(f"\nRead count table: {read_count_output}")
    print(f"  Shape: {read_count_df.shape}")
    print(f"  Samples: {len(metadata_sample_ids)}")
    print(f"\nNote:")
    print(f"  - OTUs are grouped by GENUS level (species/strain differences ignored)")
    print(f"  - Unknown/unclassified normalized to 'unclassified'")
    print(f"  - Punctuation stripped from taxonomy names")


if __name__ == '__main__':
    main()
