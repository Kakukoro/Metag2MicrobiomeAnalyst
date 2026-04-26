#!/usr/bin/env python3

import os
import glob
import pandas as pd
import re
from collections import defaultdict

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')
METADATA_FILE = os.path.join(BASE_DIR, 'metadata_bsg3.txt')
OUTPUT_READ_COUNT = os.path.join(BASE_DIR, 'read_count_table.csv')
OUTPUT_TAXONOMY = os.path.join(BASE_DIR, 'taxonomy_table.txt')


ALL_TAXONOMY_LEVELS = ['domain', 'phylum', 'class', 'subclass', 'order', 'suborder', 'family', 'genus', 'species', 'strain']

GROUPBY_LEVELS = ['domain', 'phylum', 'class', 'subclass', 'order', 'suborder', 'family', 'genus']

OUTPUT_COLS = ['#TAXONOMY', 'domain', 'phylum', 'class', 'subclass', 'order', 'suborder', 'family', 'genus']


def normalize_taxon_name(name):

    if pd.isna(name):
        return 'unclassified'
    
    name = str(name).strip()
    
    name = name.replace('&quot;', '"').replace('&amp;', '&')
    
    name = name.strip('"\'')
    

    unknown_terms = ['unclassified', 'unknown', 'uncultured', 'unidentified', 'incertae sedis', '']
    if name.lower() in unknown_terms:
        return 'unclassified'

    name = re.sub(r'^[^a-zA-Z0-9\-]+', '', name)
    name = re.sub(r'[^a-zA-Z0-9\-]+$', '', name)
    
    return name


def parse_vis_file(filepath):

    entries = []
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
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


def extract_taxonomy_and_otus():

    print("  Extracting taxonomy and OTU IDs from calc.VIS.txt files...")
    
    vis_pattern = os.path.join(DATA_DIR, 'metag*', 'calc.VIS.txt')
    all_files = glob.glob(vis_pattern)
    
    if not all_files:
        print("  ERROR: No calc.VIS.txt files found!")
        return None, None, None
    
    print(f"  Found {len(all_files)} calc.VIS.txt files")

    all_entries = []
    for f in all_files:
        entries = parse_vis_file(f)
        all_entries.extend(entries)
    
    print(f"  Total taxonomy entries: {len(all_entries)}")

    normalized_to_count = defaultdict(int)
    
    for entry in all_entries:
        count = entry[0]
        taxonomy = entry[1:]
        
        normalized_taxonomy = []
        for i, t in enumerate(taxonomy):
            if i < len(GROUPBY_LEVELS):
                normalized_taxonomy.append(normalize_taxon_name(t))
        
        normalized_taxonomy = normalized_taxonomy[:len(GROUPBY_LEVELS)]
        while len(normalized_taxonomy) < len(GROUPBY_LEVELS):
            normalized_taxonomy.append('unclassified')
        
        taxon_key = tuple(normalized_taxonomy)
        normalized_to_count[taxon_key] += count

    unique_taxonomies = sorted(normalized_to_count.keys())
    print(f"  Unique taxonomies (up to genus): {len(unique_taxonomies)}")
    
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


def create_wide_read_count_table(taxon_to_otu, metadata_sample_ids):

    print("  Creating wide read count table...")
    
    vis_pattern = os.path.join(DATA_DIR, 'metag*', 'calc.VIS.txt')
    all_files = glob.glob(vis_pattern)
    
    read_count_dict = defaultdict(lambda: defaultdict(int))

    for file_path in all_files:
        sample_id = os.path.basename(os.path.dirname(file_path))
        
        if sample_id not in metadata_sample_ids:
            continue
        
        entries = parse_vis_file(file_path)
        for entry in entries:
            count = entry[0]
            taxonomy = entry[1:]
            
            normalized_taxonomy = []
            for i, t in enumerate(taxonomy):
                if i < len(GROUPBY_LEVELS):
                    normalized_taxonomy.append(normalize_taxon_name(t))
            
            normalized_taxonomy = normalized_taxonomy[:len(GROUPBY_LEVELS)]
            while len(normalized_taxonomy) < len(GROUPBY_LEVELS):
                normalized_taxonomy.append('unclassified')
            
            taxon_key = tuple(normalized_taxonomy)
            
            if taxon_key in taxon_to_otu:
                otu_id = taxon_to_otu[taxon_key]
                read_count_dict[otu_id][sample_id] += count

    rows = []
    for otu_id in sorted(read_count_dict.keys()):
        row = {'#NAME': otu_id}
        for sample in metadata_sample_ids:
            row[sample] = read_count_dict[otu_id].get(sample, 0)
        rows.append(row)
    
    read_count_df = pd.DataFrame(rows)
    
    cols = ['#NAME'] + metadata_sample_ids
    read_count_df = read_count_df[cols]
    
    return read_count_df


def main():
    print("=" * 70)
    print("Combined Workflow - Processing calc.VIS.txt files")
    print("=" * 70)

    print("\n[Step 1] Loading metadata...")
    if os.path.exists(METADATA_FILE):
        metadata_df = pd.read_csv(METADATA_FILE, sep='\t', engine='python')
        metadata_sample_ids = metadata_df['#NAME'].tolist()
        print(f"  Loaded {len(metadata_sample_ids)} samples from metadata")
    else:
        print(f"  ERROR: Metadata file not found: {METADATA_FILE}")
        return
    
    print("\n[Step 2] Extracting taxonomy...")
    taxonomy_df, taxonomy_cols, taxon_to_otu = extract_taxonomy_and_otus()
    
    if taxonomy_df is None:
        print("  ERROR: Failed to extract taxonomy!")
        return
    
    print("\n[Step 3] Creating read count table...")
    read_count_df = create_wide_read_count_table(taxon_to_otu, metadata_sample_ids)
    
    print("\n[Step 4] Saving outputs...")
    taxonomy_df.to_csv(OUTPUT_TAXONOMY, sep='\t', index=False)
    read_count_df.to_csv(OUTPUT_READ_COUNT, index=False)
    
    print(f"\n" + "=" * 70)
    print("Workflow complete!")
    print("=" * 70)
    print(f"Taxonomy table: {OUTPUT_TAXONOMY}")
    print(f"  Shape: {taxonomy_df.shape}")
    print(f"  Columns: {list(taxonomy_df.columns)}")
    print(f"\nRead count table: {OUTPUT_READ_COUNT}")
    print(f"  Shape: {read_count_df.shape}")
    print(f"  Samples: {len(metadata_sample_ids)}")
    print(f"\nNote:")
    print(f"  - OTUs are grouped by GENUS level (species/strain differences ignored)")
    print(f"  - Unknown/unclassified normalized to 'unclassified'")
    print(f"  - Punctuation stripped from taxonomy names")
    print(f"  - Output files have '_combined_new' suffix to avoid overriding existing files")