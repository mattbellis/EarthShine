#!/usr/bin/env python3
"""
Script to combine parquet files locally.
Groups files by their base name (parameters) and combines all numbered chunks into single files.
"""

import os
import re
import shutil
from pathlib import Path
from collections import defaultdict
import pyarrow.parquet as pq
import pyarrow as pa


def get_base_name(filename):
    """Extract base name by removing the 6-digit number suffix."""
    match = re.match(r'(.+)_(\d{6})\.parquet$', filename)
    if match:
        return match.group(1)
    return None


def combine_parquet_files(input_dir, output_dir=None, move_to_subdir="NotCombined"):
    """
    Combine parquet files grouped by their base name.

    Args:
        input_dir: Directory containing the parquet files to combine
        output_dir: Directory to save combined files (defaults to input_dir)
        move_to_subdir: Subdirectory name to move original files into after combining
    """
    input_path = Path(input_dir)
    if output_dir is None:
        output_dir = input_path
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Create directory for moving original files
    if move_to_subdir:
        move_to_path = input_path / move_to_subdir
        move_to_path.mkdir(parents=True, exist_ok=True)
        print(f"Original files will be moved to: {move_to_path}")
    else:
        move_to_path = None

    # Group files by base name
    file_groups = defaultdict(list)

    print(f"Scanning directory: {input_path}")
    for file in sorted(input_path.glob("*.parquet")):
        base_name = get_base_name(file.name)
        if base_name:
            file_groups[base_name].append(file)

    print(f"\nFound {len(file_groups)} groups of files to combine:")
    for base_name, files in file_groups.items():
        print(f"  - {base_name}: {len(files)} files")

    # Combine each group
    for base_name, files in file_groups.items():
        print(f"\nCombining {len(files)} files for: {base_name}")

        # Read all files and combine
        tables = []
        for i, file in enumerate(files):
            try:
                table = pq.read_table(file)
                tables.append(table)
                if (i + 1) % 10 == 0:
                    print(f"  Loaded {i + 1}/{len(files)} files...")
            except Exception as e:
                print(f"  Error reading {file.name}: {e}")
                continue

        if not tables:
            print(f"  No valid tables to combine for {base_name}")
            continue

        # Concatenate all tables
        print(f"  Concatenating {len(tables)} tables...")
        combined_table = pa.concat_tables(tables)

        # Write combined file
        output_file = output_dir / f"{base_name}_combined.parquet"
        print(f"  Writing combined file: {output_file.name}")
        print(f"  Total rows: {len(combined_table)}")

        pq.write_table(combined_table, output_file)

        print(f"  ✓ Successfully created {output_file.name}")

        # Move original files to NotCombined directory
        if move_to_path:
            print(f"  Moving {len(files)} original files to {move_to_subdir}/...")
            for file in files:
                try:
                    shutil.move(str(file), str(move_to_path / file.name))
                except Exception as e:
                    print(f"    Warning: Could not move {file.name}: {e}")
            print(f"  ✓ Original files moved to {move_to_subdir}/")

    print(f"\n{'='*60}")
    print(f"All files combined successfully!")
    print(f"Combined files saved to: {output_dir}")


def main():
    input_dir = "/home/vamitamas/EarthShine/EarthShine/background_and_signal_acceptance_studies/OUTPUT_FILES"

    print("=" * 60)
    print("Parquet File Combiner")
    print("=" * 60)

    combine_parquet_files(input_dir)


if __name__ == "__main__":
    main()
