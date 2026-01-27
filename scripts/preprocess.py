#!/usr/bin/env python3
"""Preprocess ZFIN ortholog files using DuckDB.

This script combines the fly, mouse, and human ortholog files into a single
TSV file with normalized gene IDs for the orthology transform.
"""

import sys
from pathlib import Path

import duckdb


def main():
    """Run the preprocessing SQL script with DuckDB."""
    # Get the project root directory (parent of scripts/)
    script_dir = Path(__file__).parent
    project_root = script_dir.parent

    # Path to the SQL script
    sql_file = script_dir / "preprocess.sql"

    if not sql_file.exists():
        print(f"Error: SQL file not found: {sql_file}", file=sys.stderr)
        sys.exit(1)

    # Check that required input files exist
    data_dir = project_root / "data"
    required_files = ["fly_orthos.txt", "human_orthos.txt", "mouse_orthos.txt"]

    for filename in required_files:
        filepath = data_dir / filename
        if not filepath.exists():
            print(f"Error: Required input file not found: {filepath}", file=sys.stderr)
            print("Run 'just download' first to download the source files.", file=sys.stderr)
            sys.exit(1)

    # Read the SQL script
    sql_content = sql_file.read_text()

    # Execute the SQL script using DuckDB
    # We need to change to the project root so relative paths work correctly
    import os

    original_dir = os.getcwd()
    try:
        os.chdir(project_root)
        con = duckdb.connect(":memory:")
        con.execute(sql_content)
        con.close()
        print(f"Successfully created {data_dir / 'zfin_orthologs.tsv'}")
    finally:
        os.chdir(original_dir)


if __name__ == "__main__":
    main()
