#!/usr/bin/env python3
"""Download source data files for zfin-ingest.

This script reads download.yaml and fetches all required data files.
"""

import sys
from pathlib import Path
from urllib.request import urlretrieve

import yaml


def main():
    """Download all source data files."""
    # Get the project root directory (parent of scripts/)
    script_dir = Path(__file__).parent
    project_root = script_dir.parent

    # Path to the download config
    config_file = project_root / "download.yaml"

    if not config_file.exists():
        print(f"Error: Config file not found: {config_file}", file=sys.stderr)
        sys.exit(1)

    # Read the config
    with open(config_file) as f:
        config = yaml.safe_load(f)

    downloads = config.get("downloads", [])
    if not downloads:
        print("No downloads configured in download.yaml")
        return

    # Download each file
    for item in downloads:
        url = item["url"]
        local_name = item["local_name"]
        local_path = project_root / local_name

        # Create parent directory if needed
        local_path.parent.mkdir(parents=True, exist_ok=True)

        print(f"Downloading {url}")
        print(f"  -> {local_path}")

        try:
            urlretrieve(url, local_path)
            print(f"  Done ({local_path.stat().st_size:,} bytes)")
        except Exception as e:
            print(f"  Error: {e}", file=sys.stderr)
            sys.exit(1)

    print(f"\nAll {len(downloads)} files downloaded successfully.")


if __name__ == "__main__":
    main()
