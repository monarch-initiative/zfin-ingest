# zfin-ingest justfile

# Explicitly enumerate transforms (add new ingests here)
TRANSFORMS := "gene_to_phenotype genotype_to_phenotype orthology"

# List all commands
_default:
    @just --list

# Initialize a new project
[group('project management')]
setup: _git-init install _git-add
    git commit -m "Initialize zfin-ingest"

# Install dependencies
[group('project management')]
install:
    uv sync --group dev

# Download source data
[group('ingest')]
download:
    uv run python scripts/download.py

# Preprocess orthology data (combine fly, mouse, human orthologs)
[group('ingest')]
preprocess:
    uv run python scripts/preprocess.py

# Run all transforms
[group('ingest')]
transform-all: preprocess
    #!/usr/bin/env bash
    set -euo pipefail
    for t in {{TRANSFORMS}}; do
        if [ -n "$t" ]; then
            echo "Transforming $t..."
            uv run koza transform $t.yaml
        fi
    done

# Run specific transform
[group('ingest')]
transform NAME:
    uv run koza transform {{NAME}}.yaml

# Run tests
[group('development')]
test:
    uv run pytest

# Run tests with coverage
[group('development')]
test-cov:
    uv run pytest --cov=. --cov-report=term-missing

# Lint code
[group('development')]
lint:
    uv run ruff check .

# Format code
[group('development')]
format:
    uv run ruff format .

# Clean output directory
[group('ingest')]
clean:
    rm -rf output/

# Hidden recipes
_git-init:
    git init

_git-add:
    git add .
