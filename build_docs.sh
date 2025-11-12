#!/bin/bash
# Build PyGEL3D documentation

echo "Building PyGEL3D documentation..."

# Check if mkdocs is installed
if ! command -v mkdocs &> /dev/null
then
    echo "mkdocs not found. Installing dependencies..."
    pip install -r doc/docs/requirements.txt
fi

# Build documentation
echo "Building documentation site..."
cd doc && mkdocs build

echo "Documentation built successfully!"
echo "Output is in the doc/site/ directory"
echo ""
echo "To view locally, run: cd doc && mkdocs serve"
echo "To deploy to GitHub Pages, run: cd doc && mkdocs gh-deploy"
