# PyGEL3D Documentation

This directory contains the MkDocs documentation for PyGEL3D.

## Building the Documentation

### Install Dependencies

```bash
pip install -r doc/docs/requirements.txt
```

### Build and Serve Locally

```bash
# Serve with live reload (for development)
cd doc && mkdocs serve

# Build static HTML
cd doc && mkdocs build
```

The documentation will be available at `http://127.0.0.1:8000/`

### Build for Deployment

```bash
mkdocs build
```

This creates a `site/` directory with static HTML files.

## Documentation Structure

- `index.md` - Homepage
- `getting-started/` - Installation and quick start guides
- `api/` - API reference documentation
- `tutorials/` - In-depth tutorials
- `examples.md` - Complete working examples
- `about.md` - About PyGEL3D and GEL

## Adding New Pages

1. Create a new `.md` file in the appropriate directory
2. Add it to the `nav` section in `doc/mkdocs.yml`
3. Test with `cd doc && mkdocs serve`

## Theme

The documentation uses the Material for MkDocs theme with:
- Dark/light mode toggle
- Search functionality
- Code syntax highlighting
- Navigation tabs and sections

## Auto-generated API Docs

API documentation is automatically generated from Python docstrings using mkdocstrings.
To update:
1. Edit docstrings in the Python source files
2. Rebuild the documentation

## Deployment

To deploy to GitHub Pages:

```bash
cd doc && mkdocs gh-deploy
```

This builds the docs and pushes them to the `gh-pages` branch.
