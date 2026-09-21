# Makefile — QI-Bench-Pareto-v1
# External review round 7 (author-approved full scope): a single-command
# verification pipeline, requested across multiple rounds of review.
#
# Assumes the `qibench` conda environment (pymoo, numpy, scipy, networkx,
# matplotlib, cartopy, pytest) is already active, or PYTHON is overridden:
#   make test PYTHON=/path/to/python

PYTHON ?= python3
CODE := code
PAPER := ../paper2

.PHONY: test figures tables paper verify package clean-pycache

test:
	cd $(CODE)/.. && $(PYTHON) -m pytest tests/ -v

figures:
	cd $(CODE) && $(PYTHON) generate_figures_v3.py
	cd $(CODE) && $(PYTHON) generate_topology_maps_v3.py
	cd $(CODE) && $(PYTHON) generate_convergence_figure_v3.py

tables:
	cd $(CODE) && $(PYTHON) generate_tables_v3.py

paper:
	cd $(PAPER) && pdflatex -interaction=nonstopmode -halt-on-error main.tex
	cd $(PAPER) && bibtex main
	cd $(PAPER) && pdflatex -interaction=nonstopmode -halt-on-error main.tex
	cd $(PAPER) && pdflatex -interaction=nonstopmode -halt-on-error main.tex

clean-pycache:
	find $(CODE) -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true

checksums: clean-pycache
	$(PYTHON) ../code/make_checksums.py

# `make verify` is the full gate: tests, then a checksum check against
# whatever CHECKSUMS.sha256 currently says (does NOT regenerate it --
# run `make checksums` first if you've changed files and want to
# re-baseline, per the project's established "checksums last" convention).
verify: test
	sha256sum -c CHECKSUMS.sha256
	@echo "All verification checks passed."

package: clean-pycache checksums
	@echo "Rebuild the REVIEW_PACKAGE/ + zip via the packaging steps in README.md / CHANGELOG.md."
	@echo "(Kept as a manual step, not automated here, since it also copies from ../paper2/"
	@echo " and the exact zip name is dated per round -- see CHANGELOG.md for the convention.)"
