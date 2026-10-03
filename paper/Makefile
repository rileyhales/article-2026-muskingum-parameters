# Render manuscript.md as an Elsevier (elsarticle) manuscript.
#
#   make tex         build/manuscript.tex, with the bibliography and figures beside it (needs pandoc)
#   make pdf         build/manuscript.pdf (needs a TeX distribution with xelatex, latexmk, and elsarticle)
#   make submission  manuscript-submission.zip: every file Elsevier needs at one folder level, no subfolders
#   make clean       remove build/ and the zip
#
# build/ is flat on purpose: Elsevier's submission system cannot process LaTeX in subfolders.

SHELL := /bin/bash
export PATH := /Library/TeX/texbin:$(PATH)

BUILD := build
SOURCES := manuscript.md metadata.yaml references.bib latex/elsarticle.latex latex/manuscript.lua

.PHONY: tex pdf submission clean

tex: $(BUILD)/manuscript.tex

$(BUILD)/manuscript.tex: $(SOURCES)
	mkdir -p $(BUILD)
	pandoc manuscript.md \
		--from markdown --to latex --standalone --natbib \
		--shift-heading-level-by=-1 \
		--template latex/elsarticle.latex \
		--lua-filter latex/manuscript.lua \
		--metadata-file metadata.yaml \
		--metadata flat-figures=true \
		--metadata figure-list=$(BUILD)/figures.txt \
		--output $@
	cp references.bib $(BUILD)/
	while read -r figure; do cp "$$figure" $(BUILD)/; done < $(BUILD)/figures.txt

pdf: $(BUILD)/manuscript.pdf

$(BUILD)/manuscript.pdf: $(BUILD)/manuscript.tex
	@command -v latexmk >/dev/null || { echo 'latexmk not found: install MacTeX (brew install --cask mactex-no-gui)'; exit 1; }
	cp "$$(kpsewhich elsarticle-harv.bst)" $(BUILD)/
	cd $(BUILD) && latexmk -xelatex -interaction=nonstopmode -halt-on-error manuscript.tex

submission: $(BUILD)/manuscript.pdf
	rm -f manuscript-submission.zip
	cd $(BUILD) && zip -q ../manuscript-submission.zip manuscript.tex manuscript.bbl references.bib \
		elsarticle-harv.bst $$(cat figures.txt | xargs -n1 basename)

clean:
	rm -rf $(BUILD) manuscript-submission.zip
