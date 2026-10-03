# Manuscript

`manuscript.md` is the source of the manuscript. It renders as an Elsevier article with the `elsarticle` class, following
Elsevier's LaTeX instructions: author-year citations with `elsarticle-harv`, and a flat submission folder, since
Elsevier's submission system cannot process LaTeX files in subfolders.

| File | What it holds |
|---|---|
| `manuscript.md` | the text, with citations as pandoc keys (`[@cunge1969]`, `@nash1959`) |
| `references.bib` | the bibliography |
| `metadata.yaml` | front matter (authors, affiliations, keywords, highlights) and class options |
| `latex/elsarticle.latex` | the pandoc template for `elsarticle` |
| `latex/manuscript.lua` | the pandoc filter that maps the Markdown's conventions onto LaTeX |
| `Makefile` | the build |

## Build

Needs pandoc and a TeX distribution with `xelatex`, `latexmk`, and `elsarticle` (MacTeX has all three):

```bash
brew install --cask mactex-no-gui   # once; asks for an administrator password
make tex          # build/manuscript.tex with the bibliography and figures beside it (pandoc only)
make pdf          # build/manuscript.pdf
make submission   # manuscript-submission.zip, every file at one folder level
```

## Writing conventions the filter relies on

- The level-one heading is the title, the `## Abstract` section the abstract; `## References` is ignored, since the
  bibliography comes from `references.bib`.
- Section numbers typed in headings (`## 2 Theory`) are dropped; LaTeX numbers sections, figures, and tables in order of
  appearance, so the numbers typed in the text must follow that order.
- A figure is an image followed by an italic paragraph that starts `*Figure N.`; a table is preceded by an italic
  paragraph that starts `*Table N.`. Figures are taken as the PDF beside each PNG.
- Symbols typed as text (c₁, S₀, x_N, 𝐈) are set as math; anything longer belongs in `$...$`.
- Class options in `metadata.yaml`: `preprint` for submission, or `1p`, `3p`, `5p` for a typeset preview, and `review`
  for double spacing. Line numbers are on by default (`linenumbers: true`).
