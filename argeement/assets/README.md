# Static figures of the web page

`fig_workflow.svg` (processing workflow) and `fig_hierarchy.svg` (part of the NCBI gene hierarchy)
do not depend on the release; they are drawn with TikZ (the same code as Figs 1 and 2 of the
paper) and converted with `pdftocairo -svg`. The other figures of the page are made for every
release by `manuscript/figures.R` and come with the package (`site/img/`).
