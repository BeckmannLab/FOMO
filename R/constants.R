# requireNamespace("igraph", quietly = TRUE)
# requireNamespace("gtools", quietly = TRUE)
EMPTY_RELABELS <- data.frame(
    relabel_from = character(0),
    relabel_to = character(0)
)
VISNETWORK_LABEL_DOMAIN_SHAPES <- c(
    "dot",
    "square",
    "triangle",
    "diamond",
    "star"
)
LABEL_NOT_FOUND <- "LABELNOTFOUND"

## Colors used by the package's graph-plotting helpers (.generate_graph() in
## helpers-solve.R and .generate_corrections_graph() in helpers.R), pulled
## out as named constants (previously repeated string literals in both
## places) so the plotting palette can be changed in one place. Addresses a
## long-standing TODO ("Make the colors used in plotting global variables")
## at the top of MislabelSolver.R.
PLOT_COLOR_REGULAR_SAMPLE <- "orange"
PLOT_COLOR_ANCHOR_SAMPLE <- "forestgreen"
PLOT_COLOR_LABEL_NOT_FOUND <- "firebrick"
PLOT_COLOR_GHOST <- "lightgrey"
PLOT_COLOR_CONCORDANT_EDGE <- "forestgreen"
PLOT_COLOR_GENOTYPE_EDGE <- "orange"
PLOT_COLOR_NON_GENOTYPE_EDGE <- "cornflowerblue"
PLOT_COLOR_DEFAULT_EDGE <- "black"
