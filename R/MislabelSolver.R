## TODOs
## If given a genotype_matrix, plot the graph derived directly from the matrix rather than the factor
## In global search, handle the case where there are more genotype groups than subjects
## Collapse cycles in the relabel graph if a single cycle involves 2 samples from the same genotype group (we specifically avoid this case in majority and global search, but this can still occur from a swap applied in local search)
## Nice features to have::
## Add a feature where researcher provides samples known to be mislabeled ahead of time, so relabeling them isn't penalized
## Add a feature where researcher groups samples together that can be mislabeled together (for example samples from the same blood draw), so they aren't penalized as separate sample relabels
## (Done -- see R/constants.R's PLOT_COLOR_* constants: "Make the colors used in plotting global variables")
## (Done -- see Solved_By in .solve_state$relabel_data / writeOutput()'s Sample and Component sheets: "Output which part of solver was able to solve each component")

#' The MislabelSolver class
#'
#' The MislabelSolver class stores the sample metadata required for mislabel detection and correction.
#'
#' @slot sample_metadata A data.frame containing sample metadata, with one row per sample.
#'                       Must include columns for Sample_ID, Subject_ID. If 'genotype_matrix'
#'                       is not provided, must also include a Genotype_Group_ID column. Any
#'                       Sample_ID(s) with NA in the Genotype_Group_ID column will be
#'                       tagged as phantom samples.
#' @slot genotype_matrix (Optional) A numeric or logical matrix specifiying whether a pair of
#'                       samples came from the same person. Row and column names must come from
#'                       Sample_ID column in 'sample_metadata'. Must be square, must be symmetric.
#'                       Any Sample_ID(s) that don't have row/columns in this matrix will be
#'                       tagged as phantom samples.
#' @slot label_domains (Optional) A data.frame with one row per sample specifying the Label_Domain,
#'                 where by experimental design only samples with the same Label_Domain may be
#'                 swapped for one another. For example, assay type or batch ID information
#'                 can be used to categorize Sample_ID(s) into Label_Domain(s)
#' @slot anchor_samples (Optional) A character vector of Sample_ID(s) where the label is known to be correct
#' @slot .solve_state A purely internal slot, used to keep track of sample relabels
#'
#' @return NULL
#'
#' @export
#'
setClass(
    "MislabelSolver",
    representation(
        sample_metadata = "data.frame",
        genotype_matrix = "ANY",
        label_domains = "data.frame",
        anchor_samples = "character",
        .solve_state = "list"
    ),
    prototype(
        genotype_matrix = NULL,
        label_domains = NULL,
        anchor_samples = character(0)
    )
)


#' Constructor for the MislabelSolver class
#'
#' @param sample_metadata A data.frame containing sample metadata, with one row per sample.
#'                        Must include columns for Sample_ID, Subject_ID. If 'genotype_matrix'
#'                        is not provided, must also include a Genotype_Group_ID column.
#' @param genotype_matrix (Optional) A numeric or logical matrix specifiying whether a pair of
#'                       samples came from the same person. Row and column names must come from
#'                       Sample_ID column in 'sample_metadata'. Must be square, must be symmetric
#' @param label_domains (Optional) A data.frame with one row per sample specifying the Label_Domain,
#'                  where by experimental design only samples with the same Label_Domain may be
#'                  swapped for one another. For example, assay type or batch ID information
#'                  can be used to categorize Sample_ID(s) into Label_Domain(s)
#' @param anchor_samples (Optional) A character vector of Sample_ID(s) where the label is known to be correct
#'
#' @details
#' \code{sample_metadata} and \code{label_domains} are sorted by \code{Sample_ID}
#' internally, so the constructed object (and everything later solved from
#' it) does not depend on the row order they happen to be provided in. One
#' consequence of this: when a sample's mislabel can only be resolved by
#' treating it as a duplicate of a sample that isn't actually present in the
#' data, the placeholder ID generated for it is deterministic given the same
#' input -- constructing a MislabelSolver on the same data twice always
#' generates the same placeholder ID(s), while constructing one on different
#' data generates different ones (rather than, say, either changing on every
#' run regardless of input, or always being the same regardless of input).
#'
#' @return A MislabelSolver object
#'
#' @export
#'
MislabelSolver <- function(
    sample_metadata,
    genotype_matrix = NULL,
    label_domains = NULL,
    anchor_samples = character(0)
) {
    ## Convert and validate inputs
    sample_metadata <- as.data.frame(lapply(sample_metadata, as.character))
    .validate_sample_metadata(
        sample_metadata,
        has_genotype_matrix = !is.null(genotype_matrix)
    )

    if (!is.null(genotype_matrix)) {
        .validate_genotype_matrix(genotype_matrix, sample_metadata)
        genotype_df <- .genotype_matrix_to_genotype_df(genotype_matrix)
        sample_metadata <- sample_metadata |>
            left_join(genotype_df, by = "Sample_ID")
    }

    if (is.null(label_domains)) {
        label_domains <- sample_metadata[, "Sample_ID", drop = FALSE]
        label_domains$Label_Domain <- "Label_Domain1"
    }
    label_domains <- as.data.frame(lapply(label_domains, as.character))
    .validate_label_domains(sample_metadata, label_domains)

    anchor_samples <- unique(as.character(anchor_samples))
    .validate_anchor_samples(sample_metadata, anchor_samples)

    return(new(
        "MislabelSolver",
        sample_metadata,
        genotype_matrix,
        label_domains,
        anchor_samples
    ))
}

setMethod(
    "initialize",
    "MislabelSolver",
    function(
        .Object,
        sample_metadata,
        genotype_matrix = NULL,
        label_domains = NULL,
        anchor_samples = character(0)
    ) {
        assert_that(
            is.data.frame(sample_metadata),
            all(c("Sample_ID", "Subject_ID") %in% colnames(sample_metadata))
        )
        ## Sort deterministically by Sample_ID, so that construction (and
        ## everything derived from it -- relabel_data's row order, and the
        ## placeholder IDs generated just below) does not depend on the row
        ## order 'sample_metadata'/'label_domains' happened to be provided in.
        sample_metadata <- sample_metadata[
            order(sample_metadata$Sample_ID),
            ,
            drop = FALSE
        ]
        rownames(sample_metadata) <- NULL
        if (!is.null(label_domains)) {
            assert_that(
                is.data.frame(label_domains),
                all(c("Sample_ID", "Label_Domain") %in% colnames(label_domains))
            )
            label_domains <- label_domains[
                order(label_domains$Sample_ID),
                ,
                drop = FALSE
            ]
            rownames(label_domains) <- NULL
        }
        ## Pre-generate a random placeholder Sample_ID for every sample, once,
        ## up front, for .find_relabel_cycles_from_putative_subjects() to use
        ## later if (and only if) that specific sample is determined to need one
        ## -- i.e. no real or ghost sample is available to relabel it to, so its
        ## mislabel can only be resolved by treating it as a duplicate of a
        ## sample that doesn't actually exist ("unknown"/LABELNOTFOUND labels;
        ## see helpers-solve.R). A sample whose pre-generated ID never ends up
        ## needed is simply never referenced again.
        ##
        ## Seeded from a hash of this (now sorted) input via a dedicated RNG
        ## stream (with_seed(), which restores the prior RNG state afterward),
        ## so that: (a) the same input always yields the same placeholder IDs,
        ## (b) different input yields different ones rather than colliding on
        ## one fixed global seed, and (c) this doesn't consume from or interfere
        ## with the RNG stream solveMajoritySearch()/
        ## solveGlobalSearch()/solveLocalSearch() use for their own purposes.
        input_seed <- .hash_to_seed(list(
            sample_metadata = sample_metadata,
            label_domains = label_domains
        ))
        placeholder_ids <- with_seed(
            input_seed,
            .generate_placeholder_ids(nrow(sample_metadata))
        )
        names(placeholder_ids) <- sample_metadata$Sample_ID

        ## Provided there are enough shapes, assign a unique shape to
        ## each Label_Domain
        all_label_domain_ids <- names(sort(
            table(label_domains$Label_Domain),
            decreasing = TRUE
        ))
        label_domain_shapes <- data.frame(
            Label_Domain = all_label_domain_ids,
            Label_Domain_Shape = "dot",
            vertex_size_scalar = 1
        )
        if (
            length(all_label_domain_ids) <=
                length(VISNETWORK_LABEL_DOMAIN_SHAPES)
        ) {
            label_domain_shapes$Label_Domain_Shape <- VISNETWORK_LABEL_DOMAIN_SHAPES[seq_along(
                all_label_domain_ids
            )]
        }
        label_domains <- label_domains |>
            left_join(label_domain_shapes, by = "Label_Domain")

        ## Initialize object 'solve_state'
        relabel_data <- sample_metadata |>
            mutate(
                Init_Sample_ID = .data$Sample_ID,
                Init_Subject_ID = .data$Subject_ID,
                Is_Ghost = is.na(.data$Genotype_Group_ID),
                Is_Anchor = .data$Init_Sample_ID %in% anchor_samples,
                Solved = FALSE,
                Relabeled_By = NA_character_,
                Placeholder_ID = placeholder_ids[.data$Sample_ID]
            ) |>
            left_join(label_domains, by = "Sample_ID")
        unsolved_relabel_data <- relabel_data |>
            filter(!is.na(.data$Genotype_Group_ID))
        unsolved_ghost_data <- relabel_data |>
            filter(is.na(.data$Genotype_Group_ID))
        putative_subjects <- data.frame(
            Genotype_Group_ID = character(0),
            Subject_ID = character(0)
        )
        lnf_counts <- data.frame(
            Subject_ID = character(0),
            Label_Domain = character(0),
            count = integer(0)
        )
        ambiguous_subjects <- list()
        solve_state <- list(
            relabel_data = relabel_data,
            unsolved_relabel_data = unsolved_relabel_data,
            unsolved_ghost_data = unsolved_ghost_data,
            putative_subjects = putative_subjects,
            lnf_counts,
            ambiguous_subjects = ambiguous_subjects
        )

        .Object@sample_metadata <- sample_metadata
        .Object@genotype_matrix <- genotype_matrix
        .Object@label_domains <- label_domains
        .Object@anchor_samples <- anchor_samples
        .Object@.solve_state <- solve_state

        .Object <- .update_solve_state(.Object, initialization = TRUE)
        return(.Object)
    }
)
