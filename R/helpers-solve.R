## Shared validation for the `ghost_penalty`/`deletion_penalty` arguments used
## by solveGlobalSearch(). The relabel
## penalty is fixed at 1 (all three penalties are relative to each other, so
## fixing one just sets the scale); `ghost_penalty` and `deletion_penalty`
## default to the package's current weights (1.5 and 4) but can be
## overridden. (The deletion default was doubled from 2 to 4 -- see
## solveGlobalSearch()'s documentation -- when a double-counting bug in the
## deletion-penalty scoring was fixed, to keep the effective per-sample
## deletion cost unchanged from earlier package versions in the common case.)
##
## Recommended limits, and why: global search scores a candidate permutation
## as 1 point per ordinary relabel, `ghost_penalty` points per relabel to a
## ghost (placeholder) sample, and `deletion_penalty` points per sample whose
## label or genotype it gives up on entirely. `ghost_penalty` must be
## strictly greater than the relabel penalty (1), or the algorithm may prefer
## inventing a ghost-sample relabel even when a real, non-ghost sample swap
## is already available.
##
## `deletion_penalty` needs to clear two *separate* bars, both derived from
## minimal worked examples (and confirmed against the real scoring code, not
## just argued abstractly -- see the git history of this comment for the
## analysis): it must be strictly greater than 2*relabel_penalty (not just
## "at least" -- at exactly 2 the two options below are an exact tie, decided
## only by incidental permutation-enumeration order, which is not a real
## guarantee of anything), AND strictly greater than `ghost_penalty`.
##
##   1. > 2*relabel_penalty: with two samples that simply need to swap
##      identities with each other (the textbook case -- see
##      toy_swap_scenario() in tests/testthat/helper-toy-scenarios.R, also
##      the scenario that first surfaced the double-counting bug above), the
##      honest fix costs 2 (one relabel each way) while treating the pair as
##      one deletion instead costs exactly `deletion_penalty`. At
##      deletion_penalty <= 2, giving up scores as good as or better than the
##      honest fix.
##   2. > ghost_penalty: separately, when a leftover sample could be
##      resolved either by relabeling it to an available ghost (cost
##      `ghost_penalty`) or by counting it as a deletion (cost
##      `deletion_penalty`), deletion_penalty <= ghost_penalty makes the
##      algorithm prefer giving up over using a ghost that was right there.
##      This can happen even when deletion_penalty already clears the "> 2"
##      bar above -- the two conditions are independent, and both are needed;
##      confirmed with a small hand-built scenario (one real conflicting
##      pair plus one available ghost) where deletion_penalty = 2.5 (> 2) and
##      ghost_penalty = 3 (> 1) still made the algorithm delete rather than
##      use the ghost, purely because 2.5 < 3.
##
## Neither bar is a hard guarantee of always preferring the honest fix in
## *every* possible component, for every combination of component size and
## shape -- stress-testing this scoring formula against randomly generated
## components found configurations (larger, more lopsided ones, beyond the
## two-samples-and-a-neighbor minimal case above) that need a considerably
## higher deletion_penalty than max(2, ghost_penalty) to still prefer the
## fewer-deletion reading, and no fixed constant works for every component
## shape and size. What's checked here is the necessary bar cleared by the
## simplest, most common cases (matching the spirit of the pre-existing
## ghost_penalty check, which has the same "warn, don't enforce" limitation).
## Non-positive values are always rejected, since a non-positive penalty
## would make the search prefer *more* ghost relabels or deletions, the
## opposite of a penalty's purpose.
.validate_search_penalties <- function(ghost_penalty, deletion_penalty) {
    relabel_penalty <- 1

    assert_that(
        is_bare_numeric(ghost_penalty, n = 1) && !is.na(ghost_penalty),
        msg = "'ghost_penalty' must be a single (length-1) numeric value"
    )
    assert_that(
        is_bare_numeric(deletion_penalty, n = 1) && !is.na(deletion_penalty),
        msg = "'deletion_penalty' must be a single (length-1) numeric value"
    )
    assert_that(
        ghost_penalty > 0,
        msg = "'ghost_penalty' must be positive"
    )
    assert_that(
        deletion_penalty > 0,
        msg = "'deletion_penalty' must be positive"
    )

    if (ghost_penalty <= relabel_penalty) {
        warning(
            "'ghost_penalty' (",
            ghost_penalty,
            ") is not greater than the relabel penalty (",
            relabel_penalty,
            "); the algorithm may prefer relabeling to a ghost sample ",
            "even when a valid non-ghost sample swap is already available.",
            call. = FALSE
        )
    }
    min_deletion_penalty <- max(2 * relabel_penalty, ghost_penalty)
    if (deletion_penalty <= min_deletion_penalty) {
        warning(
            "'deletion_penalty' (",
            deletion_penalty,
            ") is not greater than max(twice the ",
            "relabel penalty, ghost_penalty) = ",
            min_deletion_penalty,
            "; the algorithm may ",
            "prefer inserting/deleting samples even when a valid non-ghost sample swap or an ",
            "available ghost sample already exists.",
            call. = FALSE
        )
    }

    invisible(NULL)
}

## Returns the sorted, unique Sample_ID(s) (from both unsolved_relabel_data
## and unsolved_ghost_data) that currently belong to a component
## solveGlobalSearch() would actually attempt to process (not skip) if
## called right now, restricted to the specific real samples whose
## Genotype_Group_ID is not already locked via putative_subjects (i.e. is
## genuinely still "free" for solveGlobalSearch() to decide), plus any
## ghost samples in that same component (ghosts have no Genotype_Group_ID
## of their own to lock, but can still be newly relevant to a component
## that hasn't changed status otherwise). Mirrors the per-component skip
## conditions in solveGlobalSearch()'s own loop: a component is excluded
## here if it has more Genotype_Group_ID(s) than Subject_ID(s) (an
## unrelated, locking-independent reason solveGlobalSearch() always skips
## it), or if its *free* genotype/subject counts exceed max_genotypes.
##
## Used by solveEnsemble() to detect when calling solveGlobalSearch() again
## would be futile -- see the comment above where it's called, in
## solveEnsemble() itself, for the full reasoning. Deliberately keyed on
## free (not-yet-locked) genotypes rather than plain component membership:
## solveMajoritySearch() pays no attention to component boundaries or
## max_genotypes and can lock a component's individual genotypes without
## changing which samples belong to it at all, and that newly-available
## information is exactly what a membership-only comparison would miss.
##
## This does duplicate part of solveGlobalSearch()'s own per-component
## logic (deliberately: it's cheap, and keeping solveGlobalSearch() itself
## completely unmodified was judged the lower-risk option here), so if
## solveGlobalSearch()'s eligibility conditions ever change, this needs to
## be updated to match.
.global_search_available_samples <- function(object, max_genotypes) {
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
    unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data
    putative_subjects <- object@.solve_state$putative_subjects

    if (nrow(unsolved_relabel_data) == 0) {
        return(character(0))
    }

    component_ids <- sort(unique(unsolved_relabel_data$Component_ID))
    available_sample_ids <- character(0)
    for (component_id in component_ids) {
        cc_unsolved_relabel_data <- unsolved_relabel_data |>
            filter(.data$Component_ID == component_id)
        cc_unsolved_ghost_data <- unsolved_ghost_data |>
            filter(.data$Component_ID == component_id)
        cc_genotypes <- unique(cc_unsolved_relabel_data$Genotype_Group_ID)
        cc_subjects <- unique(cc_unsolved_relabel_data$Subject_ID)

        ## Same locking-independent structural skip solveGlobalSearch() uses
        if (length(cc_genotypes) > length(cc_subjects)) {
            next
        }

        free_genotypes <- setdiff(
            cc_genotypes,
            putative_subjects$Genotype_Group_ID
        )
        free_subjects <- setdiff(cc_subjects, putative_subjects$Subject_ID)

        ## Same size-based skip solveGlobalSearch() uses, but based on the
        ## *free* (locking-adjusted) counts, not the raw component size
        if (
            length(free_genotypes) > max_genotypes ||
                length(free_subjects) > max_genotypes
        ) {
            next
        }

        available_sample_ids <- c(
            available_sample_ids,
            cc_unsolved_relabel_data$Sample_ID[
                cc_unsolved_relabel_data$Genotype_Group_ID %in%
                    free_genotypes
            ],
            cc_unsolved_ghost_data$Sample_ID
        )
    }
    sort(unique(available_sample_ids))
}

## Component_ID(s) that solveMajoritySearch()'s max_genotypes argument says
## to leave alone this round: components where the number of not-yet-locked
## (i.e. not already present in putative_subjects) Genotype_Group_ID(s) or
## Subject_ID(s) exceeds max_genotypes. This mirrors the free-genotype/
## free-subject size check solveGlobalSearch() uses for its own
## max_genotypes (see solveGlobalSearch() and
## .global_search_available_samples() above), but deliberately omits the
## additional "more genotypes than subjects" structural skip those use --
## that skip is specific to global search's permutation-based approach and
## has nothing to do with how expensive a component is to vote/cycle-search
## over.
.majority_search_skip_components <- function(
    unsolved_relabel_data,
    putative_subjects,
    max_genotypes
) {
    component_ids <- sort(unique(unsolved_relabel_data$Component_ID))
    skip_component_ids <- character(0)
    for (component_id in component_ids) {
        cc_unsolved_relabel_data <- unsolved_relabel_data |>
            filter(.data$Component_ID == component_id)
        cc_genotypes <- unique(cc_unsolved_relabel_data$Genotype_Group_ID)
        cc_subjects <- unique(cc_unsolved_relabel_data$Subject_ID)

        free_genotypes <- setdiff(
            cc_genotypes,
            putative_subjects$Genotype_Group_ID
        )
        free_subjects <- setdiff(cc_subjects, putative_subjects$Subject_ID)

        if (
            length(free_genotypes) > max_genotypes ||
                length(free_subjects) > max_genotypes
        ) {
            skip_component_ids <- c(skip_component_ids, component_id)
        }
    }
    skip_component_ids
}

## Mirrors .global_search_available_samples() above, but for
## solveMajoritySearch()'s own max_genotypes: the set of Sample_ID(s)
## belonging to a not-yet-locked (free) Genotype_Group_ID within a component
## small enough for majority search to actually consider this round. Used by
## solveEnsemble() to skip a redundant solveMajoritySearch() call when
## nothing new has become available to it since it last ran -- see
## .global_search_available_samples()'s own comment for why that's safe to
## do. Deliberately does not apply the "more genotypes than subjects"
## structural skip .global_search_available_samples() uses (see
## .majority_search_skip_components() above for why), and does not include
## ghost samples, since solveMajoritySearch() never relabels them.
.majority_search_available_samples <- function(object, max_genotypes) {
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
    putative_subjects <- object@.solve_state$putative_subjects

    if (nrow(unsolved_relabel_data) == 0) {
        return(character(0))
    }

    component_ids <- sort(unique(unsolved_relabel_data$Component_ID))
    available_sample_ids <- character(0)
    for (component_id in component_ids) {
        cc_unsolved_relabel_data <- unsolved_relabel_data |>
            filter(.data$Component_ID == component_id)
        cc_genotypes <- unique(cc_unsolved_relabel_data$Genotype_Group_ID)
        cc_subjects <- unique(cc_unsolved_relabel_data$Subject_ID)

        free_genotypes <- setdiff(
            cc_genotypes,
            putative_subjects$Genotype_Group_ID
        )
        free_subjects <- setdiff(cc_subjects, putative_subjects$Subject_ID)

        if (
            length(free_genotypes) > max_genotypes ||
                length(free_subjects) > max_genotypes
        ) {
            next
        }

        available_sample_ids <- c(
            available_sample_ids,
            cc_unsolved_relabel_data$Sample_ID[
                cc_unsolved_relabel_data$Genotype_Group_ID %in%
                    free_genotypes
            ]
        )
    }
    sort(unique(available_sample_ids))
}

## solver_name attributes newly-solved samples to whoever solved them, for
## the "Solved_By" column surfaced in writeOutput()'s Sample/Component
## sheets (see solveEnsemble()'s call sites for the actual solver names
## used: "majority", "global", "local"/"local_old"). Not
## meaningful when initialization=TRUE (no solver has run yet); any sample
## that comes out solved at construction time -- e.g. because its whole
## component was already internally consistent, or anchor_samples pinned it
## -- is attributed to "initial" instead, regardless of what (if anything)
## solver_name was passed.
.update_solve_state <- function(
    object,
    initialization = FALSE,
    solver_name = NA_character_
) {
    if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
        return(object)
    }
    effective_solver_name <- if (initialization) {
        "initial"
    } else {
        solver_name
    }

    ## 1. Assign a Component_ID for each unsolved Sample_ID
    combined_graph <- .generate_graph(
        object@.solve_state$unsolved_relabel_data,
        graph_type = "combined",
        object@.solve_state$unsolved_ghost_data
    )
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data |>
        mutate(
            Component_ID = as.character(components(combined_graph)$membership[
                .data$Sample_ID
            ])
        )
    unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data |>
        mutate(
            Component_ID = as.character(components(combined_graph)$membership[
                .data$Sample_ID
            ])
        )

    ## 2. Determine which Sample_ID(s) are newly solved
    ##    A Sample_ID is solved if it belongs to a Component_ID that includes only one
    ##    Genotype_Group_ID and one Subject_ID
    component_data <- rbind(unsolved_relabel_data, unsolved_ghost_data) |>
        group_by(.data$Component_ID) |>
        summarise(
            n_Genotype_Group_ID = n_distinct(.data$Genotype_Group_ID) -
                anyNA(.data$Genotype_Group_ID),
            n_Subject_ID = n_distinct(.data$Subject_ID),
            n_Sample_ID = length(.data$Sample_ID)
        ) |>
        mutate(
            Solved = .data$n_Genotype_Group_ID <= 1 & .data$n_Subject_ID == 1
        )

    solved_components <- component_data |>
        filter(.data$Solved) |>
        pull(.data$Component_ID)

    unsolved_relabel_data <- unsolved_relabel_data |>
        mutate(
            Solved = .data$Component_ID %in% solved_components,
            Solved_By = if_else(
                .data$Solved,
                effective_solver_name,
                NA_character_
            )
        )
    unsolved_ghost_data <- unsolved_ghost_data |>
        mutate(
            Solved = .data$Component_ID %in% solved_components,
            Solved_By = if_else(
                .data$Solved,
                effective_solver_name,
                NA_character_
            )
        )

    ## 3. Re-rank Component_ID(s) in order of size (so that Component1 is the largest unsolved component)
    # component_data <- component_data |> filter(!Solved)
    n_components <- nrow(component_data)
    component_data <- component_data |>
        arrange(.data$Solved, desc(.data$n_Sample_ID)) |>
        mutate(
            new_Component_ID = if (n_components != 0) {
                seq_len(n_components)
            } else {
                character(0)
            },
            new_Component_ID = str_c(
                "Component_",
                formatC(
                    .data$new_Component_ID,
                    width = str_length(n_components),
                    format = "d",
                    flag = "0"
                )
            )
        )
    unsolved_relabel_data <- unsolved_relabel_data |>
        left_join(
            component_data[, c("Component_ID", "new_Component_ID")],
            by = "Component_ID"
        ) |>
        mutate(Component_ID = .data$new_Component_ID) |>
        select(-"new_Component_ID")
    unsolved_ghost_data <- unsolved_ghost_data |>
        left_join(
            component_data[, c("Component_ID", "new_Component_ID")],
            by = "Component_ID"
        ) |>
        mutate(Component_ID = .data$new_Component_ID) |>
        select(-"new_Component_ID")

    ## 4. Update putative_subjects
    ##    During initialization, lock Subject_ID/Genotype_Group_ID pairs for anchor_samples
    if (initialization) {
        anchor_putative_subjects <- object@sample_metadata |>
            filter(
                !is.na(.data$Genotype_Group_ID),
                !is.na(.data$Subject_ID),
                .data$Sample_ID %in% object@anchor_samples
            ) |>
            select("Genotype_Group_ID", "Subject_ID") |>
            distinct()
        object <- .update_putative_subjects(object, anchor_putative_subjects)
    }
    solved_putative_subjects <- unsolved_relabel_data |>
        filter(.data$Solved) |>
        select("Genotype_Group_ID", "Subject_ID") |>
        distinct()
    object <- .update_putative_subjects(object, solved_putative_subjects)

    ## 5. Update relabel_data, and unsolved_relabel_data
    if (initialization) {
        col_order <- c(
            "Init_Sample_ID",
            "Init_Subject_ID",
            "Genotype_Group_ID",
            "Component_ID",
            "Sample_ID",
            "Subject_ID",
            "Solved",
            "Solved_By",
            "Relabeled_By",
            "Is_Ghost",
            "Is_Anchor",
            "Label_Domain",
            "Label_Domain_Shape",
            "vertex_size_scalar",
            "Placeholder_ID"
        )
        unsolved_relabel_data <- unsolved_relabel_data |>
            select(all_of(col_order)) |>
            mutate(Init_Component_ID = .data$Component_ID) |>
            relocate("Init_Component_ID", .before = "Component_ID")
        unsolved_ghost_data <- unsolved_ghost_data |>
            select(all_of(col_order)) |>
            mutate(Init_Component_ID = .data$Component_ID) |>
            relocate("Init_Component_ID", .before = "Component_ID")
        relabel_data <- rbind(unsolved_relabel_data, unsolved_ghost_data)
    } else {
        ## Update 'relabel_data' with new sample labels in 'unsolved_relabel_data' and 'unsolved_ghost_data'
        unsolved_data <- rbind(unsolved_relabel_data, unsolved_ghost_data)
        relabel_data <- rows_update(
            object@.solve_state$relabel_data,
            unsolved_data,
            by = "Init_Sample_ID"
        )
    }
    unsolved_relabel_data <- unsolved_relabel_data |> filter(!.data$Solved)
    unsolved_ghost_data <- unsolved_ghost_data |> filter(!.data$Solved)

    ## 6. Overwrite .solve_state
    object@.solve_state$relabel_data <- relabel_data
    object@.solve_state$unsolved_relabel_data <- unsolved_relabel_data
    object@.solve_state$unsolved_ghost_data <- unsolved_ghost_data

    return(object)
}

.update_putative_subjects <- function(object, proposed_putative_subjects) {
    if (nrow(proposed_putative_subjects) == 0) {
        return(object)
    }
    ## Only add Genotype_Group_ID/Subject_ID combinations if neither the
    ## Genotype_Group_ID nor the Subject_ID are already in putative_subjects
    existing_genotypes <- na.omit(
        object@.solve_state$putative_subjects$Genotype_Group_ID
    )
    existing_subjects <- na.omit(
        object@.solve_state$putative_subjects$Subject_ID
    )
    proposed_putative_subjects <- proposed_putative_subjects |>
        filter(
            !(.data$Subject_ID %in% existing_subjects),
            !(.data$Genotype_Group_ID %in% existing_genotypes)
        )
    putative_subjects <- rbind(
        object@.solve_state$putative_subjects,
        proposed_putative_subjects
    )
    # .validate_putative_subjects(object@sample_genotype_data, putative_subjects)
    object@.solve_state$putative_subjects <- putative_subjects
    return(object)
}

.generate_graph <- function(
    relabel_data,
    graph_type = c("label", "genotype", "combined"),
    ghost_data = NULL,
    genotype_matrix = NULL,
    populate_plotting_attributes = FALSE,
    collapse_samples = FALSE
) {
    graph_type_mapping <- list(
        label = "Subject_ID",
        genotype = "Genotype_Group_ID",
        combined = NA_character_
    )

    ## Collapse samples mode should only run if the graph is being plotted
    if (!populate_plotting_attributes) {
        collapse_samples <- FALSE
    }

    graph_type <- as.character(graph_type)
    graph_type <- match.arg(graph_type)

    ## Ignore ghost samples when constructing genotype graph
    if (graph_type == "genotype") {
        ghost_data <- NULL
    }

    ## 'collapse_samples' renames grouped samples to synthetic vertex names
    ## (e.g. "3 samples\nA\nRNA"), which cannot be matched back to any
    ## rowname in 'genotype_matrix'. Rather than silently building a
    ## genotype-half graph with mismatched or missing vertices, fall back to
    ## not collapsing when a genotype_matrix is in play.
    if (collapse_samples && !is.null(genotype_matrix)) {
        warning(
            "'collapse_samples' is not supported for a MislabelSolver built ",
            "with a 'genotype_matrix' (as opposed to a 'Genotype_Group_ID' ",
            "column); plotting without collapsing samples."
        )
        collapse_samples <- FALSE
    }

    if (collapse_samples) {
        relabel_data <- relabel_data |>
            group_by(
                .data$Subject_ID,
                .data$Genotype_Group_ID,
                .data$Label_Domain
            ) |>
            mutate(
                count = n(),
                vertex_size_scalar = sqrt(sum(.data$vertex_size_scalar)),
                Is_Ghost = FALSE,
                Is_Anchor = any(.data$Is_Anchor),
                ## Include Genotype_Group_ID in the synthetic label so that two
                ## different groups that happen to share Subject_ID, Label_Domain,
                ## and collapsed size don't collapse to the same vertex name
                ## (which igraph rejects with "Duplicate vertex names").
                Sample_ID = if_else(
                    .data$count == 1,
                    .data$Sample_ID,
                    paste(
                        paste(.data$count, "samples"),
                        .data$Subject_ID,
                        .data$Genotype_Group_ID,
                        .data$Label_Domain,
                        sep = "\n"
                    )
                )
            ) |>
            select(
                "Sample_ID",
                "Subject_ID",
                "Genotype_Group_ID",
                "Label_Domain",
                "Label_Domain_Shape",
                "count",
                "vertex_size_scalar",
                "Is_Ghost",
                "Is_Anchor"
            ) |>
            distinct()
    }

    all_data <- relabel_data

    if (!is.null(ghost_data)) {
        if (collapse_samples) {
            ghost_data <- ghost_data |>
                group_by(
                    .data$Subject_ID,
                    .data$Genotype_Group_ID,
                    .data$Label_Domain
                ) |>
                mutate(
                    count = n(),
                    vertex_size_scalar = sqrt(sum(.data$vertex_size_scalar)),
                    Is_Ghost = TRUE,
                    Is_Anchor = FALSE,
                    Sample_ID = if_else(
                        .data$count == 1,
                        .data$Sample_ID,
                        paste(
                            paste(.data$count, "samples"),
                            .data$Subject_ID,
                            .data$Label_Domain,
                            sep = "\n"
                        )
                    )
                ) |>
                select(
                    "Sample_ID",
                    "Subject_ID",
                    "Genotype_Group_ID",
                    "Label_Domain",
                    "Label_Domain_Shape",
                    "count",
                    "vertex_size_scalar",
                    "Is_Ghost",
                    "Is_Anchor"
                ) |>
                distinct()
        }
        all_data <- rbind(relabel_data, ghost_data)
    }

    if (graph_type == "combined") {
        genotype_graph <- .generate_graph(
            relabel_data,
            "genotype",
            genotype_matrix = genotype_matrix
        )
        E(genotype_graph)$genotypes <- TRUE
        label_graph <- .generate_graph(relabel_data, "label", ghost_data)
        E(label_graph)$labels <- TRUE
        graph <- igraph::union(genotype_graph, label_graph, byname = TRUE)
        E(graph)[is.na(E(graph)$genotypes)]$genotypes <- FALSE
        E(graph)[is.na(E(graph)$labels)]$labels <- FALSE
        E(graph)$concordant <- E(graph)$genotypes & E(graph)$labels
    } else {
        if (graph_type == "genotype" & !is.null(genotype_matrix)) {
            vertices <- all_data[, "Sample_ID", drop = FALSE]
            ## Restrict the genotype matrix to the samples actually present in
            ## 'all_data' (e.g. the caller's 'unsolved'-only subset). Without
            ## this, the genotype half of the graph always contains every
            ## sample in the object, while the label half only contains the
            ## requested subset, so the two halves end up with mismatched
            ## vertex sets when unioned below.
            genotype_matrix <- genotype_matrix[
                rownames(genotype_matrix) %in% vertices$Sample_ID,
                colnames(genotype_matrix) %in% vertices$Sample_ID,
                drop = FALSE
            ]
            graph <- graph_from_adjacency_matrix(
                genotype_matrix,
                mode = "undirected"
            )
        } else {
            group_col <- graph_type_mapping[[graph_type]]
            edges <- all_data |>
                group_by_at(group_col) |>
                mutate(
                    sample_a = .data$Sample_ID,
                    sample_b = list(.data$Sample_ID)
                ) |>
                ungroup() |>
                unnest("sample_b") |>
                transmute(
                    sample1 = pmin(.data$sample_a, .data$sample_b),
                    sample2 = pmax(.data$sample_a, .data$sample_b)
                ) |>
                filter(.data$sample1 != .data$sample2) |>
                distinct()
            vertices <- all_data[, "Sample_ID", drop = FALSE]
            graph <- graph_from_data_frame(
                edges,
                vertices = vertices,
                directed = FALSE
            )
        }
    }

    if (!populate_plotting_attributes) {
        return(graph)
    }

    ## Specify vertex and edge attributes for plotting
    vertex_shapes <- data.frame(Sample_ID = names(V(graph))) |>
        left_join(all_data, by = "Sample_ID") |>
        pull(.data$Label_Domain_Shape)
    vertex_size_scalars <- data.frame(Sample_ID = names(V(graph))) |>
        left_join(all_data, by = "Sample_ID") |>
        pull(.data$vertex_size_scalar)
    V(graph)$shape <- vertex_shapes
    V(graph)$size <- 12 * vertex_size_scalars
    V(graph)$label.cex <- 0.5

    anchor_samples <- relabel_data |>
        filter(.data$Is_Anchor) |>
        pull(.data$Sample_ID)
    anchor_samples <- intersect(anchor_samples, V(graph)$name)
    label_not_found_samples <- relabel_data |>
        filter(str_detect(.data$Sample_ID, LABEL_NOT_FOUND)) |>
        pull(.data$Sample_ID)
    label_not_found_samples <- intersect(label_not_found_samples, V(graph)$name)
    ghost_samples <- NULL
    if (!is.null(ghost_data)) {
        ghost_samples <- ghost_data |> pull(.data$Sample_ID)
        ghost_samples <- intersect(ghost_samples, V(graph)$name)
    }
    V(graph)$color <- PLOT_COLOR_REGULAR_SAMPLE
    V(graph)[anchor_samples]$color <- PLOT_COLOR_ANCHOR_SAMPLE
    V(graph)[label_not_found_samples]$color <- PLOT_COLOR_LABEL_NOT_FOUND
    V(graph)[ghost_samples]$color <- PLOT_COLOR_GHOST

    if (graph_type == "combined") {
        E(graph)$color <- if_else(
            E(graph)$concordant,
            PLOT_COLOR_CONCORDANT_EDGE,
            if_else(
                E(graph)$genotypes,
                PLOT_COLOR_GENOTYPE_EDGE,
                PLOT_COLOR_NON_GENOTYPE_EDGE
            )
        )
        E(graph)[.from(ghost_samples)]$color <- PLOT_COLOR_GHOST
    } else if (graph_type == "label") {
        E(graph)$color <- PLOT_COLOR_NON_GENOTYPE_EDGE
        E(graph)[.from(ghost_samples)]$color <- PLOT_COLOR_GHOST
    } else {
        E(graph)$color <- PLOT_COLOR_GENOTYPE_EDGE
    }
    E(graph)$width <- 6

    return(graph)
}

.relabel_samples <- function(object, relabels, solver_name = NA_character_) {
    if (nrow(relabels) == 0) {
        return(object)
    }

    relabels <- relabels |>
        rename("Sample_ID" = "relabel_to") |>
        left_join(
            object@sample_metadata[, c("Sample_ID", "Subject_ID")],
            by = "Sample_ID"
        ) |>
        mutate(
            Subject_ID = if_else(
                is.na(.data$Subject_ID),
                vapply(
                    .data$Sample_ID,
                    \(x) str_split_1(x, "#")[2],
                    character(1)
                ),
                .data$Subject_ID
            ),
            Deleted_relabel_from = str_detect(
                .data$relabel_from,
                LABEL_NOT_FOUND
            )
        ) |>
        filter(!.data$Deleted_relabel_from) |>
        select(-"Deleted_relabel_from")

    ## Call it relabeled sample ID instead
    unsolved_all_data <- rbind(
        object@.solve_state$unsolved_relabel_data,
        object@.solve_state$unsolved_ghost_data
    )
    unsolved_all_data <- unsolved_all_data |>
        left_join(
            relabels,
            by = c("Sample_ID" = "relabel_from"),
            suffix = c(".x", ".y")
        ) |>
        mutate(
            ## coalesce(a, b): a if it's not NA, else b -- exactly what these
            ## two ifelse(!is.na(a), a, b) calls were doing.
            Sample_ID = coalesce(.data$Sample_ID.y, .data$Sample_ID),
            Subject_ID = coalesce(.data$Subject_ID.y, .data$Subject_ID.x),
            ## Unlike Solved_By (set once, permanently, when a component is
            ## solved), Relabeled_By is overwritten every time a sample is
            ## actually relabeled -- Sample_ID.y is only non-NA for rows that
            ## matched a real entry in 'relabels' (i.e. were just relabeled
            ## by this call), so untouched rows simply carry forward
            ## whatever value they already had (NA if never relabeled).
            Relabeled_By = if_else(
                !is.na(.data$Sample_ID.y),
                solver_name,
                .data$Relabeled_By
            )
        ) |>
        select(-ends_with(".x"), -ends_with(".y"))

    object@.solve_state$unsolved_relabel_data <- unsolved_all_data |>
        filter(!.data$Is_Ghost)
    object@.solve_state$unsolved_ghost_data <- unsolved_all_data |>
        filter(.data$Is_Ghost)
    object <- .update_solve_state(object, solver_name = solver_name)

    tsmsg(paste(nrow(relabels), "samples relabeled"))

    return(object)
}

# TODO: drop duplicate cycles
.find_directed_cycles <- function(graph, cutoff = 1) {
    assert_that(is_directed(graph), msg = "param 'graph' must be directed")
    cycles <- list()
    for (vertex in V(graph)) {
        in_neighbors <- names(neighbors(graph, vertex, mode = "in"))
        for (in_neighbor in in_neighbors) {
            simple_paths <- all_simple_paths(
                graph,
                vertex,
                in_neighbor,
                mode = "out",
                cutoff = cutoff
            )
            cycles <- append(cycles, simple_paths)
        }
    }
    cycles <- lapply(cycles, names)
    return(cycles)
}

.find_all_relabel_cycles <- function(relabels_graph) {
    relabels <- data.frame(
        relabel_from = character(0),
        relabel_to = character(0)
    )
    if (vcount(relabels_graph) == 0) {
        return(relabels)
    }

    ## A cycle can never span more than one weakly-connected component (a
    ## cycle is itself connected), so every component's disjoint cycles can
    ## be found independently of every other component's. The
    ## pre-optimization version of this function searched the whole graph
    ## at once: at every increasing cutoff, it re-scanned every remaining
    ## vertex in *every* component, up to the size of the single largest
    ## remaining component. On a graph with many small components -- the
    ## common case in practice, since most mislabel events are simple
    ## pairwise swaps, each its own tiny 2-vertex component -- that means
    ## thousands of irrelevant vertices get rescanned at every cutoff level
    ## for no benefit, just because one large component elsewhere hasn't
    ## finished yet. Splitting into components up front and searching each
    ## one on its own, only up to *its own* size (see
    ## .find_relabel_cycles_in_component()), finds exactly the same cycles
    ## -- just far faster. The only observable difference is the row order
    ## of the returned data frame (component-by-component here, vs.
    ## interleaved by cutoff level across all components before), which
    ## does not affect correctness.
    comp <- components(relabels_graph, mode = "weak")
    relabels_by_component <- vector("list", comp$no)
    for (comp_id in seq_len(comp$no)) {
        comp_vertices <- V(relabels_graph)[comp$membership == comp_id]
        if (length(comp_vertices) < 2) {
            next
        }
        subgraph <- induced_subgraph(relabels_graph, comp_vertices)
        relabels_by_component[[comp_id]] <- .find_relabel_cycles_in_component(
            subgraph
        )
    }

    return(do.call(rbind, c(list(relabels), relabels_by_component)))
}

## Finds the same greedy, shortest-cycles-first set of vertex-disjoint
## relabel cycles within a single weakly-connected component that the
## pre-optimization .find_all_relabel_cycles() found by repeatedly
## searching the *whole* graph -- see the comment there for why searching
## one component at a time gives identical results. Split out so it can be
## called once per component instead.
.find_relabel_cycles_in_component <- function(relabels_graph) {
    relabels <- data.frame(
        relabel_from = character(0),
        relabel_to = character(0)
    )

    ## Fast path for by far the most common case: a component that is
    ## exactly two vertices with an edge in each direction between them,
    ## i.e. a simple 2-sample swap. This is exactly the (only possible)
    ## cycle the general search below would otherwise find at cutoff = 1,
    ## so this is not a behavior change -- just a shortcut around the
    ## overhead of building a tiny igraph subgraph and calling
    ## all_simple_paths() on it, for what is typically the overwhelming
    ## majority of components.
    if (vcount(relabels_graph) == 2 && ecount(relabels_graph) == 2) {
        if (all(which_mutual(relabels_graph))) {
            vertex_names <- V(relabels_graph)$name
            return(data.frame(
                relabel_from = vertex_names,
                relabel_to = rev(vertex_names)
            ))
        }
    }

    all_relabeled_samples <- NULL
    all_cycles <- list()
    cutoff <- 1
    while (
        vcount(relabels_graph) > 0 &&
            cutoff < max(table(components(relabels_graph)$membership))
    ) {
        curr_cycles <- .find_directed_cycles(relabels_graph, cutoff = cutoff)
        for (curr_cycle in curr_cycles) {
            if (!any(curr_cycle %in% all_relabeled_samples)) {
                all_cycles <- append(all_cycles, list(curr_cycle))
                all_relabeled_samples <- c(all_relabeled_samples, curr_cycle)
                relabels_graph <- delete_vertices(relabels_graph, curr_cycle)
            }
        }
        cutoff <- cutoff + 1
    }

    ## From the found cycles, construct relabels dataframe
    for (curr_cycle in all_cycles) {
        n <- length(curr_cycle)
        curr_relabels <- data.frame(
            relabel_from = curr_cycle,
            relabel_to = c(curr_cycle[2:n], curr_cycle[1])
        )
        relabels <- rbind(relabels, curr_relabels)
    }

    return(relabels)
}

.find_relabel_cycles_from_putative_subjects <- function(
    unsolved_relabel_data,
    putative_subjects,
    unsolved_ghost_data = NULL,
    unambiguous_only = FALSE,
    allow_unknowns = FALSE
) {
    allow_ghosts <- !is.null(unsolved_ghost_data)

    if (unambiguous_only) {
        allow_ghosts <- FALSE
        allow_unknowns <- FALSE
    }

    mislabel_data <- unsolved_relabel_data |>
        left_join(
            putative_subjects |> rename(Putative_Subject_ID = "Subject_ID"),
            by = "Genotype_Group_ID"
        ) |>
        mutate(
            Inferred_Correctly_Labeled = .data$Putative_Subject_ID ==
                .data$Subject_ID,
            ## Only include samples where the current label is a Subject_ID
            ## that also has a Genotype_Group_ID assigned
            Curr_Subject_ID_Genotyped = .data$Subject_ID %in%
                putative_subjects$Subject_ID
        ) |>
        filter(
            !is.na(.data$Putative_Subject_ID) &
                !.data$Inferred_Correctly_Labeled &
                .data$Curr_Subject_ID_Genotyped
        )

    relabels <- EMPTY_RELABELS

    if (nrow(mislabel_data) == 0) {
        return(relabels)
    }

    ## Every mislabeled sample with the same genotype and same label
    ## domain must also have the same potential relabels. We search
    ## for relabels at at the Label_Domain/Genotype_Group_ID level
    mislabeled_genotype_label_domains <- mislabel_data |>
        select("Label_Domain", "Genotype_Group_ID") |>
        distinct() |>
        left_join(putative_subjects, by = "Genotype_Group_ID")

    directed_edge_mats <- vector(
        "list",
        length = nrow(mislabeled_genotype_label_domains)
    )

    all_ghost_labels <- character(0)
    all_unknown_labels <- character(0)

    ## The loop below used to re-filter 'mislabel_data', 'unsolved_relabel_data',
    ## and 'unsolved_ghost_data' from scratch (via dplyr::filter()) on every
    ## iteration -- fine for the handful of rows in the package's toy
    ## scenarios, but O(number of mislabeled Label_Domain/Genotype_Group_ID
    ## pairs * table size) on realistically large datasets, where both the
    ## number of pairs and the table sizes can be large simultaneously.
    ## Splitting each table by its grouping key(s) once up front, and doing
    ## plain list lookups inside the loop, computes exactly the same
    ## per-iteration values (same rows, same order, since split() preserves
    ## each group's original row order just like filter() does) in a small
    ## fraction of the time. The composite split keys below are built by
    ## pasting the grouping columns together with a control character
    ## ("\x1f", ASCII unit separator) that is not expected to ever appear in
    ## a real Label_Domain/Genotype_Group_ID/Subject_ID value, so there is no
    ## practical risk of two distinct groups colliding onto the same key.
    key_sep <- "\x1f"

    mislabel_group_key <- paste(
        mislabel_data$Label_Domain,
        mislabel_data$Genotype_Group_ID,
        sep = key_sep
    )
    mislabeled_samples_by_group <- split(
        mislabel_data$Sample_ID,
        mislabel_group_key
    )
    ## Pre-sort by Sample_ID (once, globally) so that splitting preserves,
    ## within each group, the same order `arrange(Sample_ID) |>
    ## pull(Placeholder_ID)` produced before.
    mislabel_sample_order <- order(mislabel_data$Sample_ID)
    mislabel_placeholders_by_group <- split(
        mislabel_data$Placeholder_ID[mislabel_sample_order],
        mislabel_group_key[mislabel_sample_order]
    )

    eligible_group_key <- paste(
        unsolved_relabel_data$Label_Domain,
        unsolved_relabel_data$Subject_ID,
        sep = key_sep
    )
    eligible_samples_by_group <- split(
        unsolved_relabel_data$Sample_ID,
        eligible_group_key
    )
    eligible_genotypes_by_group <- split(
        unsolved_relabel_data$Genotype_Group_ID,
        eligible_group_key
    )

    if (allow_ghosts) {
        ghost_group_key <- paste(
            unsolved_ghost_data$Label_Domain,
            unsolved_ghost_data$Subject_ID,
            sep = key_sep
        )
        ghost_samples_by_group <- split(
            unsolved_ghost_data$Sample_ID,
            ghost_group_key
        )
    }

    ## Iterate over all Label_Domain/Genotype_Group_ID pairs and search for
    ## potential relabels
    for (i in seq_len(nrow(mislabeled_genotype_label_domains))) {
        label_domain_id <- mislabeled_genotype_label_domains[i, "Label_Domain"]
        genotype_group_id <- mislabeled_genotype_label_domains[
            i,
            "Genotype_Group_ID"
        ]
        subject_id <- mislabeled_genotype_label_domains[i, "Subject_ID"]

        ## These are all the mislabeled samples for a Label_Domain/Genotype_Group_ID pair
        curr_mislabel_key <- paste(
            label_domain_id,
            genotype_group_id,
            sep = key_sep
        )
        mislabeled_samples <- mislabeled_samples_by_group[[curr_mislabel_key]]
        if (is.null(mislabeled_samples)) {
            mislabeled_samples <- character(0)
        }
        n_mislabeled_samples <- length(mislabeled_samples)

        ## The eligible relabels have the putative Subject_ID but are in a different Genotype_Group
        curr_subject_key <- paste(label_domain_id, subject_id, sep = key_sep)
        candidate_labels <- eligible_samples_by_group[[curr_subject_key]]
        candidate_genotypes <- eligible_genotypes_by_group[[curr_subject_key]]
        eligible_labels <- if (is.null(candidate_labels)) {
            character(0)
        } else {
            candidate_labels[candidate_genotypes != genotype_group_id]
        }
        n_eligible_labels <- length(eligible_labels)

        n_label_deficit <- n_mislabeled_samples - n_eligible_labels

        ghost_labels <- character(0)
        n_ghost_labels <- 0
        ## If there aren't enough eligible labels, try using ghosts to plug the gap
        if (allow_ghosts && n_label_deficit > 0) {
            ghost_labels <- ghost_samples_by_group[[curr_subject_key]]
            if (is.null(ghost_labels)) {
                ghost_labels <- character(0)
            }
            n_ghost_labels <- length(ghost_labels)

            ## If we have more ghost labels than needed, select a subset
            if (n_label_deficit < n_ghost_labels) {
                ghost_labels <- ghost_labels[seq_len(n_label_deficit)]
                n_ghost_labels <- length(ghost_labels)
            }

            n_label_deficit <- n_mislabeled_samples -
                n_eligible_labels -
                n_ghost_labels
        }

        unknown_labels <- character(0)
        n_unknown_labels <- 0
        ## If there still aren't enough eligible labels, resort to plugging the gap with
        ## unknowns. Reuse the Placeholder_ID(s) pre-generated at MislabelSolver()
        ## construction time (see MislabelSolver.R/.generate_placeholder_ids()) for a
        ## deterministic, Sample_ID-sorted subset of this group's own mislabeled samples,
        ## rather than minting new ones now -- so the resulting label is reproducible given
        ## the same input. Every mislabeled sample already has its own pre-generated
        ## Placeholder_ID regardless of whether it ends up used here.
        if (allow_unknowns && n_label_deficit > 0) {
            chosen_placeholder_ids <- head(
                mislabel_placeholders_by_group[[curr_mislabel_key]],
                n_label_deficit
            )
            unknown_labels <- str_c(
                LABEL_NOT_FOUND,
                "#",
                subject_id,
                "#",
                label_domain_id,
                "#",
                chosen_placeholder_ids
            )
            n_unknown_labels <- length(unknown_labels)
        }

        ## We should have n_mislabeled_samples == n_eligible_labels + n_ghost_labels + n_unknown_labels

        directed_edge_mats[[i]] <- expand.grid(
            relabel_from = mislabeled_samples,
            relabel_to = c(eligible_labels, ghost_labels, unknown_labels)
        )
        all_ghost_labels <- c(all_ghost_labels, ghost_labels)
        all_unknown_labels <- c(all_unknown_labels, unknown_labels)
    }

    ## This graph doesn't have outgoing edges from ghosts or unknowns
    ## So there are no cycles that include ghosts or unknowns yet
    directed_edge_df <- as.data.frame(do.call(rbind, directed_edge_mats))
    relabels_graph <- graph_from_data_frame(directed_edge_df, directed = TRUE)

    ## If unambiguous_only, only include samples with exactly one incoming and one outgoing edge
    if (unambiguous_only) {
        samples_with_one_incoming <- V(relabels_graph)[
            degree(relabels_graph, mode = "in") == 1
        ]$name
        samples_with_one_outgoing <- V(relabels_graph)[
            degree(relabels_graph, mode = "out") == 1
        ]$name
        samples_to_filter <- intersect(
            samples_with_one_incoming,
            samples_with_one_outgoing
        )
        relabels_graph <- subgraph(relabels_graph, samples_to_filter)
    }

    ## 1. Find relabel cycles without using ghosts or unknowns
    new_relabels <- .find_all_relabel_cycles(relabels_graph)
    relabels_graph <- relabels_graph - new_relabels$relabel_from
    relabels <- rbind(relabels, new_relabels)

    ## 2. Find relabel cycles allowing for ghosts
    if (allow_ghosts) {
        all_excess_labels <- names(V(relabels_graph)[
            degree(relabels_graph, mode = "in") == 0
        ])
        V(relabels_graph)$relabel_component_id <- components(
            relabels_graph
        )$membership
        all_relabel_component_ids <- unique(
            V(relabels_graph)$relabel_component_id
        )
        ## For each component, connect ghost labels with labels that have no more incoming edges,
        ## indicating that there are no genotyped samples that can take their label
        for (curr_relabel_component_id in all_relabel_component_ids) {
            component_labels <- names(V(relabels_graph)[
                V(relabels_graph)$relabel_component_id ==
                    curr_relabel_component_id
            ])
            component_ghost_labels <- intersect(
                all_ghost_labels,
                component_labels
            )
            component_excess_labels <- intersect(
                all_excess_labels,
                component_labels
            )
            new_edges <- expand.grid(
                relabel_from = component_ghost_labels,
                relabel_to = component_excess_labels
            )
            if (nrow(new_edges) > 0) {
                new_relabels_graph <- graph_from_data_frame(new_edges)
                relabels_graph <- igraph::union(
                    relabels_graph,
                    new_relabels_graph
                )
            }
        }
        new_relabels <- .find_all_relabel_cycles(relabels_graph)
        relabels_graph <- relabels_graph - new_relabels$relabel_from
        relabels <- rbind(relabels, new_relabels)
    }

    ## 3. Find relabel cycles allowing for unknowns
    if (allow_unknowns) {
        all_excess_labels <- names(V(relabels_graph)[
            degree(relabels_graph, mode = "in") == 0
        ])
        V(relabels_graph)$relabel_component_id <- components(
            relabels_graph
        )$membership
        all_relabel_component_ids <- unique(
            V(relabels_graph)$relabel_component_id
        )
        ## For each component, connect unknown labels with labels that have no more incoming edges,
        ## indicating that there are no genotyped samples that can take their label
        for (curr_relabel_component_id in all_relabel_component_ids) {
            component_labels <- names(V(relabels_graph)[
                V(relabels_graph)$relabel_component_id ==
                    curr_relabel_component_id
            ])
            component_unknown_labels <- intersect(
                all_unknown_labels,
                component_labels
            )
            component_excess_labels <- intersect(
                all_excess_labels,
                component_labels
            )
            new_edges <- expand.grid(
                relabel_from = component_unknown_labels,
                relabel_to = component_excess_labels
            )
            if (nrow(new_edges) > 0) {
                new_relabels_graph <- graph_from_data_frame(new_edges)
                relabels_graph <- igraph::union(
                    relabels_graph,
                    new_relabels_graph
                )
            }
        }
        new_relabels <- .find_all_relabel_cycles(relabels_graph)
        relabels_graph <- relabels_graph - new_relabels$relabel_from
        relabels <- rbind(relabels, new_relabels)
    }

    return(relabels)
}

.find_neighbors <- function(
    object,
    include_ghost = FALSE,
    filter_concordant_vertices = FALSE
) {
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
    unsolved_ghost_data <- NULL
    if (include_ghost) {
        unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data
    }
    combined_graph <- .generate_graph(
        unsolved_relabel_data,
        graph_type = "combined",
        unsolved_ghost_data
    )
    unsolved_all_data <- rbind(unsolved_relabel_data, unsolved_ghost_data)
    putative_subjects <- object@.solve_state$putative_subjects

    ## Criteria 1: filter only pairs of vertices that are within at
    ## exactly 2 edges of each other
    adj_matrix_sparse <- Matrix(as_adjacency_matrix(
        combined_graph,
        sparse = TRUE
    ))
    adj_matrix_idx <- Matrix::which(adj_matrix_sparse > 0, arr.ind = TRUE)
    dist_within_2_sparse <- adj_matrix_sparse %*% adj_matrix_sparse
    dist_within_2_idx <- Matrix::which(dist_within_2_sparse > 0, arr.ind = TRUE)
    unique_pairs_exact1 <- data.frame(
        Row = rownames(adj_matrix_sparse)[adj_matrix_idx[, 1]],
        Col = colnames(adj_matrix_sparse)[adj_matrix_idx[, 2]]
    ) |>
        transmute(
            Sample_A = pmin(.data$Row, .data$Col),
            Sample_B = pmax(.data$Row, .data$Col)
        ) |>
        filter(.data$Sample_A != .data$Sample_B) |>
        distinct()
    unique_pairs_within2 <- data.frame(
        Row = rownames(dist_within_2_sparse)[dist_within_2_idx[, 1]],
        Col = colnames(dist_within_2_sparse)[dist_within_2_idx[, 2]]
    ) |>
        transmute(
            Sample_A = pmin(.data$Row, .data$Col),
            Sample_B = pmax(.data$Row, .data$Col)
        ) |>
        filter(.data$Sample_A != .data$Sample_B) |>
        distinct()
    unique_pairs <- anti_join(
        unique_pairs_within2,
        unique_pairs_exact1,
        by = c("Sample_A", "Sample_B")
    )

    ## Criteria 2: filter out pairs of vertices that include elements
    ## in anchor_samples
    unique_pairs <- unique_pairs[
        !(unique_pairs$Sample_A %in%
            object@anchor_samples |
            unique_pairs$Sample_B %in% object@anchor_samples),
    ]

    ## Criteria 3: filter only pairs of vertices that are within the
    ## same label domain
    unique_pairs <- unique_pairs |>
        left_join(
            unsolved_all_data[, c("Sample_ID", "Label_Domain")],
            by = c("Sample_A" = "Sample_ID")
        ) |>
        rename("Label_Domain_A" = "Label_Domain") |>
        left_join(
            unsolved_all_data[, c("Sample_ID", "Label_Domain")],
            by = c("Sample_B" = "Sample_ID")
        ) |>
        rename("Label_Domain_B" = "Label_Domain") |>
        filter(.data$Label_Domain_A == .data$Label_Domain_B) |>
        select("Sample_A", "Sample_B")

    ## Criteria 4: filter out vertices that have at least 1 concordant edge
    if (filter_concordant_vertices) {
        concordant_edges <- E(combined_graph)[E(combined_graph)$concordant]
        concordant_vertices <- unique(c(ends(combined_graph, concordant_edges)))
        unique_pairs <- unique_pairs |>
            filter(
                !(.data$Sample_A %in% concordant_vertices),
                !(.data$Sample_B %in% concordant_vertices)
            )
    }

    ## Criteria 5: filter out pairs of vertices where either side will violate putative_subjects
    unique_pairs <- unique_pairs |>
        left_join(
            unsolved_all_data[, c(
                "Sample_ID",
                "Subject_ID",
                "Genotype_Group_ID"
            )],
            by = c("Sample_A" = "Sample_ID")
        ) |>
        rename(
            "Subject_A" = "Subject_ID",
            "Genotype_Group_A" = "Genotype_Group_ID"
        ) |>
        left_join(
            putative_subjects |> filter(!is.na(.data$Genotype_Group_ID)),
            by = c("Genotype_Group_A" = "Genotype_Group_ID")
        ) |>
        rename("Putative_Subject_A" = "Subject_ID") |>
        left_join(
            unsolved_all_data[, c(
                "Sample_ID",
                "Subject_ID",
                "Genotype_Group_ID"
            )],
            by = c("Sample_B" = "Sample_ID")
        ) |>
        rename(
            "Subject_B" = "Subject_ID",
            "Genotype_Group_B" = "Genotype_Group_ID"
        ) |>
        left_join(
            putative_subjects |> filter(!is.na(.data$Genotype_Group_ID)),
            by = c("Genotype_Group_B" = "Genotype_Group_ID")
        ) |>
        rename("Putative_Subject_B" = "Subject_ID") |>
        mutate(
            # The swap is invalid if
            # 1. Sample_A is already in its Putative_Subject
            # 2. Sample_B is already in its Putative_Subject
            # For swaps where the putative subjects of both samples are assigned, the swap is also invalid if
            # 3. Neither Sample_A nor Sample_B match their Putative_Subject after
            Invalid_Swap = (!is.na(.data$Putative_Subject_A) &
                .data$Subject_A == .data$Putative_Subject_A) |
                (!is.na(.data$Putative_Subject_B) &
                    .data$Subject_B == .data$Putative_Subject_B) |
                (!is.na(.data$Putative_Subject_A) &
                    !is.na(.data$Putative_Subject_B) &
                    .data$Putative_Subject_A != .data$Subject_B &
                    .data$Putative_Subject_B != .data$Subject_A)
        ) |>
        filter(!.data$Invalid_Swap) |>
        select("Sample_A", "Sample_B")

    return(unique_pairs)
}

## Extracted from solveLocalSearch()'s body (where both were previously
## defined as inline closures) since neither needs access to solveLocalSearch()'s
## other local variables beyond what's now passed explicitly as arguments.
##
## calc_scaled_entropy() only ever used its own argument `x`, so it moves out
## unchanged. calc_swapped_delta_entropy() referenced `votes` and
## `base_entropies` from solveLocalSearch()'s enclosing scope (recomputed
## each outer iteration from the object's current state); it now takes both
## as explicit arguments instead, matching the pattern already used by
## .calc_swapped_delta_entropy_fast() below.
.calc_scaled_entropy <- function(x) {
    # n <- sum(x)
    # return(-n*sum(x/n *log(x/n), na.rm=TRUE))
    # n <- sum(x)
    # return(n*sum(log(x/n), na.rm=TRUE))
    return(sum(x * log(x / sum(x)), na.rm = TRUE))
}

.calc_swapped_delta_entropy <- function(
    votes,
    base_entropies,
    swap_from_subject,
    swap_from_genotype,
    swap_to_subject,
    swap_to_genotype
) {
    delta <- 0
    if (!is.na(swap_from_genotype)) {
        genotype_base_entropy <- base_entropies[swap_from_genotype]
        genotype_votes_vec <- votes[swap_from_genotype, ]
        genotype_votes_vec[swap_from_subject] <- genotype_votes_vec[
            swap_from_subject
        ] -
            1
        genotype_votes_vec[swap_to_subject] <- genotype_votes_vec[
            swap_to_subject
        ] +
            1
        genotype_new_entropy <- .calc_scaled_entropy(genotype_votes_vec)
        delta <- delta + genotype_new_entropy - genotype_base_entropy
    }
    if (!is.na(swap_to_genotype)) {
        genotype_base_entropy <- base_entropies[swap_to_genotype]
        genotype_votes_vec <- votes[swap_to_genotype, ]
        genotype_votes_vec[swap_to_subject] <- genotype_votes_vec[
            swap_to_subject
        ] -
            1
        genotype_votes_vec[swap_from_subject] <- genotype_votes_vec[
            swap_from_subject
        ] +
            1
        genotype_new_entropy <- .calc_scaled_entropy(genotype_votes_vec)
        delta <- delta + genotype_new_entropy - genotype_base_entropy
    }
    return(delta)
}

## Closed-form, algebraically exact replacement for the per-swap entropy delta
## used by solveLocalSearchOld()'s calc_scaled_entropy()-based computation.
## "Algebraically" is an important qualifier here, not a hedge -- see the
## floating-point note below, and solveLocalSearch()'s documentation.
##
## For a genotype's vote vector x, calc_scaled_entropy(x) = sum_i x_i*log(x_i/n)
## = sum_i x_i*log(x_i) - n*log(n), where n = sum(x). A candidate swap moves
## one vote out of bucket `a`'s count and into bucket `b`'s count within the
## same genotype's row, leaving n unchanged, so n*log(n) cancels in the
## before/after difference. The change is exactly
##   xlogx(a-1) - xlogx(a) + xlogx(b+1) - xlogx(b)
## with xlogx(v) = v*log(v) for v>0 and 0 otherwise -- an O(1) lookup and a
## handful of scalar log() calls, regardless of how many other subjects exist
## in that row. This matters because solveLocalSearch()'s current cost scales
## with the *width* of the votes table (one column per distinct Subject_ID
## across the whole unsolved dataset, not just the candidate swap's own
## component), while this closed form does not.
##
## One case needs explicit handling: if a candidate swap's two samples already
## report the *same* Subject_ID, the original code's sequential
## decrement-then-increment on that single cell is a true no-op (it nets back
## to the same value); applying the two-cell formula naively there would
## treat one cell as two and give a nonzero answer, so that side's
## contribution is forced to 0 instead.
##
## Floating point: this is exact algebra (true for any real n, a, b under
## infinite precision) but NOT a bit-exact re-derivation of the original
## computation, because it evaluates log() on different arguments. The
## original computes log(x_i/n) for every bucket in the row and subtracts two
## full row-sums; this closed form computes log(a-1), log(a), log(b+1), log(b)
## directly and never materializes n at all (it cancels symbolically before
## any log() call happens, not numerically after). log(x/n) and
## log(x) - log(n) round differently in general, so the two paths can and do
## disagree in their last few bits -- confirmed empirically (candidate swaps
## built from real per-sample records, as .find_neighbors() would produce):
## about 3 in 5 evaluations differ from the original at the ~1e-14 to 1e-15
## (relative) level. That is normally inconsequential, but because
## solveLocalSearchOld() picks the single best swap via `delta == max(delta)`
## within a component, a difference this small can occasionally flip which
## swap wins a near-tie, which is why solveLocalSearch() is documented as
## a close-but-not-bit-exact replacement rather than an exact one. See that
## function's documentation for the practical implications.
##
## Verified against the original mapply(calc_swapped_delta_entropy, ...)
## computation on realistic synthetic data (candidate swaps derived from
## actual per-sample records, matching how .find_neighbors() produces them,
## rather than arbitrary label combinations) across component sizes from 60
## to 3,000 candidate subjects: identical results (all.equal tolerance 1e-9)
## in every case tested, 9-380x faster, with the gap widening as the subject
## pool grows. See solveLocalSearch().
.xlogx <- function(v) if_else(v <= 0, 0, v * log(v))

.calc_swapped_delta_entropy_fast <- function(
    votes_mat,
    swap_from_subject,
    swap_from_genotype,
    swap_to_subject,
    swap_to_genotype
) {
    n <- length(swap_from_subject)
    delta <- numeric(n)

    has_from <- !is.na(swap_from_genotype)
    same_from <- has_from & (swap_from_subject == swap_to_subject)
    do_from <- has_from & !same_from
    if (any(do_from)) {
        a <- votes_mat[cbind(
            swap_from_genotype[do_from],
            swap_from_subject[do_from]
        )]
        b <- votes_mat[cbind(
            swap_from_genotype[do_from],
            swap_to_subject[do_from]
        )]
        delta[do_from] <- delta[do_from] +
            (.xlogx(a - 1) - .xlogx(a)) +
            (.xlogx(b + 1) - .xlogx(b))
    }

    has_to <- !is.na(swap_to_genotype)
    same_to <- has_to & (swap_from_subject == swap_to_subject)
    do_to <- has_to & !same_to
    if (any(do_to)) {
        a <- votes_mat[cbind(swap_to_genotype[do_to], swap_to_subject[do_to])]
        b <- votes_mat[cbind(swap_to_genotype[do_to], swap_from_subject[do_to])]
        delta[do_to] <- delta[do_to] +
            (.xlogx(a - 1) - .xlogx(a)) +
            (.xlogx(b + 1) - .xlogx(b))
    }

    return(delta)
}

## Base-R vectorized replacement for the per-swap-category permutation-scoring
## block in solveGlobalSearch() (the loop that melts perm_genotypes,
## left_joins it against label/ghost/genotype/concordant count tables, and
## group_by/summarizes per Permutation_ID). Two independent changes, both
## exact (not approximations):
##
## 1. Locked/free separation: perm_genotypes has one column per genotype in
##    the component, but only the (at most max_genotypes) *free* columns
##    actually vary across permutations -- every *locked* column repeats the
##    same subject ID in all ~n_perms rows, because that genotype's identity
##    was already resolved before global search reached this
##    component. Locked columns can outnumber free ones considerably in a
##    large component. Their contribution to each permutation's score is
##    therefore a constant, computable once from a single row rather than
##    melted and joined ~n_perms times.
## 2. left_join()/group_by()/summarize() is replaced with direct
##    named-vector lookups and rowsum(): profiling showed the original
##    block's cost dominated by dplyr/vctrs per-call overhead (data-mask
##    construction, hash-join bookkeeping), not the arithmetic itself, so
##    removing that overhead matters at least as much as reducing row count.
##
## Unlike .calc_swapped_delta_entropy_fast() (see its own floating-point note
## above), this one really is bit-exact, not merely close, and that is
## provable rather than just empirically observed: every quantity summed here
## (n_labels, n_ghost_labels, n_in_genotype, n_samples_correct, and everything
## pmin()/pmax()/subtraction-derived from them) is a small non-negative
## integer count, and IEEE-754 addition/subtraction of integers is exact
## (no rounding at all) as long as the running total stays under 2^53 --
## utterly unremarkable for realistic sample sizes. With no rounding ever
## occurring in the count arithmetic, it does not matter that this function
## sums free and locked genotypes' contributions separately (rowsum() then
## sweep()) where solveGlobalSearch() sums them together in one pass:
## reordering exact values changes nothing. The only non-integer quantities
## are the final `ghost_penalty * n_samples_to_relabel_ghost` and
## `deletion_penalty * pmax(n_genotype_deletions, n_label_deletions)` terms,
## but those multiply the *same* exact integer by the *same* penalty value
## passed to both functions, combined via the identical left-to-right
## `a + b + c` expression and the identical across-swap-category accumulation
## order (both loop over cc_label_domain_ids, derived the same way from the same
## input, and accumulate via `permutation_stats <- permutation_stats + ...`)
## -- so there is no step anywhere in either function where the same
## arithmetic operation is ever applied to different operands, or the same
## operands combined in a different order. That holds for *any* valid
## ghost_penalty/deletion_penalty, not just the package defaults.
##
## Verified against the original dplyr-based computation on real captured
## mid-solve state (a component with 8 free / 16 locked genotypes, 40,320
## permutations, 6 swap categories): identical scores and identical best
## permutation in every case tested, ~16x faster on that case (~1.2x on a
## small synthetic case with 0 locked genotypes, where fixed overhead
## dominates either way and there is nothing to separate out). See
## solveGlobalSearch().
.score_permutations_fast <- function(
    perm_genotypes,
    free_genotypes,
    locked_genotypes,
    cc_label_domain_ids,
    label_counts,
    ghost_label_counts,
    genotype_counts,
    genotype_subject_concordant_counts,
    ghost_penalty = 1.5,
    deletion_penalty = 4
) {
    n_perms <- nrow(perm_genotypes)
    permutation_ids <- rownames(perm_genotypes)
    sum_cols <- c(
        "n_samples_correct",
        "n_samples_to_relabel",
        "n_samples_to_relabel_ghost",
        "n_genotype_deletions",
        "n_label_deletions"
    )

    free_cols <- intersect(colnames(perm_genotypes), free_genotypes)
    locked_cols <- setdiff(colnames(perm_genotypes), free_cols)
    has_free <- length(free_cols) > 0
    has_locked <- length(locked_cols) > 0

    if (has_free) {
        long_free <- melt(perm_genotypes[, free_cols, drop = FALSE])
        colnames(long_free) <- c(
            "Permutation_ID",
            "Genotype_Group_ID",
            "Subject_ID"
        )
        long_free$Permutation_ID <- as.character(long_free$Permutation_ID)
        long_free$Genotype_Group_ID <- as.character(long_free$Genotype_Group_ID)
        long_free$Subject_ID <- as.character(long_free$Subject_ID)
    }
    if (has_locked) {
        locked_row <- perm_genotypes[1, locked_cols, drop = FALSE]
        long_locked <- data.frame(
            Genotype_Group_ID = colnames(locked_row),
            Subject_ID = as.character(locked_row[1, ]),
            stringsAsFactors = FALSE
        )
    }

    ## key -> value lookup as a named vector; missing keys resolve to 0 (matching
    ## the original's coalesce(., 0) after a non-matching left_join)
    vec_lookup <- function(df, key_cols, val_col) {
        key <- if (length(key_cols) == 1) {
            df[[key_cols]]
        } else {
            do.call(paste, c(df[key_cols], sep = "\x1f"))
        }
        setNames(df[[val_col]], key)
    }
    lookup0 <- function(vec, key) {
        v <- unname(vec[key])
        v[is.na(v)] <- 0
        v
    }

    permutation_stats <- matrix(
        0,
        nrow = n_perms,
        ncol = length(sum_cols) + 1,
        dimnames = list(permutation_ids, c(sum_cols, "perm_score"))
    )

    for (label_domain_id in cc_label_domain_ids) {
        v_label <- vec_lookup(
            label_counts[label_counts$Label_Domain == label_domain_id, ],
            "Subject_ID",
            "n_labels"
        )
        v_ghost <- vec_lookup(
            ghost_label_counts[
                ghost_label_counts$Label_Domain == label_domain_id,
            ],
            "Subject_ID",
            "n_ghost_labels"
        )
        v_geno <- vec_lookup(
            genotype_counts[genotype_counts$Label_Domain == label_domain_id, ],
            "Genotype_Group_ID",
            "n_in_genotype"
        )
        v_conc <- vec_lookup(
            genotype_subject_concordant_counts[
                genotype_subject_concordant_counts$Label_Domain ==
                    label_domain_id,
            ],
            c("Subject_ID", "Genotype_Group_ID"),
            "n_samples_correct"
        )

        score_rows <- function(df) {
            n_labels <- lookup0(v_label, df$Subject_ID)
            n_ghost_labels <- lookup0(v_ghost, df$Subject_ID)
            n_in_genotype <- lookup0(v_geno, df$Genotype_Group_ID)
            n_samples_correct <- lookup0(
                v_conc,
                str_c(df$Subject_ID, df$Genotype_Group_ID, sep = "\x1f")
            )
            n_samples_to_relabel <- pmin(n_in_genotype, n_labels) -
                n_samples_correct
            n_samples_to_relabel_ghost <- pmin(
                n_in_genotype - n_samples_correct - n_samples_to_relabel,
                n_ghost_labels
            )
            n_label_deletions <- pmax(
                0,
                n_in_genotype - n_labels - n_ghost_labels
            )
            n_genotype_deletions <- pmax(0, n_labels - n_in_genotype)
            cbind(
                n_samples_correct,
                n_samples_to_relabel,
                n_samples_to_relabel_ghost,
                n_genotype_deletions,
                n_label_deletions
            )
        }

        if (has_free) {
            free_scored <- score_rows(long_free)
            free_sums <- rowsum(
                free_scored,
                group = long_free$Permutation_ID,
                reorder = FALSE
            )
            free_sums <- free_sums[permutation_ids, , drop = FALSE]
        } else {
            free_sums <- matrix(
                0,
                nrow = n_perms,
                ncol = length(sum_cols),
                dimnames = list(permutation_ids, sum_cols)
            )
        }

        if (has_locked) {
            locked_totals <- colSums(score_rows(long_locked))
        } else {
            locked_totals <- setNames(rep(0, length(sum_cols)), sum_cols)
        }

        total <- sweep(free_sums, 2, locked_totals[colnames(free_sums)], "+")
        n_samples_to_relabel <- total[, "n_samples_to_relabel"] +
            pmin(
                total[, "n_genotype_deletions"],
                total[, "n_samples_to_relabel_ghost"]
            )
        n_genotype_deletions <- pmax(
            0,
            total[, "n_genotype_deletions"] -
                total[, "n_samples_to_relabel_ghost"]
        )
        perm_score <- n_samples_to_relabel +
            ghost_penalty * total[, "n_samples_to_relabel_ghost"] +
            deletion_penalty *
                pmax(n_genotype_deletions, total[, "n_label_deletions"])
        label_domain_total <- cbind(
            n_samples_correct = total[, "n_samples_correct"],
            n_samples_to_relabel,
            n_samples_to_relabel_ghost = total[, "n_samples_to_relabel_ghost"],
            n_genotype_deletions,
            n_label_deletions = total[, "n_label_deletions"],
            perm_score
        )
        ## When n_perms == 1 (a component with a single free genotype, so only
        ## one arrangement is possible), R's `[` silently drops rownames on
        ## every single-row extraction above (matrix[, "col"] returns an
        ## unnamed scalar when the matrix has exactly one row, unlike the
        ## multi-row case where rownames carry through as names) -- so
        ## label_domain_total ends up with no rownames at all despite having the
        ## right values, and the reindex below then fails with "subscript out
        ## of bounds". Row order is never changed by sweep()/cbind() above, so
        ## it's always safe to relabel directly from permutation_ids.
        ##
        ## Considered using `drop = FALSE` on the six `total[, "col"]`
        ## extractions above instead of patching rownames here, but that
        ## turns out to be *more* fragile, not less: `pmax()`/`pmin()` silently
        ## drop matrix-ness (and with it, dimnames) whenever one argument is a
        ## bare scalar -- e.g. `pmax(0, total[, "x", drop = FALSE])` returns a
        ## plain unnamed vector, not a 1-column matrix, even though the input
        ## was explicitly kept 2-D -- and both `n_samples_to_relabel` and
        ## `n_genotype_deletions` above are combined with a bare `0` via
        ## `pmax()`. Getting drop=FALSE to actually work end-to-end here would
        ## require rewriting those pmax()/pmin() calls too, for no benefit
        ## over just reasserting the (known-correct) rownames once, in one
        ## place, right before they matter.
        rownames(label_domain_total) <- permutation_ids
        permutation_stats <- permutation_stats +
            label_domain_total[
                rownames(permutation_stats),
                colnames(permutation_stats)
            ]
    }

    permutation_stats |>
        as.data.frame() |>
        rownames_to_column("Permutation_ID") |>
        arrange(.data$perm_score)
}
