#' Majority-based Sample Relabeling
#'
#' This heuristic algorithm aims to quickly correct "low-hanging fruit" mislabel
#' situations, in which a majority of samples for each subject are labeled
#' correctly. First, the algorithm "locks" in a Genotype_Group_ID <-> Subject_ID
#' pairing wherever the two agree on who the majority partner is: a majority of
#' the samples currently claiming a given Genotype_Group_ID must share the same
#' Subject_ID, *and* a majority of the samples currently claiming that
#' Subject_ID must share that same Genotype_Group_ID. Samples that disagree with
#' a locked pairing are then relabeled to match it.
#'
#' This is the quickest solver, because it never needs to search or score
#' alternatives, only count votes. Conversely, it is the least flexible, so it
#' relies on other solvers (*via* [solveEnsemble()]) to handle the complex cases
#' that it cannot.
#'
#' Use via [solveEnsemble()] rather than calling this directly, unless you are
#' composing your own solver loop.
#'
#' @param object A MislabelSolver object
#' @param unambiguous_only If TRUE, restrict relabeling to only the clearest
#'   cases. Candidate relabels that would require inventing a placeholder
#'   ("duplicate that doesn't exist") label or using up a ghost (ungenotyped)
#'   sample are skipped, as are any samples that have more than one possible
#'   real relabel target to choose between. Anything skipped this way is left
#'   for a later solver to attempt instead.
#' @param max_genotypes Components of the mislabel network with more
#'   not-yet-resolved Genotype_Group_ID(s) or Subject_ID(s) than this are
#'   skipped entirely for this call, the same way [solveGlobalSearch()]'s own
#'   `max_genotypes` argument works. Unlike global search, majority search's
#'   cost per component doesn't grow factorially with its size, so it can
#'   afford a much larger default; this mainly guards against wasting time on
#'   the rare intractably-large component in extremely large inputs.
#'
#' @return A MislabelSolver object
#'
#' @seealso [solveEnsemble()], [solveGlobalSearch()], [solveLocalSearch()]
#'
#' @export
#'
solveMajoritySearch <- function(
    object,
    unambiguous_only = FALSE,
    max_genotypes = 100
) {
    tsmsg("Starting majority search")
    if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
        tsmsg("0 samples relabeled")
        return(object)
    }

    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
    putative_subjects <- object@.solve_state$putative_subjects
    skip_component_ids <- .majority_search_skip_components(
        unsolved_relabel_data,
        putative_subjects,
        max_genotypes
    )
    if (length(skip_component_ids) > 0) {
        unsolved_relabel_data <- unsolved_relabel_data |>
            filter(!(.data$Component_ID %in% skip_component_ids))
    }
    if (nrow(unsolved_relabel_data) == 0) {
        tsmsg("0 samples relabeled")
        return(object)
    }

    ## 1. Update putative subjects
    votes <- table(
        unsolved_relabel_data$Genotype_Group_ID,
        unsolved_relabel_data$Subject_ID
    )
    votes_by_genotype <- data.frame(
        Genotype_Group_ID = rownames(votes),
        Max_Subject_ID = colnames(votes)[apply(votes, 1, which.max)],
        n = rowSums(votes),
        n_Max_Subject_ID = apply(votes, 1, max)
    ) |>
        filter(
            # Multiple subjects to choose from
            .data$n_Max_Subject_ID >= 2,
            # One subject has a majority
            .data$n_Max_Subject_ID > n / 2
        ) |>
        rename(Subject_ID = "Max_Subject_ID") |>
        select("Genotype_Group_ID", "Subject_ID")
    votes_by_subject <- data.frame(
        Subject_ID = colnames(votes),
        Max_Genotype_Group_ID = rownames(votes)[apply(votes, 2, which.max)],
        n = colSums(votes),
        n_Max_Genotype_Group_ID = apply(votes, 2, max)
    ) |>
        filter(
            # Multiple genotypes to choose from
            .data$n_Max_Genotype_Group_ID >= 2,
            # One genotype has a majority
            .data$n_Max_Genotype_Group_ID > n / 2
        ) |>
        rename(Genotype_Group_ID = "Max_Genotype_Group_ID") |>
        select("Subject_ID", "Genotype_Group_ID")
    # Inner join gives mutual best matches
    new_putative_subjects <- inner_join(
        votes_by_subject,
        votes_by_genotype,
        by = c("Genotype_Group_ID", "Subject_ID")
    ) |>
        anti_join(
            object@.solve_state$putative_subjects,
            by = c("Genotype_Group_ID", "Subject_ID")
        )
    object <- .update_putative_subjects(object, new_putative_subjects)

    ## 2. Find relabel cycles
    putative_subjects <- object@.solve_state$putative_subjects
    relabels <- .find_relabel_cycles_from_putative_subjects(
        unsolved_relabel_data,
        putative_subjects,
        unambiguous_only = unambiguous_only,
        allow_unknowns = FALSE
    )

    ## 3. Relabel samples and update solve state
    object <- .relabel_samples(object, relabels, solver_name = "majority")
    return(object)
}

#' Global Search Sample Relabeling
#'
#' This global search function permutes over all possible assignments of
#' Subject_IDs to Genotype_Group_IDs, then picks the assignment that implies the
#' fewest number of sample mislabels and deletions. As such, it is guaranteed to
#' find the optimal assignment given the selected scoring penalties. However,
#' since the running time is factorial in the number of genotype groups, it is
#' limited to only small connected components (by default, those involving no
#' more than 8 genotype groups).
#'
#' When choosing penalty values, keep in mind that the penalty for a relabeled
#' normal sample is always 1, and the other penalties are set relative to this.
#'
#' Use via [solveEnsemble()] rather than calling this directly, unless you are
#' composing your own solver loop.
#'
#' @param object A MislabelSolver object
#' @param max_genotypes The number of combinations scales in factorial with the
#'   number of genotypes in the largest connected component. The algorithm will
#'   skip over all components that exceed this size.
#' @param ghost_penalty The penalty for relabeling one ghost (placeholder)
#'   sample rather than a real one. Must be positive number. The ghost penalty
#'   should be strictly greater than 1, since otherwise the algorithm would
#'   prefer relabeling to a ghost sample even when a valid non-ghost sample swap
#'   is already available.
#' @param deletion_penalty The penalty for a label or genotype deletion (giving
#'   up on reconciling that sample entirely). Must be a positive number. The
#'   deletion penalty should be greater than twice the relabel penalty (i.e.
#'   greater than 2) *and* greater `ghost_penalty`, since otherwise the
#'   algorithm may prefer deleting samples even when a valid non-ghost sample
#'   swap or ghost sample is already available.
#'
#' @return A MislabelSolver object, potentially with some samples relabeled.
#'
#' @export
#'
solveGlobalSearch <- function(
    object,
    max_genotypes = 8,
    ghost_penalty = 1.5,
    deletion_penalty = 4
) {
    tsmsg("Starting global search")
    .validate_search_penalties(ghost_penalty, deletion_penalty)
    if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
        tsmsg("0 samples relabeled")
        return(object)
    }

    putative_subjects <- object@.solve_state$putative_subjects

    ## 1. Update putative subjects
    component_ids <- sort(unique(
        object@.solve_state$unsolved_relabel_data$Component_ID
    ))
    for (component_id in component_ids) {
        cc_unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data |>
            filter(.data$Component_ID == component_id)
        cc_unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data |>
            filter(.data$Component_ID == component_id)
        cc_sample_ids <- c(
            cc_unsolved_relabel_data$Sample_ID,
            cc_unsolved_ghost_data$Sample_ID
        )
        cc_label_domain_ids <- unique(c(
            cc_unsolved_relabel_data$Label_Domain,
            cc_unsolved_ghost_data$Label_Domain
        ))
        cc_genotypes <- unique(cc_unsolved_relabel_data$Genotype_Group_ID)
        cc_subjects <- unique(cc_unsolved_relabel_data$Subject_ID)

        ## For now, pass out of components where number of Genotype_Group(s) is
        ## greater than number of Subject_ID(s). In this case there will be one
        ## GG with no assigned subject ID.
        if (length(cc_genotypes) > length(cc_subjects)) {
            next
        }

        ## Lock genotypes that already have a putative subject assigned, and
        ## find all possible permutations for free genotypes
        locked_genotypes <- intersect(
            putative_subjects$Genotype_Group_ID,
            cc_genotypes
        )
        locked_subjects <- intersect(putative_subjects$Subject_ID, cc_subjects)
        free_genotypes <- setdiff(cc_genotypes, locked_genotypes)
        free_subjects <- setdiff(cc_subjects, locked_subjects)

        ## Pass out of components if they have too many genotypes
        if (
            max(length(free_genotypes), length(free_subjects)) > max_genotypes
        ) {
            next
        }

        if (length(free_genotypes) > 0 && length(free_subjects) > 0) {
            n <- length(free_subjects)
            r <- length(free_genotypes)
            perm_genotypes <- permutations(n, r, free_subjects)
            colnames(perm_genotypes) <- sort(free_genotypes)
            n_perm <- nrow(perm_genotypes)
            for (locked_genotype_id in locked_genotypes) {
                locked_subject_id <- putative_subjects[
                    putative_subjects$Genotype_Group_ID == locked_genotype_id,
                    "Subject_ID"
                ][[1]]
                new_perm_col <- matrix(
                    data = locked_subject_id,
                    ncol = 1,
                    nrow = n_perm,
                    dimnames = list(NULL, locked_genotype_id)
                )
                perm_genotypes <- cbind(perm_genotypes, new_perm_col)
            }
        } else {
            locked_putative_subjects <- putative_subjects[
                putative_subjects$Genotype_Group_ID %in% locked_genotypes,
            ]
            perm_genotypes <- t(locked_putative_subjects$Subject_ID)
            colnames(
                perm_genotypes
            ) <- locked_putative_subjects$Genotype_Group_ID
        }
        perm_genotypes <- as.matrix(
            perm_genotypes,
            dimnames = c("Permutation_ID", "Genotype_Group_ID")
        )
        n_perms <- nrow(perm_genotypes)
        permutation_ids <- str_c(
            "Permutation",
            formatC(
                seq_len(n_perms),
                width = str_length(n_perms),
                format = "d",
                flag = "0"
            )
        )
        rownames(perm_genotypes) <- permutation_ids

        ## For each Genotype_Group_ID/Subject_ID permutation, determine
        ## 1. The number of existing samples to relabel
        ## 2. The number of ghost samples needed to add
        ## 3. The number of indels required after ghost samples are included
        label_counts <- cc_unsolved_relabel_data |>
            select(
                "Sample_ID",
                "Subject_ID",
                "Genotype_Group_ID",
                "Label_Domain"
            ) |>
            group_by(.data$Subject_ID, .data$Label_Domain) |>
            summarize(n_labels = n(), .groups = "drop")
        ghost_label_counts <- cc_unsolved_ghost_data |>
            select(
                "Sample_ID",
                "Subject_ID",
                "Genotype_Group_ID",
                "Label_Domain"
            ) |>
            group_by(.data$Subject_ID, .data$Label_Domain) |>
            summarize(n_ghost_labels = n(), .groups = "drop")
        genotype_counts <- cc_unsolved_relabel_data |>
            select(
                "Sample_ID",
                "Subject_ID",
                "Genotype_Group_ID",
                "Label_Domain"
            ) |>
            group_by(.data$Genotype_Group_ID, .data$Label_Domain) |>
            summarize(n_in_genotype = n(), .groups = "drop")
        genotype_subject_concordant_counts <- cc_unsolved_relabel_data |>
            select(
                "Sample_ID",
                "Subject_ID",
                "Genotype_Group_ID",
                "Label_Domain"
            ) |>
            group_by(
                .data$Subject_ID,
                .data$Genotype_Group_ID,
                .data$Label_Domain
            ) |>
            summarize(n_samples_correct = n(), .groups = "drop")

        ## Base-R vectorized locked/free-aware scoring; see
        ## .score_permutations_fast() in helpers-solve.R
        permutation_stats <- .score_permutations_fast(
            perm_genotypes,
            free_genotypes,
            locked_genotypes,
            cc_label_domain_ids,
            label_counts,
            ghost_label_counts,
            genotype_counts,
            genotype_subject_concordant_counts,
            ghost_penalty = ghost_penalty,
            deletion_penalty = deletion_penalty
        )

        best_permutation <- perm_genotypes[
            permutation_stats$Permutation_ID[1],
            ,
            drop = FALSE
        ]

        if (nrow(permutation_stats) > 1) {
            best_score <- permutation_stats$perm_score[1]
            tied_permutation_stats <- permutation_stats |>
                filter(.data$perm_score == best_score)
            if (nrow(tied_permutation_stats) > 1) {
                tied_permutations <- perm_genotypes[
                    tied_permutation_stats$Permutation_ID,
                    ,
                    drop = FALSE
                ]
                for (curr_genotype_group in colnames(tied_permutations)) {
                    if (
                        !(curr_genotype_group %in%
                            names(object@.solve_state$ambiguous_subjects))
                    ) {
                        object@.solve_state$ambiguous_subjects[[
                            curr_genotype_group
                        ]] <- tied_permutations[, curr_genotype_group]
                    } else {
                        object@.solve_state$ambiguous_subjects[[
                            curr_genotype_group
                        ]] <-
                            unique(c(
                                object@.solve_state$ambiguous_subjects[[
                                    curr_genotype_group
                                ]],
                                tied_permutations[curr_genotype_group]
                            ))
                    }
                }
            }
        }

        new_putative_subjects <- best_permutation |>
            t() |>
            as.data.frame() |>
            rownames_to_column("X") |>
            relocate("X")
        colnames(new_putative_subjects) <- c("Genotype_Group_ID", "Subject_ID")
        rownames(new_putative_subjects) <- NULL

        ## Also update putative_subjects when a Subject_ID in the component
        ## doesn't have a Genotype_Group_ID, or vice versa
        if (length(cc_genotypes) > length(cc_subjects)) {
            unmatched_genotypes <- setdiff(
                cc_genotypes,
                new_putative_subjects$Genotype_Group_ID
            )
            new_putative_subjects <- rbind(
                new_putative_subjects,
                data.frame(
                    Genotype_Group_ID = unmatched_genotypes,
                    Subject_ID = NA_character_
                )
            )
        }
        if (length(cc_genotypes) < length(cc_subjects)) {
            unmatched_subjects <- setdiff(
                cc_subjects,
                new_putative_subjects$Subject_ID
            )
            new_putative_subjects <- rbind(
                new_putative_subjects,
                data.frame(
                    Genotype_Group_ID = NA_character_,
                    Subject_ID = unmatched_subjects
                )
            )
        }
        new_putative_subjects <- new_putative_subjects |>
            anti_join(
                object@.solve_state$putative_subjects,
                by = c("Genotype_Group_ID", "Subject_ID")
            )
        object <- .update_putative_subjects(object, new_putative_subjects)
    }

    ## Find relabel cycles
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
    unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data
    putative_subjects <- object@.solve_state$putative_subjects
    relabels <- .find_relabel_cycles_from_putative_subjects(
        unsolved_relabel_data,
        putative_subjects,
        unsolved_ghost_data,
        allow_unknowns = TRUE
    )

    ## Relabel samples and update solve state
    object <- .relabel_samples(object, relabels, solver_name = "global")
    return(object)
}

#' Local Search Sample Relabeling
#'
#' This search function looks through all possible swaps of 2 samples, and
#' selects the swap that minimizes the sum of within-genotype entropies within
#' each connected component of the mislabel network. Each iteration relabels up
#' to one pair of samples within each component.
#'
#' Each candidate swap's entropy delta is evaluated with a closed-form update.
#' The old, much slower algorithm, which rebuilds and rescans the full vote
#' vector for the affected genotype(s), is still available in
#' [solveLocalSearchOld()]. The new algorithm is 9-380x faster on synthetic
#' components ranging from 60 to 3,000 candidate subjects in testing.
#'
#' The closed-form update is algebraically equivalent to
#' [solveLocalSearchOld()]'s entropy formula, but it can produce numerically
#' different resuts due to the finite presicion of floating point numbers.
#' Generally the difference is less than 1e-9. This can cause the two functions
#' to choose different "best" swaps when two or more candidates are tied or
#' nearly tied. If you need results that exactly match [solveLocalSearchOld()],
#' use that function instead.
#'
#' Use via [solveEnsemble()] rather than calling this directly, unless you are
#' composing your own solver loop.
#'
#' @param object A MislabelSolver object
#' @param n_iter The number of hill-climbing iterations to run within this
#'   single call. Each iteration recomputes candidate swaps from the current
#'   (possibly already-updated, by an earlier iteration in this same call)
#'   state, and applies the single best non-conflicting swap per connected
#'   component. This returns early, before using up every requested iteration,
#'   only once every sample becomes solved, or no further entropy-improving
#'   swaps can be found.
#' @param include_ghost If TRUE, allow swaps to include ghost samples
#' @param filter_concordant_vertices If TRUE, filter out samples with at least
#'   one concordant edge. This reduces the search space but misses cases where
#'   multiple pairs of samples are swapped between the same subjects.
#' @param min_genotypes Components of the mislabel network with fewer
#'   not-yet-resolved Genotype_Group_ID(s) than this are skipped entirely for
#'   this call. Defaults to 1, which has no effect, since every component has
#'   at least 1 unsolved genotype group by definition. [solveEnsemble()]
#'   raises this above 1 so that local search only works on components too
#'   large for the other active solvers to have already handled themselves,
#'   rather than duplicating (and potentially disturbing) their work; see its
#'   documentation for the exact value used.
#'
#' @return A MislabelSolver object, potentially with some samples relabeled.
#'
#' @seealso [solveLocalSearchOld()], the original (exact) version of this
#'   algorithm, retained for cases needing bit-exact reproducibility
#'
#' @export
#'
solveLocalSearch <- function(
    object,
    n_iter = 1,
    include_ghost = FALSE,
    filter_concordant_vertices = FALSE,
    min_genotypes = 1
) {
    tsmsg("Starting local search")

    for (i in 1:n_iter) {
        tsmsg(paste(
            "Local search iteration (",
            i,
            " of ",
            n_iter,
            "):: 'include_ghost'=",
            include_ghost,
            ", 'filter_concordant_vertices'=",
            filter_concordant_vertices,
            sep = ""
        ))
        unsolved_all_data <- rbind(
            object@.solve_state$unsolved_relabel_data,
            object@.solve_state$unsolved_ghost_data
        )
        if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
            # All components already solved
            tsmsg("0 samples relabeled")
            return(object)
        }
        votes <- table(
            object@.solve_state$unsolved_relabel_data$Genotype_Group_ID,
            object@.solve_state$unsolved_relabel_data$Subject_ID
        )

        neighbors <- .find_neighbors(
            object,
            include_ghost,
            filter_concordant_vertices
        ) |>
            left_join(
                unsolved_all_data[, c(
                    "Sample_ID",
                    "Subject_ID",
                    "Genotype_Group_ID"
                )],
                by = c("Sample_A" = "Sample_ID")
            ) |>
            rename(
                Subject_A = "Subject_ID",
                Genotype_Group_A = "Genotype_Group_ID"
            ) |>
            left_join(
                unsolved_all_data[, c(
                    "Sample_ID",
                    "Subject_ID",
                    "Genotype_Group_ID",
                    "Component_ID"
                )],
                by = c("Sample_B" = "Sample_ID")
            ) |>
            rename(
                Subject_B = "Subject_ID",
                Genotype_Group_B = "Genotype_Group_ID"
            )

        tsmsg(paste(nrow(neighbors), "candidate swaps being evaluated..."))
        all_component_ids <- sort(unique(
            object@.solve_state$unsolved_relabel_data$Component_ID
        ))
        ## Precomputed once per iteration so the min_genotypes check below is
        ## a cheap lookup rather than a fresh filter() per component.
        unsolved_genotype_ids_by_component <- split(
            object@.solve_state$unsolved_relabel_data$Genotype_Group_ID,
            object@.solve_state$unsolved_relabel_data$Component_ID
        )
        relabels <- data.frame(matrix(
            data = NA,
            nrow = length(all_component_ids),
            ncol = 2,
            dimnames = list(c(), c("relabel_from", "relabel_to"))
        ))
        curr_idx <- 1
        for (curr_component_id in all_component_ids) {
            n_unsolved_genotypes <- length(unique(
                unsolved_genotype_ids_by_component[[as.character(
                    curr_component_id
                )]]
            ))
            if (n_unsolved_genotypes < min_genotypes) {
                next
            }

            cc_neighbors <- neighbors |>
                filter(.data$Component_ID == curr_component_id)

            if (nrow(cc_neighbors) == 0) {
                next
            }
            ## Closed-form vectorized delta instead of mapply(.calc_swapped_delta_entropy, ...);
            ## see .calc_swapped_delta_entropy_fast() in helpers-solve.R for the derivation.
            cc_neighbor_objectives <- cc_neighbors |>
                mutate(
                    delta = .calc_swapped_delta_entropy_fast(
                        votes,
                        .data$Subject_A,
                        .data$Genotype_Group_A,
                        .data$Subject_B,
                        .data$Genotype_Group_B
                    )
                )
            cc_relabels <- cc_neighbor_objectives |>
                filter(.data$delta > 0, .data$delta == max(.data$delta))
            if (nrow(cc_relabels) == 0) {
                next
            }
            cc_relabels <- cc_relabels |>
                sample_n(1) |>
                transmute(
                    relabel_from = .data$Sample_A,
                    relabel_to = .data$Sample_B
                )
            relabels[curr_idx, c("relabel_from", "relabel_to")] <- cc_relabels
            curr_idx <- curr_idx + 1
        }

        relabels <- relabels[!is.na(relabels[, 1]) & !is.na(relabels[, 2]), ]
        relabels <- rbind(
            relabels,
            data.frame(
                relabel_from = relabels$relabel_to,
                relabel_to = relabels$relabel_from
            )
        )
        object <- .relabel_samples(object, relabels, solver_name = "local")
        if (nrow(relabels) == 0) {
            # If nothing was relabeled in this iteration, then the algorithm has
            # converged and further iterations are not needed.
            break
        }
    }

    return(object)
}

#' Local Search Sample Relabeling (Original Algorithm)
#'
#' This search function looks through all possible swaps of 2 samples, and
#' selects the swap that minimizes the sum of within-genotype entropies.
#'
#' This is the original local search algorithm. It is superseded by
#' [solveLocalSearch()], which implements an optimized version of the same
#' algorithm. The new algorithm is algebraically but not numerically equivalent,
#' so the old algorithm is still provided. It should only be used if you need to
#' reproduce results generated using the old algorithm.
#'
#' @inheritParams solveLocalSearch
#'
#' @return A MislabelSolver object, potentially with some samples relabeled.
#'
#' @seealso [solveLocalSearch()], the standard (faster) version of this
#'   algorithm used by default.
#'
#' @export
#'
solveLocalSearchOld <- function(
    object,
    n_iter = 1,
    include_ghost = FALSE,
    filter_concordant_vertices = FALSE
) {
    tsmsg("Starting local search (old)")

    for (i in 1:n_iter) {
        tsmsg(paste(
            "Local search (old) iteration (",
            i,
            " of ",
            n_iter,
            "):: 'include_ghost'=",
            include_ghost,
            ", 'filter_concordant_vertices'=",
            filter_concordant_vertices,
            sep = ""
        ))
        unsolved_all_data <- rbind(
            object@.solve_state$unsolved_relabel_data,
            object@.solve_state$unsolved_ghost_data
        )
        if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
            tsmsg("0 samples relabeled")
            return(object)
        }
        votes <- table(
            object@.solve_state$unsolved_relabel_data$Genotype_Group_ID,
            object@.solve_state$unsolved_relabel_data$Subject_ID
        )
        base_entropies <- apply(votes, MARGIN = 1, .calc_scaled_entropy)

        neighbors <- .find_neighbors(
            object,
            include_ghost,
            filter_concordant_vertices
        ) |>
            left_join(
                unsolved_all_data[, c(
                    "Sample_ID",
                    "Subject_ID",
                    "Genotype_Group_ID"
                )],
                by = c("Sample_A" = "Sample_ID")
            ) |>
            rename(
                Subject_A = "Subject_ID",
                Genotype_Group_A = "Genotype_Group_ID"
            ) |>
            left_join(
                unsolved_all_data[, c(
                    "Sample_ID",
                    "Subject_ID",
                    "Genotype_Group_ID",
                    "Component_ID"
                )],
                by = c("Sample_B" = "Sample_ID")
            ) |>
            rename(
                Subject_B = "Subject_ID",
                Genotype_Group_B = "Genotype_Group_ID"
            )

        tsmsg(paste(nrow(neighbors), "candidate swaps being evaluated..."))
        all_component_ids <- sort(unique(
            object@.solve_state$unsolved_relabel_data$Component_ID
        ))
        relabels <- data.frame(matrix(
            data = NA,
            nrow = length(all_component_ids),
            ncol = 2,
            dimnames = list(c(), c("relabel_from", "relabel_to"))
        ))
        curr_idx <- 1
        for (curr_component_id in all_component_ids) {
            cc_relabel_data <- unsolved_all_data |>
                filter(.data$Component_ID == curr_component_id)
            cc_neighbors <- neighbors |>
                filter(.data$Component_ID == curr_component_id)

            if (nrow(cc_neighbors) == 0) {
                next
            }
            cc_neighbor_objectives <- cc_neighbors |>
                mutate(
                    delta = mapply(
                        .calc_swapped_delta_entropy,
                        swap_from_subject = .data$Subject_A,
                        swap_from_genotype = .data$Genotype_Group_A,
                        swap_to_subject = .data$Subject_B,
                        swap_to_genotype = .data$Genotype_Group_B,
                        MoreArgs = list(
                            votes = votes,
                            base_entropies = base_entropies
                        )
                    )
                )
            cc_relabels <- cc_neighbor_objectives |>
                filter(.data$delta > 0, .data$delta == max(.data$delta))
            if (nrow(cc_relabels) == 0) {
                next
            }
            cc_relabels <- cc_relabels |>
                sample_n(1) |>
                transmute(
                    relabel_from = .data$Sample_A,
                    relabel_to = .data$Sample_B
                )
            relabels[curr_idx, c("relabel_from", "relabel_to")] <- cc_relabels
            curr_idx <- curr_idx + 1
        }

        relabels <- relabels[!is.na(relabels[, 1]) & !is.na(relabels[, 2]), ]
        relabels <- rbind(
            relabels,
            data.frame(
                relabel_from = relabels$relabel_to,
                relabel_to = relabels$relabel_from
            )
        )
        object <- .relabel_samples(object, relabels, solver_name = "local_old")
        if (nrow(relabels) == 0) {
            # If nothing was relabeled in this iteration, then the algorithm has
            # converged and further iterations are not needed.
            break
        }
    }

    return(object)
}

#' Ensemble Sample Relabeling
#'
#' This ensemble solver uses a combination of global search, majority-search
#' heuristic, and local search to identify and correct mislabels. This allows
#' the ensemble to handle a wider range of mislabel situations than any single
#' algorithm can by itself. This is therefore the recommended way to run FOMO.
#'
#' @param object A MislabelSolver object
#' @param use_solvers A character vector giving the subset of single-method
#'   solvers to run on each iteration of the ensemble loop. Must be a non-empty
#'   subset of `"majority"`, `"global"`, `"local"`. Any solver left out of
#'   `use_solvers` is skipped entirely. For example, setting `use_solvers` to
#'   `c("global", "majority")` will skip local search. Additionally, the old
#'   local search algorithm ([solveLocalSearchOld()]) can be used instead of the
#'   new one ([solveLocalSearch()]) by replacing `"local"` with `"local_old"`.
#' @param time_limit (Default = 7200, i.e. 2 hours) The maximum time, in
#'   seconds, to let the solver run. Elapsed time is checked once per iteration
#'   of the while loop, which means that the actual run time may slightly exceed
#'   the specified limit. If the limit is reached before the solver has
#'   converged, the loop is terminated early, a warning is issued, and the
#'   MislabelSolver object is returned in its current, incompletely solved
#'   state. Note that FOMO will not come anywhere near this time limit for most
#'   real world data sets, so users generally should not need to modify this.
#' @param seed The random seed to use for this run. The seed defaults to 1,
#'   which means this function is reproducible by default. You can disable this
#'   behavior by setting the seed to `NULL`.
#' @param majority_max_genotypes Passed to [solveMajoritySearch()] as its
#'   `max_genotypes` argument (only used if `use_solvers` includes
#'   `"majority"`); see that function's documentation for what it controls.
#' @param global_max_genotypes,global_ghost_penalty,global_deletion_penalty
#'   Passed to [solveGlobalSearch()] as its `max_genotypes`, `ghost_penalty`,
#'   and `deletion_penalty` arguments respectively (only used if `use_solvers`
#'   includes `"global"`); see that function's documentation for what each one
#'   controls.
#' @param local_iter_per_cycle Passed to [solveLocalSearch()] as its `n_iter`
#'   argument (only used if `use_solvers` includes `"local"`); see that
#'   function's documentation for what it controls. Note that, unlike calling
#'   [solveLocalSearch()] directly, this does not control the *total* number of
#'   local search iterations `solveEnsemble()` runs; it controls how many local
#'   search iterations run per cycle, i.e. in between each round of the other
#'   solvers in `use_solvers`. Also unlike calling [solveLocalSearch()]
#'   directly, `solveEnsemble()` always sets its `min_genotypes` argument to
#'   one more than the smallest of `global_max_genotypes`/
#'   `majority_max_genotypes` among whichever of `"global"`/`"majority"` are
#'   also in `use_solvers` (or leaves it at its own no-op default of 1 if
#'   neither is); this is not user-configurable. This keeps local search
#'   focused on components too large for the other active solvers to have
#'   already handled, instead of redoing their work. (This does not apply
#'   when `use_solvers` includes `"local_old"` instead of `"local"`, since
#'   [solveLocalSearchOld()] has no `min_genotypes` argument.)
#'
#' @return A MislabelSolver object with all detected correctable mislabeled
#'   samples relabeled (unless the time limit is reached).
#'
#' @seealso [solveMajoritySearch()], [solveGlobalSearch()],
#'   [solveLocalSearch()], the solvers used in this ensemble algorithm.
#'
#' @export
#'
solveEnsemble <- function(
    object,
    use_solvers = c("majority", "global", "local"),
    time_limit = 2 * 60 * 60,
    seed = 1,
    majority_max_genotypes = 100,
    global_max_genotypes = 8,
    global_ghost_penalty = 1.5,
    global_deletion_penalty = 4,
    local_iter_per_cycle = 1
) {
    object <- fixup_MislabelSolver(object)
    valid_solvers <- c(
        "majority",
        "global",
        "local",
        "local_old"
    )
    assert_that(
        is.character(use_solvers) && length(use_solvers) > 0,
        msg = "'use_solvers' must be a non-empty character vector"
    )
    assert_that(
        all(use_solvers %in% valid_solvers),
        msg = paste0(
            "'use_solvers' must only contain values from: ",
            paste(valid_solvers, collapse = ", ")
        )
    )

    use_solvers <- unique(use_solvers)

    assert_that(
        !all(c("local", "local_old") %in% use_solvers),
        msg = "'use_solvers' cannot contain both 'local' and 'local_old' -- they are alternative implementations of the same step"
    )

    assert_that(
        is.numeric(time_limit) &&
            length(time_limit) == 1 &&
            !is.na(time_limit) &&
            time_limit >= 0,
        msg = "'time_limit' must be a single non-negative number (of seconds)"
    )

    if (!is.null(seed)) {
        assert_that(
            is.numeric(seed) && length(seed) == 1 && !is.na(seed),
            msg = "'seed' must be a single number"
        )
        local_seed(seed)
    }

    ## If local search is the only solver in use, then the outer while loop no
    ## longer alternates between solvers, and its only purpose is to enforce the
    ## time limit. Otherwise, returning control to the outer loop only adds
    ## additional overhead. As such, in this situation we ignore the specfied
    ## local_iter_per_cycle and instead set it to a large value (unless the
    ## caller already specfied an even larger value.)
    only_local_in_use <- all(use_solvers %in% c("local", "local_old"))
    effective_local_iter_per_cycle <- local_iter_per_cycle
    accelerated_local_iter_per_cycle <- 100L
    if (
        only_local_in_use &&
            effective_local_iter_per_cycle < accelerated_local_iter_per_cycle
    ) {
        effective_local_iter_per_cycle <- accelerated_local_iter_per_cycle
        tsmsg(
            "'use_solvers' contains only local search; overriding ",
            "'local_iter_per_cycle' (",
            local_iter_per_cycle,
            ") with ",
            effective_local_iter_per_cycle,
            "."
        )
    }

    run_global <- "global" %in% use_solvers

    ## Tracks the set of samples available for solveGlobalSearch() to analyze
    ## (see .global_search_available_samples() in helpers-solve.R) as of the
    ## last time it actually ran, so it can be skipped below when nothing new
    ## has become available to it since then -- calling it again in that case is
    ## guaranteed to be futile: every component it would look at now is one it
    ## either already fully resolved (and so has left the unsolved pool
    ## entirely) or already skipped for the exact same reason (too large, or
    ## more genotypes than subjects) last time, and neither of those can change
    ## without changing this set. Starts empty, before anything has been
    ## analyzed.
    global_available_samples <- character(0)

    start_time <- Sys.time()
    time_limit_exceeded <- FALSE
    while (TRUE) {
        if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
            tsmsg("All components fully solved. Stopping.")
            break
        }
        if (
            as.numeric(difftime(Sys.time(), start_time, units = "secs")) >
                time_limit
        ) {
            warning(
                "solveEnsemble() reached 'time_limit' of ",
                time_limit,
                " second(s) before converging; returning the object in its current, incompletely solved state."
            )
            ## Out of time; stop looping and return the best-effort partial
            ## result computed so far instead of continuing indefinitely.
            time_limit_exceeded <- TRUE
            break
        }
        prev_relabel_data <- object@.solve_state$unsolved_relabel_data

        if (run_global) {
            current_available_samples <- .global_search_available_samples(
                object,
                global_max_genotypes
            )
            if (
                !identical(current_available_samples, global_available_samples)
            ) {
                object <- solveGlobalSearch(
                    object,
                    max_genotypes = global_max_genotypes,
                    ghost_penalty = global_ghost_penalty,
                    deletion_penalty = global_deletion_penalty
                )
                global_available_samples <- current_available_samples
            } else {
                tsmsg(
                    "Skipping global search: no new samples have become available to it since it last ran."
                )
            }
        }
        if ("majority" %in% use_solvers) {
            object <- solveMajoritySearch(
                object,
                max_genotypes = majority_max_genotypes
            )
        }
        if (run_global) {
            current_available_samples <- .global_search_available_samples(
                object,
                global_max_genotypes
            )
            if (
                !identical(current_available_samples, global_available_samples)
            ) {
                object <- solveGlobalSearch(
                    object,
                    max_genotypes = global_max_genotypes,
                    ghost_penalty = global_ghost_penalty,
                    deletion_penalty = global_deletion_penalty
                )
                global_available_samples <- current_available_samples
            } else {
                tsmsg(
                    "Skipping global search: no new samples have become available to it since it last ran."
                )
            }
        }

        global_relabel_data <- object@.solve_state$unsolved_relabel_data
        if ("local" %in% use_solvers || "local_old" %in% use_solvers) {
            use_local_old <- "local_old" %in% use_solvers
            local_solver <- if (use_local_old) {
                solveLocalSearchOld
            } else {
                solveLocalSearch
            }

            ## solveLocalSearchOld() has no min_genotypes argument, and is
            ## kept as an unmodified reference implementation, so this only
            ## ever applies to solveLocalSearch(). The cap itself is based on
            ## whichever of global/majority search are actually active in
            ## use_solvers: if neither is, there is no other solver for local
            ## search to defer to, so it keeps its own no-op default of 1
            ## (nothing skipped).
            local_solver_extra_args <- list()
            if (!use_local_old) {
                active_max_genotypes <- c(
                    if ("global" %in% use_solvers) global_max_genotypes,
                    if ("majority" %in% use_solvers) majority_max_genotypes
                )
                if (length(active_max_genotypes) > 0) {
                    local_solver_extra_args <- list(
                        min_genotypes = min(active_max_genotypes) + 1
                    )
                }
            }

            object <- do.call(
                local_solver,
                c(
                    list(
                        object,
                        n_iter = effective_local_iter_per_cycle,
                        include_ghost = TRUE,
                        filter_concordant_vertices = TRUE
                    ),
                    local_solver_extra_args
                )
            )

            ## If local search found no swaps, try allowing concordant vertices
            if (
                nrow(global_relabel_data) ==
                    nrow(object@.solve_state$unsolved_relabel_data)
            ) {
                if (
                    identical(
                        global_relabel_data,
                        object@.solve_state$unsolved_relabel_data
                    )
                ) {
                    object <- do.call(
                        local_solver,
                        c(
                            list(
                                object,
                                n_iter = effective_local_iter_per_cycle,
                                include_ghost = TRUE,
                                filter_concordant_vertices = FALSE
                            ),
                            local_solver_extra_args
                        )
                    )
                }
            }
        }

        if (
            nrow(prev_relabel_data) ==
                nrow(object@.solve_state$unsolved_relabel_data)
        ) {
            if (
                identical(
                    prev_relabel_data,
                    object@.solve_state$unsolved_relabel_data
                )
            ) {
                ## This pass made no further progress at all (the unsolved
                ## data is byte-for-byte unchanged from before this
                ## iteration ran), so further iterations would loop forever
                ## without converging any further; stop here.
                tsmsg(
                    "Stopping because the solver has converged. Some samples may remain unsolved."
                )
                break
            }
        }
    }

    # We give the global search one final try (if it is allowed to run at all),
    # since it may have been skipped in the final iteration of the loop. This
    # ensures that the global search always has the last word.
    if (
        run_global &&
            !time_limit_exceeded &&
            # Skip if everything is already solved
            nrow(object@.solve_state$unsolved_relabel_data) > 0
    ) {
        tsmsg("Running final global search")
        object <- solveGlobalSearch(
            object,
            max_genotypes = global_max_genotypes,
            ghost_penalty = global_ghost_penalty,
            deletion_penalty = global_deletion_penalty
        )
        if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
            tsmsg("All components fully solved after final global search.")
        }
    }

    ## After the solve, check if cycles can be broken down

    return(object)
}
