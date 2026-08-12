#' Majority-based Sample Relabeling
#'
#' This heuristic function assigns a Subject_ID to a Genotype_Group_ID if a
#' majority of samples within that Genotype_Group_ID have the same Subject_ID label.
#' In other words, if most samples in a particular group share the same Subject_ID,
#' then that Subject_ID is assigned to the entire group.
#'
#' @param object A MislabelSolver object
#' @param unambiguous_only (Default = FALSE) If true, only correct sample mislabels if they are unambiguous.
#'
#' @return A MislabelSolver object
#'
#'
#' @export
#'
solveMajoritySearch <- function(object, unambiguous_only = FALSE) {
    set.seed(1)
    message("Starting majority search")
    if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
        message("0 samples relabeled")
        return(object)
    }

    ## 1. Update putative subjects
    votes <- table(
        object@.solve_state$unsolved_relabel_data$Genotype_Group_ID,
        object@.solve_state$unsolved_relabel_data$Subject_ID
    )
    votes_by_genotype <- data.frame(
        Genotype_Group_ID = rownames(votes),
        Max_Subject_ID = colnames(votes)[apply(votes, 1, which.max)],
        n = rowSums(votes),
        n_Max_Subject_ID = apply(votes, 1, max)
    ) |>
        filter(
            .data$n_Max_Subject_ID >= 2,
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
            .data$n_Max_Genotype_Group_ID >= 2,
            .data$n_Max_Genotype_Group_ID > n / 2
        ) |>
        rename(Genotype_Group_ID = "Max_Genotype_Group_ID") |>
        select("Subject_ID", "Genotype_Group_ID")
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
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
    putative_subjects <- object@.solve_state$putative_subjects
    relabels <- .find_relabel_cycles_from_putative_subjects(
        unsolved_relabel_data,
        putative_subjects,
        unambiguous_only = unambiguous_only,
        allow_unknowns = FALSE
    )

    ## 3. Relabel samples and update solve state
    object <- .relabel_samples(object, relabels, solver_name = "majority")
    # print(paste(nrow(relabels), "samples relabeled"))
    return(object)
}

#' Local Search Sample Relabeling (Original Algorithm)
#'
#' This search function looks through all possible swaps of 2 samples, and selects
#' the swap that minimizes the sum of within-genotype entropies.
#'
#' This is the original local search algorithm. It is no longer the
#' default -- [solveLocalSearch()] now is -- because its closed-form
#' entropy update is not always bit-exact and can therefore occasionally
#' select a different swap than this function on a near-tie (see
#' [solveLocalSearch()] for the full explanation). Use this version if you
#' need results that are guaranteed identical to this exact algorithm.
#'
#' @param object A MislabelSolver object
#' @param n_iter (Default = 1) The number of
#' @param include_ghost (Default = FALSE) If TRUE, allow swaps to include ghost samples
#' @param filter_concordant_vertices (Default = FALSE) If TRUE, filter out samples
#'                                   with at least one concordant edge
#'
#' @return A MislabelSolver object
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
    set.seed(1)
    message("Starting local search (old)")

    for (i in 1:n_iter) {
        message(paste(
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
            message("0 samples relabeled")
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

        message(paste(nrow(neighbors), "candidate swaps being evaluated..."))
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
        # print(paste(nrow(relabels), "samples relabeled"))
    }

    return(object)
}

#' Global Search Sample Relabeling
#'
#' This global search function permutes over all combinations of assigning a
#' Subject_ID to a Genotype_Group_ID, then picks the assignment that implies
#' the fewest number of sample mislabels and deletions. Scores permutations
#' efficiently: (1) a component's *locked* genotype columns (already
#' resolved before global search reached this component) are constant
#' across every permutation, so their score contribution is computed once
#' from a single row rather than recomputed on every one of up to 8! =
#' 40,320 permutation rows; (2) the per-swap-category scoring itself is
#' done with named-vector lookups and `rowsum()` rather than
#' `left_join()`/`group_by()`/`summarize()`, since profiling showed cost
#' dominated by join/data-mask overhead rather than the arithmetic itself.
#'
#' Use via `solveEnsemble(object, use_solvers = c("majority",
#' "global", "local"))` rather than calling this directly,
#' unless you are composing your own solver loop.
#'
#' @param object A MislabelSolver object
#' @param max_genotypes (Default = 8) The number of combinations scales in factorial
#'                      with the number of genotypes in the largest connected component.
#'                      The algorithm will skip over all components that exceed this size.
#' @param ghost_penalty (Default = 1.5) The score charged, per sample, for relabeling
#'                      to a ghost (placeholder) sample rather than a real one. The
#'                      ordinary relabel penalty is fixed at 1, since all penalties are
#'                      relative to it. Must be a single positive numeric value; a
#'                      warning is issued if it is not strictly greater than 1, since
#'                      otherwise the algorithm may prefer relabeling to a ghost sample
#'                      even when a valid non-ghost sample swap is already available.
#' @param deletion_penalty (Default = 4) The score charged, per sample, for a label or
#'                      genotype deletion (giving up on reconciling that sample
#'                      entirely). Must be a single positive numeric value; a warning
#'                      is issued if it is not strictly greater than *both* twice the
#'                      relabel penalty (i.e. 2) and `ghost_penalty`, since otherwise
#'                      the algorithm may prefer inserting/deleting samples even when a
#'                      valid non-ghost sample swap, or an available ghost sample, is
#'                      already available -- these are two independent comparisons
#'                      (against a plain swap, and against an available ghost) and
#'                      `deletion_penalty` needs to clear both; see the comment above
#'                      `.validate_search_penalties()` in `helpers-solve.R` for the
#'                      worked examples behind both bars, and why neither -- nor their
#'                      combination -- is a hard guarantee for every possible
#'                      component, just the bar cleared by the simplest, most common
#'                      cases. (This default was doubled from 2 to 4 when a
#'                      double-counting bug was fixed: a single orphaned sample -- e.g.
#'                      one displaced by relabeling another sample to a duplicate of it
#'                      -- was previously counted as both a genotype deletion and a
#'                      label deletion and penalized for both, so the default was
#'                      doubled to keep the effective per-sample deletion cost, and
#'                      hence overall solver behavior, unchanged from earlier package
#'                      versions in the common case.)
#'
#' @return A MislabelSolver object
#'
#' @export
#'
solveGlobalSearch <- function(
    object,
    max_genotypes = 8,
    ghost_penalty = 1.5,
    deletion_penalty = 4
) {
    set.seed(1)
    message("Starting global search")
    .validate_search_penalties(ghost_penalty, deletion_penalty)
    if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
        message("0 samples relabeled")
        return(object)
    }

    putative_subjects <- object@.solve_state$putative_subjects

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
        cc_swap_cat_ids <- unique(c(
            cc_unsolved_relabel_data$SwapCat_ID,
            cc_unsolved_ghost_data$SwapCat_ID
        ))
        cc_genotypes <- unique(cc_unsolved_relabel_data$Genotype_Group_ID)
        cc_subjects <- unique(cc_unsolved_relabel_data$Subject_ID)

        if (length(cc_genotypes) > length(cc_subjects)) {
            next
        }

        locked_genotypes <- intersect(
            putative_subjects$Genotype_Group_ID,
            cc_genotypes
        )
        locked_subjects <- intersect(putative_subjects$Subject_ID, cc_subjects)
        free_genotypes <- setdiff(cc_genotypes, locked_genotypes)
        free_subjects <- setdiff(cc_subjects, locked_subjects)

        if (
            length(free_genotypes) > max_genotypes |
                length(free_subjects) > max_genotypes
        ) {
            next
        }

        if (length(free_genotypes) > 0 & length(free_subjects) > 0) {
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

        label_counts <- cc_unsolved_relabel_data |>
            select(
                "Sample_ID",
                "Subject_ID",
                "Genotype_Group_ID",
                "SwapCat_ID"
            ) |>
            group_by(.data$Subject_ID, .data$SwapCat_ID) |>
            summarize(n_labels = n(), .groups = "drop")
        ghost_label_counts <- cc_unsolved_ghost_data |>
            select(
                "Sample_ID",
                "Subject_ID",
                "Genotype_Group_ID",
                "SwapCat_ID"
            ) |>
            group_by(.data$Subject_ID, .data$SwapCat_ID) |>
            summarize(n_ghost_labels = n(), .groups = "drop")
        genotype_counts <- cc_unsolved_relabel_data |>
            select(
                "Sample_ID",
                "Subject_ID",
                "Genotype_Group_ID",
                "SwapCat_ID"
            ) |>
            group_by(.data$Genotype_Group_ID, .data$SwapCat_ID) |>
            summarize(n_in_genotype = n(), .groups = "drop")
        genotype_subject_concordant_counts <- cc_unsolved_relabel_data |>
            select(
                "Sample_ID",
                "Subject_ID",
                "Genotype_Group_ID",
                "SwapCat_ID"
            ) |>
            group_by(
                .data$Subject_ID,
                .data$Genotype_Group_ID,
                .data$SwapCat_ID
            ) |>
            summarize(n_samples_correct = n(), .groups = "drop")

        ## Base-R vectorized locked/free-aware scoring; see .score_permutations_fast() in helpers-solve.R
        permutation_stats <- .score_permutations_fast(
            perm_genotypes,
            free_genotypes,
            locked_genotypes,
            cc_swap_cat_ids,
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

    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
    unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data
    putative_subjects <- object@.solve_state$putative_subjects
    relabels <- .find_relabel_cycles_from_putative_subjects(
        unsolved_relabel_data,
        putative_subjects,
        unsolved_ghost_data,
        allow_unknowns = TRUE
    )

    object <- .relabel_samples(object, relabels, solver_name = "global")
    return(object)
}

#' Local Search Sample Relabeling
#'
#' This search function looks through all possible swaps of 2 samples, and
#' selects the swap that minimizes the sum of within-genotype entropies.
#' Evaluates each candidate swap's entropy delta with a closed-form update
#' instead of rebuilding and rescanning the full vote vector for the
#' affected genotype(s); see [solveLocalSearchOld()], which does the
#' latter and whose cost scales with the width of the votes table -- one
#' column per distinct `Subject_ID` across the *whole* unsolved dataset,
#' not just the candidate swap's own component -- so the speedup from
#' avoiding that rescan grows with dataset size: 9-380x faster than
#' [solveLocalSearchOld()] on synthetic components ranging from 60 to
#' 3,000 candidate subjects in testing.
#'
#' The closed-form update is algebraically exact (an exact rearrangement of
#' [solveLocalSearchOld()]'s entropy formula, true under infinite-precision
#' arithmetic), but it is not a bit-exact re-derivation of it, because it
#' evaluates `log()` on different arguments than the original (see the
#' comment above `.calc_swapped_delta_entropy_fast()` in `helpers-solve.R`
#' for the full explanation). In practice this means entropy deltas from the
#' two functions typically agree to about 1e-9 or tighter, not to the last
#' bit. That is inconsequential on its own, but because [solveLocalSearchOld()]
#' selects a single best swap per component via `delta == max(delta)`, a
#' difference this small can occasionally change which swap wins a near-tie
#' -- so on rare inputs, this function's swap selections (and hence its
#' output) can differ slightly from [solveLocalSearchOld()]'s, even though
#' both are picking from among equally-good (or all but indistinguishably
#' good) options. If you need results that exactly match
#' [solveLocalSearchOld()], use that function instead.
#'
#' Use via `solveEnsemble(object, use_solvers = c("majority", "global",
#' "local"))` rather than calling this directly, unless you are
#' composing your own solver loop.
#'
#' @inheritParams solveLocalSearchOld
#'
#' @return A MislabelSolver object
#'
#' @seealso [solveLocalSearchOld()], the original (exact) version of this
#'   algorithm, retained for cases needing bit-exact reproducibility.
#'
#' @export
#'
solveLocalSearch <- function(
    object,
    n_iter = 1,
    include_ghost = FALSE,
    filter_concordant_vertices = FALSE
) {
    set.seed(1)
    message("Starting local search")

    for (i in 1:n_iter) {
        message(paste(
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
            message("0 samples relabeled")
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

        message(paste(nrow(neighbors), "candidate swaps being evaluated..."))
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
        # print(paste(nrow(relabels), "samples relabeled"))
    }

    return(object)
}

#' Ensemble Sample Relabeling
#'
#' This ensemble solver uses a combination of majority-search heuristic,
#' global search, and local search to identify and correct mislabels
#'
#' @param object A MislabelSolver object
#' @param use_solvers (Default = `c("majority", "global", "local")`) A
#'   character vector giving the subset of single-method solvers to run on
#'   each iteration of the ensemble loop. Must be a non-empty subset of
#'   `"majority"`, `"global"`, `"local"`, and `"local_old"`. Any solver left
#'   out of `use_solvers` is skipped entirely. For example, setting
#'   `use_solvers` to `c("global", "majority")` will skip local search.
#'   `"local"` is [solveLocalSearch()], the standard local search algorithm;
#'   `"local_old"` is [solveLocalSearchOld()], the original algorithm it
#'   replaced as the default. The two are alternative implementations of
#'   the same step and cannot both be requested at once: [solveLocalSearch()]
#'   is only extremely close (not bit-exact) to [solveLocalSearchOld()] --
#'   see [solveLocalSearch()] for why that distinction exists -- so use
#'   `"local_old"` if you need results that exactly match it.
#' @param time_limit (Default = 7200, i.e. 2 hours) The maximum time, in
#'   seconds, to let the solver run. Elapsed time is checked once per
#'   iteration of the while loop; if `time_limit` is reached before the
#'   solver has converged, the loop is stopped early, a warning is issued,
#'   and the object is returned in its current, incompletely solved state.
#' @param seed (Default = 1) The random seed, passed to `set.seed()`, used
#'   for reproducibility.
#' @param global_max_genotypes,global_ghost_penalty,global_deletion_penalty
#'   Passed to [solveGlobalSearch()] as its `max_genotypes`, `ghost_penalty`,
#'   and `deletion_penalty` arguments respectively (only relevant when
#'   `"global"` is in `use_solvers`); see that function's documentation for
#'   what each one controls.
#' @param local_iter_per_cycle (Default = 1) Passed to whichever of
#'   [solveLocalSearch()]/[solveLocalSearchOld()] is in use (per
#'   `use_solvers`) as its `n_iter` argument (only relevant when `"local"`
#'   or `"local_old"` is in `use_solvers`); see that function's
#'   documentation for what it controls. Note that, unlike calling
#'   [solveLocalSearch()] directly, this does not control the *total*
#'   number of local search iterations solveEnsemble() runs -- it controls
#'   how many local search iterations run per cycle, i.e. in between each
#'   round of the other solvers in `use_solvers`. Ignored (with a message)
#'   if `use_solvers` contains only `"local"`/`"local_old"` and neither
#'   `"majority"` nor `"global"`, since in that case local search is the
#'   only thing the ensemble loop does anyway; a large internal default is
#'   used for `n_iter` instead so the ensemble loop still periodically
#'   returns control to itself (to check `time_limit`, etc.) without the
#'   overhead of doing so after every single local search iteration.
#'
#' @return A MislabelSolver object
#'
#' @export
#'
solveEnsemble <- function(
    object,
    use_solvers = c("majority", "global", "local"),
    time_limit = 2 * 60 * 60,
    seed = 1,
    global_max_genotypes = 8,
    global_ghost_penalty = 1.5,
    global_deletion_penalty = 4,
    local_iter_per_cycle = 1
) {
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

    assert_that(
        is.numeric(seed) && length(seed) == 1 && !is.na(seed),
        msg = "'seed' must be a single number"
    )

    set.seed(seed)

    ## If local search is the only solver in use (no "majority" and no
    ## "global"), nothing else in the outer while loop below does anything
    ## on any iteration, so there's no benefit to returning control to it
    ## after only local_iter_per_cycle local search iterations -- doing so
    ## only adds the outer loop's own per-iteration overhead (the
    ## identical() comparisons below, mainly) with nothing to show for it.
    ## Give local search a much larger iteration budget instead in that
    ## case, regardless of what local_iter_per_cycle was set to.
    ## Deliberately a large, finite value rather than Inf or
    ## .Machine$integer.max: the outer loop is still what enforces
    ## time_limit and emits its own progress messages, and local search
    ## isn't proven to always converge within a bounded number of
    ## iterations, so control still needs to come back to the outer loop
    ## periodically -- aiming for roughly once a minute -- rather than
    ## potentially running unchecked for a very long time.
    only_local_in_use <- all(use_solvers %in% c("local", "local_old"))
    effective_local_iter_per_cycle <- local_iter_per_cycle
    if (only_local_in_use) {
        effective_local_iter_per_cycle <- 1000L
        message(
            "'use_solvers' contains only local search; overriding ",
            "'local_iter_per_cycle' (",
            local_iter_per_cycle,
            ") with ",
            effective_local_iter_per_cycle,
            " for this solveEnsemble() call."
        )
    }

    run_global <- "global" %in% use_solvers

    ## Tracks the set of samples available for solveGlobalSearch() to
    ## analyze (see .global_search_available_samples() in helpers-solve.R)
    ## as of the last time it actually ran, so it can be skipped below when
    ## nothing new has become available to it since then -- calling it
    ## again in that case is guaranteed to be futile: every component it
    ## would look at now is one it either already fully resolved (and so
    ## has left the unsolved pool entirely) or already skipped for the
    ## exact same reason (too large, or more genotypes than subjects) last
    ## time, and neither of those can change without changing this set.
    ## Starts empty, before anything has been analyzed.
    global_available_samples <- character(0)

    start_time <- Sys.time()
    time_limit_exceeded <- FALSE
    while (TRUE) {
        if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
            ## Everything has already been resolved; nothing left to solve.
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
                message(
                    "Skipping global search: no new samples have become ",
                    "available to it since it last ran."
                )
            }
        }
        if ("majority" %in% use_solvers) {
            object <- solveMajoritySearch(object)
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
                message(
                    "Skipping global search: no new samples have become ",
                    "available to it since it last ran."
                )
            }
        }

        global_relabel_data <- object@.solve_state$unsolved_relabel_data
        if ("local" %in% use_solvers || "local_old" %in% use_solvers) {
            local_solver <- if ("local_old" %in% use_solvers) {
                solveLocalSearchOld
            } else {
                solveLocalSearch
            }
            object <- local_solver(
                object,
                n_iter = effective_local_iter_per_cycle,
                include_ghost = TRUE,
                filter_concordant_vertices = TRUE
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
                    object <- local_solver(
                        object,
                        n_iter = effective_local_iter_per_cycle,
                        include_ghost = TRUE,
                        filter_concordant_vertices = FALSE
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
                break
            }
        }
    }

    ## The loop above can skip calling global search on its very last
    ## iteration (if nothing new became available to it right at the end),
    ## so run it once more here to make sure it never ends up skipped
    ## entirely -- unless there's no point: either the loop stopped because
    ## it ran out of time (in which case any further solving, global search
    ## included, should also be skipped), or the caller didn't request
    ## global search in the first place.
    if (run_global && !time_limit_exceeded) {
        object <- solveGlobalSearch(
            object,
            max_genotypes = global_max_genotypes,
            ghost_penalty = global_ghost_penalty,
            deletion_penalty = global_deletion_penalty
        )
    }

    ## After the solve, check if cycles can be broken down

    return(object)
}
