#' Produce tables summarizing the output of FOMO
#'
#' @param object An object of class \code{MislabelSolver}.
#'
#' @returns A named list of four data frames summarizing the results at
#'   different levels of granularity -- "Sample", "Genotype_Group",
#'   "Component", and "Dataset" -- each a coarser aggregation of the one
#'   before it. The "Sample" sheet is the primary, most granular table (one
#'   row per sample) and includes:
#'   \itemize{
#'     \item \code{Solved_By}: which solver step resolved this sample's
#'       *component* (\code{"majority"}, \code{"global"},
#'       \code{"local"}/\code{"local_old"}, or \code{"initial"} if it was
#'       already resolved at construction time), \code{NA} if still
#'       unsolved. Set once, permanently, the moment the component is solved
#'       -- it does not mean this specific sample was itself relabeled.
#'     \item \code{Relabeled_By}: which solver most recently changed this
#'       specific sample's own label, updated every time a relabel happens
#'       (\code{NA} if the sample was never relabeled). Can differ from
#'       \code{Solved_By} -- e.g. a sample can be relabeled by
#'       \code{"local"} and have its component later marked solved by a
#'       \code{"global"} pass that didn't need to touch that sample again.
#'     \item \code{Selected_For_Review}: a coarse triage category for
#'       manual review, one of \code{"ghost_relabeled"}, \code{"ghost"},
#'       \code{"inconsistent_genotype"}, \code{"deletion_or_duplication"},
#'       \code{"relabel_low_confidence"}, \code{"relabel_high_confidence"},
#'       \code{"singleton_no_inference"},
#'       \code{"not_relabeled_low_confidence"}, or
#'       \code{"no_review_needed"}.
#'     \item Contamination metrics (\code{Sample_Contamination_Metric} and
#'       its numerator/denominator), estimating, from \code{genotype_matrix}
#'       (when provided), what fraction of a sample's expected genetic
#'       neighbors don't actually show up as genetically related.
#'     \item \code{Mislabeling_Event_ID}: samples whose corrections are
#'       linked (e.g. a swap between two samples) share the same ID.
#'   }
#'   The "Genotype_Group", "Component", and "Dataset" sheets aggregate these
#'   same fields (counts of each \code{Selected_For_Review} category, etc.)
#'   at increasingly coarse levels; see the column names themselves (most
#'   are self-descriptive, e.g. \code{n_Samples_relabel_high_confidence}).
#' @export
#'
#' @seealso [writeOutput()]
collateOutput <- function(object) {
    sample_summary <- object@.solve_state$relabel_data |>
        group_by(.data$Genotype_Group_ID, .data$Init_Subject_ID) |>
        mutate(n_agree = n()) |>
        ungroup("Init_Subject_ID") |>
        mutate(Sample_Count_In_Genotype_Group = n()) |>
        ungroup("Genotype_Group_ID") |>
        left_join(
            object@.solve_state$putative_subjects |>
                filter(!is.na(.data$Genotype_Group_ID)),
            by = "Genotype_Group_ID",
            suffix = c("", "_putative")
        ) |>
        transmute(
            Component_ID = .data$Init_Component_ID,
            .data$Genotype_Group_ID,
            .data$SwapCat_ID,
            Ghost = is.na(.data$Genotype_Group_ID),
            .data$Init_Subject_ID,
            .data$Init_Sample_ID,
            Proposed_Final_Subject_ID = .data$Subject_ID,
            Proposed_Final_Sample_ID = .data$Sample_ID,
            Inferred_Subject_ID = .data$Subject_ID_putative,
            All_Valid_Subject_IDs = sapply(
                .data$Genotype_Group_ID,
                function(x) {
                    if (x %in% names(object@.solve_state$ambiguous_subjects)) {
                        str_c(
                            object@.solve_state$ambiguous_subjects[[x]],
                            collapse = ", "
                        )
                    } else {
                        NA_character_
                    }
                }
            ),
            Sample_Count_In_Genotype_Group = if_else(
                !.data$Ghost,
                .data$Sample_Count_In_Genotype_Group,
                NA_integer_
            ),
            Sample_Count_In_Genotype_Group_with_Same_Initial_Subject_Label = if_else(
                !.data$Ghost,
                .data$n_agree,
                NA_integer_
            ),
            Mislabeled = .data$Init_Sample_ID != .data$Proposed_Final_Sample_ID,
            ## Which solver step resolved this sample: "majority", "global",
            ## "local"/"local_old", or "initial" if it was already resolved at
            ## construction time (e.g. a singleton component, or
            ## anchor_samples), before any solver ran. NA if still unsolved
            ## (only possible if solveEnsemble() hit its time_limit before
            ## converging). TODO I think this comment may be inaccurate.
            .data$Solved_By,
            ## Which solver most recently changed this sample's label, updated
            ## every time a relabel happens (unlike Solved_By, which is set
            ## once, permanently, the moment a component is solved). NA if the
            ## sample was never relabeled. Can still be non-NA even when
            ## Init_Sample_ID == Proposed_Final_Sample_ID, in the rare case
            ## where a later solver reverts an earlier relabel back to the
            ## original label -- that's expected, not a bug.
            .data$Relabeled_By,
            Selected_For_Review = case_when(
                .data$Ghost & .data$Mislabeled ~ "ghost_relabeled",
                .data$Ghost ~ "ghost",
                is.na(.data$Inferred_Subject_ID) ~ "inconsistent_genotype",
                str_detect(
                    .data$Proposed_Final_Sample_ID,
                    LABEL_NOT_FOUND
                ) ~ "deletion_or_duplication",
                .data$Mislabeled &
                    (.data$Proposed_Final_Subject_ID !=
                        .data$Inferred_Subject_ID |
                        .data$Sample_Count_In_Genotype_Group ==
                            1) ~ "relabel_low_confidence",
                .data$Mislabeled ~ "relabel_high_confidence",
                .data$Sample_Count_In_Genotype_Group ==
                    1 ~ "singleton_no_inference",
                .data$n_agree < 2 ~ "not_relabeled_low_confidence",
                TRUE ~ "no_review_needed"
            )
        ) |>
        arrange(
            .data$Component_ID,
            .data$Genotype_Group_ID,
            .data$Proposed_Final_Subject_ID,
            .data$Proposed_Final_Sample_ID
        )

    ## Mislabeled samples in the same genotype group
    ## with identical swappable categories are
    ## are ambiguities for one another
    ambiguity_summary <- sample_summary |>
        filter(.data$Init_Subject_ID != .data$Inferred_Subject_ID) |>
        group_by(.data$Genotype_Group_ID, .data$SwapCat_ID) |>
        mutate(
            n_LABELNOTFOUND = sum(str_detect(
                .data$Proposed_Final_Sample_ID,
                LABEL_NOT_FOUND
            )),
            has_Ghost_Solution = any(.data$Ghost),
            has_LABELNOTFOUND_Solution = .data$n_LABELNOTFOUND > 0
        ) |>
        ungroup()
    ambiguity_summary$All_Valid_Sample_IDs <- NA_character_
    for (i in seq_len(nrow(ambiguity_summary))) {
        genotype_group_id <- ambiguity_summary$Genotype_Group_ID[i]
        swap_cat_id <- ambiguity_summary$SwapCat_ID[i]
        inferred_subject_id <- ambiguity_summary$Inferred_Subject_ID[i]
        has_LABELNOTFOUND <- ambiguity_summary$has_LABELNOTFOUND_Solution[i]
        sample_ambiguities <- ambiguity_summary |>
            filter(
                .data$Genotype_Group_ID != genotype_group_id,
                .data$Init_Subject_ID == inferred_subject_id,
                .data$SwapCat_ID == swap_cat_id
            ) |>
            pull(.data$Init_Sample_ID)
        if (has_LABELNOTFOUND) {
            sample_ambiguities <- c(
                sample_ambiguities,
                str_c(
                    LABEL_NOT_FOUND,
                    inferred_subject_id,
                    swap_cat_id,
                    collapse = "#"
                )
            )
        }
        if (length(sample_ambiguities) > 1) {
            ambiguity_summary[i, "All_Valid_Sample_IDs"] <- str_c(
                sample_ambiguities,
                collapse = ", "
            )
        }
    }
    ambiguity_summary <- ambiguity_summary |>
        select("Init_Sample_ID", "All_Valid_Sample_IDs")

    sample_summary <- sample_summary |>
        left_join(ambiguity_summary, by = "Init_Sample_ID") |>
        mutate(
            Multiple_Valid_Solutions = case_when(
                !is.na(.data$All_Valid_Subject_IDs) >
                    0 ~ "multiple_valid_subjects",
                !is.na(.data$All_Valid_Sample_IDs) >
                    0 ~ "one_valid_subject_multiple_valid_samples"
            )
        )

    sample_summary$Sample_Contamination_Metric_Numerator <- NA_integer_
    sample_summary$Sample_Contamination_Metric_Denominator <- NA_integer_
    sample_summary$Sample_Contamination_Metric <- NA
    if (!is.null(object@genotype_matrix)) {
        sample_summary$Sample_Contamination_Metric_Numerator <- NA_integer_
        sample_summary$Sample_Contamination_Metric_Denominator <- NA_integer_
        sample_summary$Sample_Contamination_Metric <- NA
        genotyped_sample_ids <- sample_summary |>
            filter(!.data$Ghost) |>
            pull(.data$Init_Sample_ID)
        for (sample_id in genotyped_sample_ids) {
            neighbor_samples <- names(which(
                object@genotype_matrix[sample_id, ] == 1
            ))
            ## drop = FALSE: if sample_id has exactly one neighbor, a plain
            ## `[` here would silently collapse this to a bare scalar (both
            ## dimensions have extent 1); upper.tri() below happens to treat
            ## that scalar like a 1x1 matrix and still return the right
            ## (empty) answer, but only by coincidence. Forcing a real 1x1
            ## matrix here makes that explicit and correct by construction
            ## rather than by accident, with no change in the actual numbers
            ## (verified for 0, 1, and 2+ neighbors).
            neighbor_matrix <- object@genotype_matrix[
                neighbor_samples,
                neighbor_samples,
                drop = FALSE
            ]
            existing_edges <- sum(neighbor_matrix[upper.tri(neighbor_matrix)])
            total_edges <- length(neighbor_matrix[upper.tri(neighbor_matrix)])
            missing_edges <- total_edges - existing_edges
            sample_summary[
                sample_summary$Init_Sample_ID == sample_id,
                "Sample_Contamination_Metric_Numerator"
            ] <- missing_edges
            sample_summary[
                sample_summary$Init_Sample_ID == sample_id,
                "Sample_Contamination_Metric_Denominator"
            ] <- total_edges
            if (total_edges > 0) {
                sample_summary[
                    sample_summary$Init_Sample_ID == sample_id,
                    "Sample_Contamination_Metric"
                ] <- missing_edges / total_edges
            }
        }
    }

    corrections_graph <- .generate_corrections_graph(
        object@.solve_state$relabel_data
    )
    corrections_components <- components(corrections_graph)
    Mislabeling_Event_ID_df <- data.frame
    Mislabeling_Event_ID_df <- data.frame(
        Init_Sample_ID = names(corrections_components$membership),
        Mislabeling_Event_ID = corrections_components$membership
    )
    sample_summary <- sample_summary |>
        left_join(Mislabeling_Event_ID_df, by = "Init_Sample_ID")
    Mislabeling_Event_ID_renamer_df <- sample_summary |>
        filter(!is.na(.data$Mislabeling_Event_ID)) |>
        select("Component_ID", "Mislabeling_Event_ID") |>
        distinct() |>
        group_by(.data$Component_ID) |>
        mutate(Event_Number = row_number()) |>
        ungroup() |>
        mutate(
            Mislabeling_Event_ID_New = str_c(
                .data$Component_ID,
                "_Mislabeling_Event_",
                .data$Event_Number
            )
        ) |>
        select("Mislabeling_Event_ID", "Mislabeling_Event_ID_New")
    Mislabeling_Event_ID_renamer_map <- setNames(
        Mislabeling_Event_ID_renamer_df$Mislabeling_Event_ID_New,
        Mislabeling_Event_ID_renamer_df$Mislabeling_Event_ID
    )
    ## Mislabeling_Event_ID_renamer_map's values are character (built via
    ## str_c() below), so the "missing" branch needs a typed NA_character_
    ## rather than a bare (logical) NA for if_else()'s stricter type check.
    sample_summary$Mislabeling_Event_ID <- sapply(
        sample_summary$Mislabeling_Event_ID,
        \(x) {
            if_else(
                is.na(x),
                NA_character_,
                Mislabeling_Event_ID_renamer_map[x]
            )
        }
    )

    sample_summary <- sample_summary |>
        select(
            Connected_Component_ID = "Component_ID",
            "Genotype_Group_ID",
            "SwapCat_ID",
            Is_Ghost = "Ghost",
            Initial_Subject_ID = "Init_Subject_ID",
            Initial_Sample_ID = "Init_Sample_ID",
            "Selected_For_Review",
            "Mislabeled",
            "Solved_By",
            "Relabeled_By",
            "Mislabeling_Event_ID",
            "Multiple_Valid_Solutions",
            "All_Valid_Subject_IDs",
            "All_Valid_Sample_IDs",
            "Inferred_Subject_ID",
            "Proposed_Final_Subject_ID",
            "Proposed_Final_Sample_ID",
            "Sample_Count_In_Genotype_Group",
            "Sample_Count_In_Genotype_Group_with_Same_Initial_Subject_Label",
            "Sample_Contamination_Metric",
            "Sample_Contamination_Metric_Denominator",
            "Sample_Contamination_Metric_Numerator"
        )

    genotype_group_summary <- sample_summary |>
        group_by(.data$Genotype_Group_ID, .data$Inferred_Subject_ID) |>
        filter(!is.na(.data$Genotype_Group_ID)) |>
        summarize(
            # Majority_Subject_ID = names(sort(table(Proposed_Final_Subject_ID), decreasing = TRUE)[1]),
            n_Samples_no_review_needed = sum(
                .data$Selected_For_Review == "no_review_needed"
            ),
            n_Samples_inconsistent_genotype = sum(
                .data$Selected_For_Review == "inconsistent_genotype"
            ),
            n_Samples_deletion_or_duplication = sum(
                .data$Selected_For_Review == "deletion_or_duplication"
            ),
            n_Samples_relabel_low_confidence = sum(
                .data$Selected_For_Review == "relabel_low_confidence"
            ),
            n_Samples_relabel_high_confidence = sum(
                .data$Selected_For_Review == "relabel_high_confidence"
            ),
            n_Samples_singleton_no_inference = sum(
                .data$Selected_For_Review == "singleton_no_inference"
            ),
            n_Samples_not_relabeled_low_confidence = sum(
                .data$Selected_For_Review == "not_relabeled_low_confidence"
            ),
            n_Samples_total = n(),
            n_Samples_Initially_Matching_Inferred_Subject = sum(
                .data$Initial_Subject_ID == .data$Inferred_Subject_ID
            ),
            Selected_For_Review = case_when(
                .data$n_Samples_total ==
                    .data$n_Samples_no_review_needed ~ "no_review_needed",
                .data$n_Samples_total ==
                    .data$n_Samples_singleton_no_inference ~ "singleton_no_inference",
                TRUE ~ "check_sample_table"
            )
        ) |>
        ungroup() |>
        select(
            "Genotype_Group_ID",
            "n_Samples_total",
            "Inferred_Subject_ID",
            "Selected_For_Review",
            "n_Samples_Initially_Matching_Inferred_Subject",
            everything()
        )

    genotype_group_summary$Genotype_Contamination_Metric <- NA
    genotype_group_summary$Genotype_Contamination_Metric_Denominator <- NA_integer_
    genotype_group_summary$Genotype_Contamination_Metric_Numerator <- NA_integer_
    if (!is.null(object@genotype_matrix)) {
        genotype_group_summary$Genotype_Contamination_Metric <- NA
        genotype_group_summary$Genotype_Contamination_Metric_Denominator <- NA_integer_
        genotype_group_summary$Genotype_Contamination_Metric_Numerator <- NA_integer_
        for (genotype_group_id in genotype_group_summary$Genotype_Group_ID) {
            genotype_group_init_samples <- sample_summary |>
                filter(.data$Genotype_Group_ID == genotype_group_id) |>
                pull(.data$Initial_Sample_ID)
            ## drop = FALSE: see the identical reasoning in the per-sample
            ## contamination metric loop above (this is the same pattern, one
            ## element short of matrix-ness when a genotype group has exactly
            ## one member).
            genotype_group_matrix <- object@genotype_matrix[
                genotype_group_init_samples,
                genotype_group_init_samples,
                drop = FALSE
            ]
            existing_edges <- sum(genotype_group_matrix[upper.tri(
                genotype_group_matrix
            )])
            total_edges <- length(genotype_group_matrix[upper.tri(
                genotype_group_matrix
            )])
            missing_edges <- total_edges - existing_edges
            genotype_group_summary[
                genotype_group_summary$Genotype_Group_ID == genotype_group_id,
                "Genotype_Contamination_Metric_Numerator"
            ] <- missing_edges
            genotype_group_summary[
                genotype_group_summary$Genotype_Group_ID == genotype_group_id,
                "Genotype_Contamination_Metric_Denominator"
            ] <- total_edges
            if (total_edges != 0) {
                genotype_group_summary[
                    genotype_group_summary$Genotype_Group_ID ==
                        genotype_group_id,
                    "Genotype_Contamination_Metric"
                ] <- missing_edges / total_edges
            }
        }
    }

    component_summary <- sample_summary |>
        group_by(.data$Connected_Component_ID) |>
        summarize(
            n_Genotype_Groups = length(unique(.data$Genotype_Group_ID)),
            n_Subjects = length(unique(.data$Initial_Subject_ID)),
            n_Samples_total = n(),
            Same_Number_of_Genotypes_And_Subjects = .data$n_Genotype_Groups ==
                .data$n_Subjects,
            n_Samples_no_review_needed = sum(
                .data$Selected_For_Review == "no_review_needed"
            ),
            n_Samples_ghost_relabeled = sum(
                .data$Selected_For_Review == "ghost_relabeled"
            ),
            n_Samples_ghost = sum(.data$Selected_For_Review == "ghost"),
            n_Samples_inconsistent_genotype = sum(
                .data$Selected_For_Review == "inconsistent_genotype"
            ),
            n_Samples_deletion_or_duplication = sum(
                .data$Selected_For_Review == "deletion_or_duplication"
            ),
            n_Samples_relabel_low_confidence = sum(
                .data$Selected_For_Review == "relabel_low_confidence"
            ),
            n_Samples_relabel_high_confidence = sum(
                .data$Selected_For_Review == "relabel_high_confidence"
            ),
            n_Samples_singleton_no_inference = sum(
                .data$Selected_For_Review == "singleton_no_inference"
            ),
            n_Samples_not_relabeled_low_confidence = sum(
                .data$Selected_For_Review == "not_relabeled_low_confidence"
            ),
            ## Which solver(s) contributed to resolving this component's
            ## samples -- e.g. "global" if a single solver step resolved
            ## everything in it, "local, majority" if resolving it took a
            ## combination. Excludes any still-unsolved samples (Solved_By
            ## NA); those are already visible via the n_Samples_* columns
            ## above and Same_Number_of_Genotypes_And_Subjects.
            Solved_By = str_c(
                sort(unique(.data$Solved_By[!is.na(.data$Solved_By)])),
                collapse = ", "
            ),
            Selected_For_Review = case_when(
                .data$n_Samples_total ==
                    .data$n_Samples_no_review_needed ~ "no_review_needed",
                .data$n_Samples_total ==
                    .data$n_Samples_singleton_no_inference ~ "singleton_no_inference",
                TRUE ~ "check_sample_table"
            )
        )

    dataset_summary <- sample_summary |>
        summarize(
            n_Components = length(unique(.data$Connected_Component_ID)),
            n_Genotype_Groups = length(unique(.data$Genotype_Group_ID)),
            n_Subjects = length(unique(.data$Initial_Subject_ID)),
            n_Samples_total = n(),
            n_Samples_no_review_needed = sum(
                .data$Selected_For_Review == "no_review_needed"
            ),
            n_Samples_ghost_relabeled = sum(
                .data$Selected_For_Review == "ghost_relabeled"
            ),
            n_Samples_ghost = sum(.data$Selected_For_Review == "ghost"),
            n_Samples_inconsistent_genotype = sum(
                .data$Selected_For_Review == "inconsistent_genotype"
            ),
            n_Samples_deletion_or_duplication = sum(
                .data$Selected_For_Review == "deletion_or_duplication"
            ),
            n_Samples_relabel_low_confidence = sum(
                .data$Selected_For_Review == "relabel_low_confidence"
            ),
            n_Samples_relabel_high_confidence = sum(
                .data$Selected_For_Review == "relabel_high_confidence"
            ),
            n_Samples_singleton_no_inference = sum(
                .data$Selected_For_Review == "singleton_no_inference"
            ),
            n_Samples_not_relabeled_low_confidence = sum(
                .data$Selected_For_Review == "not_relabeled_low_confidence"
            )
        )

    summary_list <- list(
        "Sample" = sample_summary,
        "Genotype_Group" = genotype_group_summary,
        "Component" = component_summary,
        "Dataset" = dataset_summary
    )

    return(summary_list)
}


#' Write the output of FOMO to an Excel file.
#'
#' @param object Either a MislabelSolver object, or a named list of data frames
#'   (i.e. the return value of [collateOutput()]).
#' @param file The file name to write to. This should end in ".xlsx".
#'
#' @export
#'
#' @seealso [collateOutput()]
writeOutput <- function(object, file) {
    if (is(object, "MislabelSolver")) {
        object <- collateOutput(object)
    }
    assert_that(
        is.list(object),
        all(vapply(object, is, logical(1), "data.frame")),
        msg = "object must be either a MislabelSolver or the result of collateOutput()."
    )
    # Freeze the header row
    formatted_tables <- lapply(object, xl_sheet, freeze = "A2")
    write_xlsx(formatted_tables, path = file)
    message(paste0("Output successfully written to ", file))
    invisible(NULL)
}
