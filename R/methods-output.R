#' Plot corrections for \code{MislabelSolver} objects
#'
#' Leverages \code{visNetwork} to generate an interactive plot of sample
#' relabels
#'
#' @param object An object of class \code{MislabelSolver}
#' @param dir_path Output directory where solver results are written. Will be
#'   created if directory doesn't already exist. If no directory provided,
#'   output will be written to the current working directory.
#'
#' @export
#'
writeOutput <- function(object, dir_path=NULL) {
    if (is.null(dir_path)) {
        dir_path <- getwd()
    }

    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive=FALSE)
        if (dir.exists(dir_path)) {
            message(paste0("Created new directory ", dir_path))
        } else {
            message(paste0("Failed to create new directory ", dir_path, ". Please check if parent directory for the provided path exists."))
        }
    }

    sample_summary <- object@.solve_state$relabel_data |>
        dplyr::group_by(Genotype_Group_ID, Init_Subject_ID) |>
        dplyr::mutate(n_agree=n()) |>
        dplyr::ungroup(Init_Subject_ID) |>
        dplyr::mutate(Sample_Count_In_Genotype_Group=n()) |>
        dplyr::ungroup(Genotype_Group_ID) |>
        dplyr::left_join(
            object@.solve_state$putative_subjects |>
                filter(!is.na(Genotype_Group_ID)),
            by = "Genotype_Group_ID",
            suffix = c("", "_putative")
        ) |>
        dplyr::transmute(
            Component_ID = Init_Component_ID,
            Genotype_Group_ID,
            SwapCat_ID,
            Ghost = is.na(Genotype_Group_ID),
            Init_Subject_ID,
            Init_Sample_ID,
            Proposed_Final_Subject_ID = Subject_ID,
            Proposed_Final_Sample_ID = Sample_ID,
            Inferred_Subject_ID = Subject_ID_putative,
            All_Valid_Subject_IDs = sapply(Genotype_Group_ID, function(x) {
                if (x %in% names(object@.solve_state$ambiguous_subjects)) {
                    paste(object@.solve_state$ambiguous_subjects[[x]], collapse=", ")
                } else {
                    NA_character_
                }
            }),
            Sample_Count_In_Genotype_Group = ifelse(!Ghost, Sample_Count_In_Genotype_Group, NA_integer_),
            Sample_Count_In_Genotype_Group_with_Same_Initial_Subject_Label = ifelse(!Ghost, n_agree, NA_integer_),
            Mislabeled = Init_Sample_ID != Proposed_Final_Sample_ID,
            Selected_For_Review = case_when(
                Ghost & Mislabeled ~ "ghost_relabeled",
                Ghost ~ "ghost",
                is.na(Inferred_Subject_ID) ~ "inconsistent_genotype",
                grepl(LABEL_NOT_FOUND, Proposed_Final_Sample_ID) ~ "deletion_or_duplication",
                Mislabeled & (Proposed_Final_Subject_ID != Inferred_Subject_ID | Sample_Count_In_Genotype_Group == 1) ~ "relabel_low_confidence",
                Mislabeled ~ "relabel_high_confidence",
                Sample_Count_In_Genotype_Group == 1 ~ "singleton_no_inference",
                n_agree < 2 ~ "not_relabeled_low_confidence",
                TRUE ~ "no_review_needed"
            )
        ) |>
        dplyr::arrange(Component_ID, Genotype_Group_ID, Proposed_Final_Subject_ID, Proposed_Final_Sample_ID)

    ## Mislabeled samples in the same genotype group
    ## with identical swappable categories are
    ## are ambiguities for one another
    ambiguity_summary <- sample_summary |>
        dplyr::filter(Init_Subject_ID != Inferred_Subject_ID) |>
        dplyr::group_by(Genotype_Group_ID, SwapCat_ID) |>
        dplyr::mutate(
            n_LABELNOTFOUND = sum(grepl(LABEL_NOT_FOUND, Proposed_Final_Sample_ID)),
            has_Ghost_Solution = any(Ghost),
            has_LABELNOTFOUND_Solution = .data$n_LABELNOTFOUND > 0
        ) |>
        dplyr::ungroup()
    ambiguity_summary$All_Valid_Sample_IDs <- NA_character_
    for (i in seq_len(nrow(ambiguity_summary))) {
        genotype_group_id <- ambiguity_summary$Genotype_Group_ID[i]
        swap_cat_id <- ambiguity_summary$SwapCat_ID[i]
        inferred_subject_id <- ambiguity_summary$Inferred_Subject_ID[i]
        has_LABELNOTFOUND <- ambiguity_summary$has_LABELNOTFOUND_Solution[i]
        sample_ambiguities <- ambiguity_summary |>
            filter(
                Genotype_Group_ID != genotype_group_id,
                Init_Subject_ID == inferred_subject_id,
                SwapCat_ID == swap_cat_id
            ) |>
            pull(Init_Sample_ID)
        if (has_LABELNOTFOUND) {
            sample_ambiguities <- c(sample_ambiguities, paste0(LABEL_NOT_FOUND, inferred_subject_id, swap_cat_id, collapse="#"))
        }
        if (length(sample_ambiguities) > 1) {
            ambiguity_summary[i, "All_Valid_Sample_IDs"] <- paste(sample_ambiguities, collapse=", ")
        }
    }
    ambiguity_summary <- ambiguity_summary |>
        dplyr::select(Init_Sample_ID, All_Valid_Sample_IDs)

    sample_summary <- sample_summary |>
        dplyr::left_join(ambiguity_summary, by="Init_Sample_ID") |>
        dplyr::mutate(Multiple_Valid_Solutions = case_when(
            !is.na(All_Valid_Subject_IDs) > 0 ~ "multiple_valid_subjects",
            !is.na(All_Valid_Sample_IDs) > 0 ~ "one_valid_subject_multiple_valid_samples"))

    sample_summary$Sample_Contamination_Metric_Numerator <- NA_integer_
    sample_summary$Sample_Contamination_Metric_Denominator <- NA_integer_
    sample_summary$Sample_Contamination_Metric <- NA
    if (!is.null(object@genotype_matrix)) {
        sample_summary$Sample_Contamination_Metric_Numerator <- NA_integer_
        sample_summary$Sample_Contamination_Metric_Denominator <- NA_integer_
        sample_summary$Sample_Contamination_Metric <- NA
        genotyped_sample_ids <- sample_summary |> dplyr::filter(!Ghost) |> dplyr::pull(Init_Sample_ID)
        for (sample_id in genotyped_sample_ids) {
            neighbor_samples <- names(which(object@genotype_matrix[sample_id, ] == 1))
            neighbor_matrix <- object@genotype_matrix[neighbor_samples, neighbor_samples]
            existing_edges <- sum(neighbor_matrix[upper.tri(neighbor_matrix)])
            total_edges <- length(neighbor_matrix[upper.tri(neighbor_matrix)])
            missing_edges <- total_edges - existing_edges
            sample_summary[sample_summary$Init_Sample_ID == sample_id, "Sample_Contamination_Metric_Numerator"] <- missing_edges
            sample_summary[sample_summary$Init_Sample_ID == sample_id, "Sample_Contamination_Metric_Denominator"] <- total_edges
            if (total_edges > 0) {
                sample_summary[sample_summary$Init_Sample_ID == sample_id, "Sample_Contamination_Metric"] <- missing_edges / total_edges
            }
        }
    }

    corrections_graph <- .generate_corrections_graph(object@.solve_state$relabel_data)
    corrections_components <- components(corrections_graph)
    Mislabeling_Event_ID_df <- data.frame
    Mislabeling_Event_ID_df <- data.frame(Init_Sample_ID = names(corrections_components$membership),
                                    Mislabeling_Event_ID = corrections_components$membership)
    sample_summary <- sample_summary |>
        dplyr::left_join(Mislabeling_Event_ID_df, by="Init_Sample_ID")
    Mislabeling_Event_ID_renamer_df <- sample_summary |>
        dplyr::filter(!is.na(Mislabeling_Event_ID)) |>
        dplyr::select(Component_ID, Mislabeling_Event_ID) |>
        dplyr::distinct() |>
        dplyr::group_by(Component_ID) |>
        dplyr::mutate(Event_Number = row_number()) |>
        dplyr::ungroup() |>
        dplyr::mutate(Mislabeling_Event_ID_New = paste0(Component_ID, "_Mislabeling_Event_", Event_Number)) |>
        dplyr::select(Mislabeling_Event_ID, Mislabeling_Event_ID_New)
    Mislabeling_Event_ID_renamer_map <- setNames(Mislabeling_Event_ID_renamer_df$Mislabeling_Event_ID_New, Mislabeling_Event_ID_renamer_df$Mislabeling_Event_ID)
    sample_summary$Mislabeling_Event_ID <- sapply(sample_summary$Mislabeling_Event_ID, \(x) ifelse(is.na(x), NA, Mislabeling_Event_ID_renamer_map[x]))

    sample_summary <- sample_summary |>
        dplyr::select(Connected_Component_ID = Component_ID, Genotype_Group_ID, SwapCat_ID, Is_Ghost = Ghost,
                      Initial_Subject_ID = Init_Subject_ID, Initial_Sample_ID = Init_Sample_ID,
                      Selected_For_Review, Mislabeled, Mislabeling_Event_ID, Multiple_Valid_Solutions, All_Valid_Subject_IDs, All_Valid_Sample_IDs,
                      Inferred_Subject_ID, Proposed_Final_Subject_ID, Proposed_Final_Sample_ID,
                      Sample_Count_In_Genotype_Group, Sample_Count_In_Genotype_Group_with_Same_Initial_Subject_Label,
                      Sample_Contamination_Metric, Sample_Contamination_Metric_Denominator, Sample_Contamination_Metric_Numerator)

    genotype_group_summary <- sample_summary |>
        dplyr::group_by(Genotype_Group_ID, Inferred_Subject_ID) |>
        dplyr::filter(!is.na(Genotype_Group_ID)) |>
        dplyr::summarize(
            # Majority_Subject_ID = names(sort(table(Proposed_Final_Subject_ID), decreasing = TRUE)[1]),
            n_Samples_no_review_needed = sum(Selected_For_Review == "no_review_needed"),
            n_Samples_inconsistent_genotype = sum(Selected_For_Review == "inconsistent_genotype"),
            n_Samples_deletion_or_duplication = sum(Selected_For_Review == "deletion_or_duplication"),
            n_Samples_relabel_low_confidence = sum(Selected_For_Review == "relabel_low_confidence"),
            n_Samples_relabel_high_confidence = sum(Selected_For_Review == "relabel_high_confidence"),
            n_Samples_singleton_no_inference = sum(Selected_For_Review == "singleton_no_inference"),
            n_Samples_not_relabeled_low_confidence = sum(Selected_For_Review == "not_relabeled_low_confidence"),
            n_Samples_total = n(),
            n_Samples_Initially_Matching_Inferred_Subject = sum(Initial_Subject_ID == Inferred_Subject_ID),
            Selected_For_Review = case_when(
                n_Samples_total == n_Samples_no_review_needed ~ "no_review_needed",
                n_Samples_total == n_Samples_singleton_no_inference ~ "singleton_no_inference",
                TRUE ~ "check_sample_table"
            )
        ) |>
        dplyr::ungroup() |>
        dplyr::select(Genotype_Group_ID, n_Samples_total, Inferred_Subject_ID, Selected_For_Review, n_Samples_Initially_Matching_Inferred_Subject, everything())

    genotype_group_summary$Genotype_Contamination_Metric <- NA
    genotype_group_summary$Genotype_Contamination_Metric_Denominator <- NA_integer_
    genotype_group_summary$Genotype_Contamination_Metric_Numerator <- NA_integer_
    if (!is.null(object@genotype_matrix)) {
        genotype_group_summary$Genotype_Contamination_Metric <- NA
        genotype_group_summary$Genotype_Contamination_Metric_Denominator <- NA_integer_
        genotype_group_summary$Genotype_Contamination_Metric_Numerator <- NA_integer_
        for (genotype_group_id in genotype_group_summary$Genotype_Group_ID) {
            genotype_group_init_samples <- sample_summary |>
                filter(Genotype_Group_ID == genotype_group_id) |>
                pull(Initial_Sample_ID)
            genotype_group_matrix <- object@genotype_matrix[genotype_group_init_samples, genotype_group_init_samples]
            existing_edges <- sum(genotype_group_matrix[upper.tri(genotype_group_matrix)])
            total_edges <- length(genotype_group_matrix[upper.tri(genotype_group_matrix)])
            missing_edges <- total_edges - existing_edges
            genotype_group_summary[genotype_group_summary$Genotype_Group_ID == genotype_group_id, "Genotype_Contamination_Metric_Numerator"] <- missing_edges
            genotype_group_summary[genotype_group_summary$Genotype_Group_ID == genotype_group_id, "Genotype_Contamination_Metric_Denominator"] <- total_edges
            if (total_edges != 0) {
                genotype_group_summary[genotype_group_summary$Genotype_Group_ID == genotype_group_id, "Genotype_Contamination_Metric"] <- missing_edges / total_edges
            }
        }
    }

    component_summary <- sample_summary |>
        dplyr::group_by(Connected_Component_ID) |>
        dplyr::summarize(
            n_Genotype_Groups = length(unique(Genotype_Group_ID)),
            n_Subjects = length(unique(Initial_Subject_ID)),
            n_Samples_total = n(),
            Same_Number_of_Genotypes_And_Subjects = n_Genotype_Groups == n_Subjects,
            n_Samples_no_review_needed = sum(Selected_For_Review == "no_review_needed"),
            n_Samples_ghost_relabeled = sum(Selected_For_Review == "ghost_relabeled"),
            n_Samples_ghost = sum(Selected_For_Review == "ghost"),
            n_Samples_inconsistent_genotype = sum(Selected_For_Review == "inconsistent_genotype"),
            n_Samples_deletion_or_duplication = sum(Selected_For_Review == "deletion_or_duplication"),
            n_Samples_relabel_low_confidence = sum(Selected_For_Review == "relabel_low_confidence"),
            n_Samples_relabel_high_confidence = sum(Selected_For_Review == "relabel_high_confidence"),
            n_Samples_singleton_no_inference = sum(Selected_For_Review == "singleton_no_inference"),
            n_Samples_not_relabeled_low_confidence = sum(Selected_For_Review == "not_relabeled_low_confidence"),
            Selected_For_Review = case_when(
                n_Samples_total == n_Samples_no_review_needed ~ "no_review_needed",
                n_Samples_total == n_Samples_singleton_no_inference ~ "singleton_no_inference",
                TRUE ~ "check_sample_table"
            )
        )

    dataset_summary <- sample_summary |>
        dplyr::summarize(
            n_Components = length(unique(Connected_Component_ID)),
            n_Genotype_Groups = length(unique(Genotype_Group_ID)),
            n_Subjects = length(unique(Initial_Subject_ID)),
            n_Samples_total = n(),
            n_Samples_no_review_needed = sum(Selected_For_Review == "no_review_needed"),
            n_Samples_ghost_relabeled = sum(Selected_For_Review == "ghost_relabeled"),
            n_Samples_ghost = sum(Selected_For_Review == "ghost"),
            n_Samples_inconsistent_genotype = sum(Selected_For_Review == "inconsistent_genotype"),
            n_Samples_deletion_or_duplication = sum(Selected_For_Review == "deletion_or_duplication"),
            n_Samples_relabel_low_confidence = sum(Selected_For_Review == "relabel_low_confidence"),
            n_Samples_relabel_high_confidence = sum(Selected_For_Review == "relabel_high_confidence"),
            n_Samples_singleton_no_inference = sum(Selected_For_Review == "singleton_no_inference"),
            n_Samples_not_relabeled_low_confidence = sum(Selected_For_Review == "not_relabeled_low_confidence")
        )

    summary_list <- list(
        "Sample" = sample_summary,
        "Genotype_Group" = genotype_group_summary,
        "Component" = component_summary,
        "Dataset" = dataset_summary)
    return(summary_list)

    excel_filename <- file.path(dir_path, "corrections_summary.xlsx")
    openxlsx::write.xlsx(summary_list, file=excel_filename)
    message(paste0("Output successfully written to ", excel_filename))
}

