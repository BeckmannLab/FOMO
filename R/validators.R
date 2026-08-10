.validate_sample_metadata <- function(sample_metadata, has_genotype_matrix=FALSE) {
    required_columns <- c("Sample_ID", "Subject_ID")
    if (!has_genotype_matrix) {
        required_columns <- c(required_columns, "Genotype_Group_ID")
    }
    missing_columns <- setdiff(required_columns, colnames(sample_metadata))
    assert_that(
        length(missing_columns) == 0,
        msg = paste0("'sample_metadata' is missing required column(s) ", paste(missing_columns, collapse=", "))
    )

    if (has_genotype_matrix) {
        assert_that(
            !("Genotype_Group_ID" %in% colnames(sample_metadata)),
            msg = paste0("'sample_metadata' should not contain a Genotype_Group_ID column if 'genotype_matrix' is provided")
        )
    }

    duplicated_samples <- sample_metadata$Sample_ID[duplicated(sample_metadata$Sample_ID)]
    assert_that(
        length(duplicated_samples) == 0,
        msg = paste0("'sample_metadata' has non-unique 'Sample_ID'(s) ", paste(duplicated_samples, collapse=", "))
    )

    assert_that(
        all(!is.na(sample_metadata$Sample_ID)),
        msg = "'sample_metadata' missing value(s) in 'Sample_ID' column"
    )

    assert_that(
        all(!is.na(sample_metadata$Subject_ID)),
        msg = "'sample_metadata' missing value(s) in 'Subject_ID' column"
    )
}

.validate_genotype_matrix <- function(genotype_matrix, sample_metadata) {
    rows_missing_in_cols <- setdiff(rownames(genotype_matrix), colnames(genotype_matrix))
    assert_that(
        length(rows_missing_in_cols) == 0,
        msg = paste0("'genotype_matrix' contains rowname(s) missing in colname(s) ", paste(rows_missing_in_cols, collapse=", "))
    )

    cols_missing_in_rows <- setdiff(colnames(genotype_matrix), rownames(genotype_matrix))
    assert_that(
        length(cols_missing_in_rows) == 0,
        msg = paste0("'genotype_matrix' contains colname(s) missing in rowname(s) ", paste(cols_missing_in_rows, collapse=", "))
    )

    assert_that(
        nrow(genotype_matrix) == ncol(genotype_matrix),
        msg = paste0("'genotype_matrix' is not square, nrow=", nrow(genotype_matrix), " and ncol=", ncol(genotype_matrix))
    )

    assert_that(
        isSymmetric(genotype_matrix),
        msg = paste0("'genotype_matrix' must be symmetric")
    )

    assert_that(
        sum(is.na(genotype_matrix)) == 0,
        msg = paste0("'genotype_matrix' contains NA value(s). If they represent samples with different genotypes, fill explicitly with zero or FALSE.")
    )

    assert_that(
        all(is.numeric(genotype_matrix)) | all(is.logical(genotype_matrix)),
        msg = paste0("Values in 'genotype_matrix' should all be logical or numeric")
    )

    assert_that(
        sum(abs(diag(genotype_matrix))) == 0,
        msg = paste0("'genotype_matrix' must have all zero or FALSE in the diagonal")
    )

    missing_samples <- setdiff(rownames(genotype_matrix), sample_metadata$Sample_ID)
    assert_that(
        length(missing_samples) == 0,
        msg = paste0("Row and column names of 'genotype_matrix' must come from Sample_ID(s) found in 'sample_metadata'. Check row/colname(s) ", paste(missing_samples, collapse=", "))
    )
}

.validate_swap_cats <- function(sample_metadata, swap_cats) {
    required_columns <- c("Sample_ID", "SwapCat_ID")
    missing_columns <- setdiff(required_columns, colnames(swap_cats))
    assert_that(
        length(missing_columns) == 0,
        msg = paste0("'swap_cats' is missing required column(s) ", paste(missing_columns, collapse=", "))
    )

    missing_values <- swap_cats[is.na(swap_cats$SwapCat_ID), "Sample_ID"]
    assert_that(
        length(missing_values) == 0,
        msg = paste0("'swap_cats' is missing SwapCat_ID for Sample_ID(s) ", paste(missing_values, collapse=", "))
    )

    duplicated_samples <- swap_cats$Sample_ID[duplicated(swap_cats$Sample_ID)]
    assert_that(
        length(duplicated_samples) == 0,
        msg = paste0("'swap_cats' has non-unique Sample_ID(s) ", paste(duplicated_samples, collapse=", "))
    )

    missing_samples <- setdiff(sample_metadata$Sample_ID, swap_cats$Sample_ID)
    assert_that(
        length(missing_samples) == 0,
        msg = paste0("'swap_cats' is missing Sample_ID(s) that exist in 'sample_metadata', check Sample_ID(s) ", paste(missing_samples, collapse=", "))
    )
}

.validate_anchor_samples <- function(sample_metadata, anchor_samples) {
    extra_samples <- setdiff(anchor_samples, sample_metadata$Sample_ID)
    assert_that(
        length(extra_samples) == 0,
        msg = paste0("'anchor_samples' contains Sample_ID(s) not in 'sample_metadata', check ", paste(extra_samples, collapse=", "))
    )

    ## Check that there are no cases where either
    ## 1. Two or more anchor samples with the different Subject_ID(s) have same Genotype_Group_ID
    ## 2. Two or more anchor samples with the same Subject_ID have different Genotype_Group_ID(s)
    anchor_samples_consistency <- data.frame(Sample_ID = anchor_samples) |>
        left_join(
            sample_metadata[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID")],
            by="Sample_ID"
        ) |>
        filter(!is.na(Genotype_Group_ID)) |>
        group_by(Genotype_Group_ID) |>
        mutate(
            n_Subject_ID = length(unique(Subject_ID))
        ) |>
        ungroup() |>
        group_by(Subject_ID) |>
        mutate(
            n_Genotype_Group_ID = length(unique(Genotype_Group_ID))
        ) |>
        ungroup()
    subject_inconsistent_samples <- anchor_samples_consistency |>
        filter(n_Subject_ID != 1) |>
        arrange(Sample_ID) |>
        pull(Sample_ID)
    genotype_inconsistent_samples <- anchor_samples_consistency |>
        filter(n_Genotype_Group_ID != 1) |>
        arrange(Sample_ID) |>
        pull(Sample_ID)
    assert_that(
        length(subject_inconsistent_samples) == 0,
        msg = paste0("'anchor_samples' contains Sample_ID(s) that have the same Genotype_Group_ID
                   but different Subject_ID(s), check ", paste(subject_inconsistent_samples, collapse=", "))
    )
    assert_that(
        length(genotype_inconsistent_samples) == 0,
        msg = paste0("'anchor_samples' contains Sample_ID(s) that have the same Subject_ID
                   but different Genotype_Group_ID(s), check ", paste(genotype_inconsistent_samples, collapse=", "))
    )
}

# --------------------------------------------------------------------------
# Don't need validators for internally created variables 'putative subjects' and 'relabels'


# .validate_putative_subjects <- function(sample_metadata, putative_subjects) {
#     required_columns <- c("Genotype_Group_ID", "Subject_ID")
#     missing_columns <- setdiff(required_columns, colnames(putative_subjects))
#     assert_that(
#         length(missing_columns) == 0,
#         msg = glue("'putative_subjects' is missing required column(s) {paste(missing_columns, collapse=\", \")}")
#     )
#
#     extra_genotype_groups <- setdiff(na.omit(putative_subjects$Genotype_Group_ID), sample_metadata$Genotype_Group_ID)
#     assert_that(
#         length(extra_genotype_groups) == 0,
#         msg = glue("'putative_subjects' has 'Genotype_Group_ID'(s) not found in 'sample_metadata', check {paste(extra_genotype_groups, collapse=\", \")}")
#     )
#     extra_subjects <- setdiff(na.omit(putative_subjects$Subject_ID), sample_metadata$Subject_ID)
#     assert_that(
#         length(extra_subjects) == 0,
#         msg = glue("'putative_subjects' has 'Subject_ID'(s) not found in 'sample_metadata', check {paste(extra_subjects, collapse=\", \")}")
#     )
#
#     duplicated_genotype_groups <- putative_subjects$Genotype_Group_ID[duplicated(na.omit(putative_subjects$Genotype_Group_ID))]
#     assert_that(
#         length(duplicated_genotype_groups) == 0,
#         msg = glue("'putative_subjects' does not map 'Subject_ID' to 'Genotype_Group_ID' one-to-one, check 'Genotype_Group'(s) {paste(duplicated_genotype_groups, collapse=\", \")}")
#     )
#     duplicated_subjects <- putative_subjects$Subject_ID[duplicated(na.omit(putative_subjects$Subject_ID))]
#     assert_that(
#         length(duplicated_subjects) == 0,
#         msg = glue("'putative_subjects' does not map 'Subject_ID' to 'Genotype_Group_ID' one-to-one, 'check Subject_ID'(s) {paste(duplicated_samples, collapse=\", \")}")
#     )
# }
#
# .validate_relabels <- function(sample_metadata, relabels) {
#     required_columns <- c("relabel_from", "relabel_to")
#     missing_columns <- setdiff(required_columns, names(relabels))
#     assert_that(
#         length(missing_columns) == 0,
#         msg = glue("'relabels' is missing required column(s) {paste(missing_columns, collapse=\", \")}")
#     )
#
#     missing_relabel_from <- setdiff(relabels$relabel_from, sample_metadata$Sample_ID)
#     assert_that(
#         length(missing_relabel_from) == 0,
#         msg = glue("'relabels$relabel_from' has Sample_ID(s) not in 'sample_metadata', check {paste(missing_relabel_from, collapse=\", \")}")
#     )
#     missing_relabel_to <- setdiff(relabels$relabel_to, sample_metadata$Sample_ID)
#     assert_that(
#         length(missing_relabel_to) == 0,
#         msg = glue("'relabels$relabel_to' has Sample_ID(s) not in 'sample_metadata', check {paste(missing_relabel_to, collapse=\", \")}")
#     )
#
#     duplicated_relabel_from <- relabels$relabel_from[duplicated(relabels$relabel_from)]
#     assert_that(
#         length(duplicated_relabel_from) == 0,
#         msg = glue("'relabels$relabel_from' has non-unique Sample_ID(s), check Sample_ID(s) {paste(duplicated_relabel_from, collapse=\", \")}")
#     )
#     duplicated_relabel_to <- relabels$relabel_to[duplicated(relabels$relabel_to)]
#     assert_that(
#         length(duplicated_relabel_to) == 0,
#         msg = glue("'relabels$relabel_to' has non-unique Sample_ID(s), check Sample_ID(s) {paste(duplicated_genotype_groups, collapse=\", \")}")
#     )
#
#     creating_duplicates <- setdiff(relabels$relabel_to, relabels$relabel_from)
#     assert_that(
#         length(creating_duplicates) == 0,
#         msg = glue("'relabels$relabel_to' contains Sample_ID(s) that don't exist in 'relabels$relabel_from', check Sample_ID(s) {paste(creating_duplicates, collapse=\", \")}")
#     )
# }
