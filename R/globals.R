# This is a workaround for the note "No visible binding for global variable"
# https://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
#' @import utils
#'
global_vars <- c("Sample_ID", "Subject_ID", "Genotype_Group_ID", "Component_ID",
                 "SwapCat_ID", "Permutation_ID", "n_genotype_deletions", "n_ghost_labels",
                 "n_in_genotype", "n_label_deletions", "n_labels", "n_samples_correct",
                 "n_samples_to_relabel", "n_samples_to_relabel_ghost", "perm_score", ".from",
                 "Init_Component_ID", "Is_Ghost", "Max_Genotype_Group_ID", "Max_Subject_ID",
                 "n_Max_Genotype_Group_ID", "n_Max_Subject_ID", ".", "Col",
                 "Curr_Subject_ID_Genotyped", "Genotype_Group_A", "Genotype_Group_B",
                 "Inferred_Correctly_Labeled", "Invalid_Swap", "Putative_Subject_A",
                 "Putative_Subject_B", "Putative_Subject_ID", "Row", "Sample_A", "Sample_B",
                 "Subject_A", "Subject_B", "SwapCat_A", "SwapCat_B", "X", "delta",
                 "Is_Anchor", "SwapCat_Shape", "sample1", "sample2", "sample_a", "sample_b",
                 "vertex_size_scalar", "Deleted_relabel_from", "relabel_from", "Sample_ID.y",
                 "Subject_ID.x", "Subject_ID.y", "Init_Sample_ID", "relabel_from",
                 "Final_Sample_ID", "Final_Subject_ID", "Ghost", "Inferred_Subject_ID",
                 "Init_Fraction_Match", "Init_Subject_ID", "Solved", "Status", "Subject_ID_putative",
                 "n_Genotype_Group", "n_Genotype_Group_ID", "n_Genotype_Groups",
                 "n_Sample_ID", "n_Samples_ignored", "n_Samples_total", "n_Samples_validated",
                 "n_Subject_ID", "n_Subjects", "n_agree", "new_Component_ID")
utils::globalVariables(global_vars)
