# This is a workaround for the note "No visible binding for global variable"
# https://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
#
# The list below is generated, not hand-maintained: it's every name reported
# by running codetools::checkUsagePackage("fomo", all = TRUE) (after
# pkgload::load_all()) under either of the two NOTE categories
# globalVariables() actually suppresses -- "no visible binding for global
# variable" (bare column names used inside dplyr/NSE contexts) and "no
# visible global function definition for" (this second category is just
# ".from", igraph's edge-selector syntax used in .generate_graph() --
# confirmed it isn't a real function at all: `.from` doesn't exist anywhere
# in igraph's namespace, exported or internal, yet `E(g)[.from(1)]` still
# works, because `[.igraph.es` recognizes the literal, unevaluated call as
# special selector syntax rather than ever actually calling a function named
# .from(). That's a well-known false-positive pattern when statically
# analyzing igraph-based code, and globalVariables() is the standard,
# documented way to suppress it).
#
# If this list goes stale again, regenerate it the same way rather than
# hand-editing -- that's what happened to the previous version of this file:
# it had drifted enough that ~20 more bindings had accumulated in the
# package with no corresponding entries here, while globalVariables() itself
# had been commented out at some point and stopped running entirely, so
# neither the old list nor a fresh one was actually being applied.
#
# The more thorough fix -- eliminating the need for this file entirely by
# following the dplyr non-standard-evaluation vignette's .data$-prefix
# pattern throughout the package instead --
# https://cran.r-project.org/web/packages/dplyr/vignettes/in-packages.html
# -- remains a larger, separate undertaking (on the order of a hundred call
# sites across most of the package's R files) and hasn't been done.
global_vars <- c(".from", "All_Valid_Sample_IDs", "All_Valid_Subject_IDs", "Col",
                  "Component_ID", "Connected_Component_ID", "Curr_Subject_ID_Genotyped", "Deleted_relabel_from",
                  "delta", "Event_Number", "Genotype_Group_A", "Genotype_Group_B",
                  "Genotype_Group_ID", "Ghost", "Inferred_Correctly_Labeled", "Inferred_Subject_ID",
                  "Init_Component_ID", "Init_Sample_ID", "Init_Subject_ID", "Initial_Sample_ID",
                  "Initial_Subject_ID", "Invalid_Swap", "Is_Anchor", "Is_Ghost",
                  "Max_Genotype_Group_ID", "Max_Subject_ID", "Mislabeled", "Mislabeling_Event_ID",
                  "Mislabeling_Event_ID_New", "Multiple_Valid_Solutions", "n_agree", "n_genotype_deletions",
                  "n_Genotype_Group_ID", "n_Genotype_Groups", "n_ghost_labels", "n_in_genotype",
                  "n_label_deletions", "n_labels", "n_Max_Genotype_Group_ID", "n_Max_Subject_ID",
                  "n_Sample_ID", "n_samples_correct", "n_Samples_Initially_Matching_Inferred_Subject", "n_samples_to_relabel",
                  "n_samples_to_relabel_ghost", "n_Samples_total", "n_Subject_ID", "n_Subjects",
                  "new_Component_ID", "Permutation_ID", "Placeholder_ID", "Proposed_Final_Sample_ID",
                  "Proposed_Final_Subject_ID", "Putative_Subject_A", "Putative_Subject_B", "Putative_Subject_ID",
                  "relabel_from", "Row", "sample_a", "Sample_A",
                  "sample_b", "Sample_B", "Sample_Contamination_Metric", "Sample_Contamination_Metric_Denominator",
                  "Sample_Contamination_Metric_Numerator", "Sample_Count_In_Genotype_Group", "Sample_Count_In_Genotype_Group_with_Same_Initial_Subject_Label", "Sample_ID",
                  "Sample_ID.y", "sample1", "sample2", "Selected_For_Review",
                  "Solved", "Subject_A", "Subject_B", "Subject_ID",
                  "Subject_ID_putative", "Subject_ID.x", "Subject_ID.y", "SwapCat_A",
                  "SwapCat_B", "SwapCat_ID", "SwapCat_Shape", "vertex_size_scalar",
                  "X")
globalVariables(global_vars)
