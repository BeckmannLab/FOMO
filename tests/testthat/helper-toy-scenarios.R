## Shared toy-scenario constructors used across test files.

## G1 is torn between two real samples (reporting S1 and S2); G2 has one
## real sample reporting S1, linking G1 and G2 into a single component. The
## lowest-cost resolution swaps sample2 (G1, currently S2) with sample3
## (G2, currently S1): sample2 -> S1, sample3 -> S2. This is also the
## scenario originally used to diagnose and fix the global-search
## deletion-penalty double-counting bug: the *other* permutation (G1->S2,
## G2->S1) double-counted the resulting orphaned sample before that fix.
toy_swap_scenario <- function(ids = c("sample1", "sample2", "sample3")) {
    sample_metadata <- data.frame(
        Sample_ID = ids,
        Subject_ID = c("S1", "S2", "S1"),
        Genotype_Group_ID = c("G1", "G1", "G2"),
        stringsAsFactors = FALSE
    )
    label_domains <- data.frame(
        Sample_ID = ids,
        Label_Domain = "omic1",
        stringsAsFactors = FALSE
    )
    list(sample_metadata = sample_metadata, label_domains = label_domains)
}

## G1 has two real samples, both incorrectly reporting S2; G2 has one real
## sample correctly reporting S2 (so S2 "has a genotype" per G2). No real or
## ghost sample anywhere reports S1 (G1's true subject), so BOTH of G1's
## samples can only be resolved by inventing a placeholder ID for a
## duplicate that doesn't exist ("unknown"/LABELNOTFOUND labels).
toy_deficit_scenario <- function(
    ids = c("sample1", "sample2", "sample3"),
    shuffle = FALSE
) {
    sample_metadata <- data.frame(
        Sample_ID = ids,
        Subject_ID = c("S2", "S2", "S2"),
        Genotype_Group_ID = c("G1", "G1", "G2"),
        stringsAsFactors = FALSE
    )
    label_domains <- data.frame(
        Sample_ID = ids,
        Label_Domain = "omic1",
        stringsAsFactors = FALSE
    )
    if (shuffle) {
        ord <- rev(seq_len(nrow(sample_metadata)))
        sample_metadata <- sample_metadata[ord, ]
        label_domains <- label_domains[ord, ]
    }
    list(sample_metadata = sample_metadata, label_domains = label_domains)
}

## Two genotype groups (G1, G2) each with a 2-1 majority pointing at one
## subject: G1's 3 samples mostly report S1 (2 real + 1 mislabeled
## reporting S2), and G2's 3 samples mostly report S2 (2 real + 1
## mislabeled reporting S1). Both directions of the majority vote agree (a
## majority of G1's samples report S1 *and* a majority of S1's samples
## report G1, and likewise for G2/S2), so solveMajoritySearch() alone can
## lock G1<->S1 and G2<->S2 without needing global or local search, then
## resolve the two mislabeled samples by swapping their Subject_ID labels
## with each other. All 6 samples share a single connected component with
## exactly 2 (as yet unlocked) genotype groups and 2 subjects, which
## doubles this scenario as a minimal case for exercising
## solveMajoritySearch()'s max_genotypes argument: max_genotypes = 1 must
## skip the whole component (2 > 1), while the default (100) must not.
toy_majority_scenario <- function(
    ids = c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6")
) {
    sample_metadata <- data.frame(
        Sample_ID = ids,
        Subject_ID = c("S1", "S1", "S2", "S2", "S2", "S1"),
        Genotype_Group_ID = c("G1", "G1", "G1", "G2", "G2", "G2"),
        stringsAsFactors = FALSE
    )
    label_domains <- data.frame(
        Sample_ID = ids,
        Label_Domain = "omic1",
        stringsAsFactors = FALSE
    )
    list(sample_metadata = sample_metadata, label_domains = label_domains)
}
