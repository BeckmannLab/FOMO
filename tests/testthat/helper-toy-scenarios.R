## Shared toy-scenario constructors used across test files.

## G1 is torn between two real samples (reporting S1 and S2); G2 has one real
## sample reporting S1, linking G1 and G2 into a single component. The
## lowest-cost resolution swaps sample2 (G1, currently S2) with sample3 (G2,
## currently S1): sample2 -> S1, sample3 -> S2. This is also the scenario
## originally used to diagnose and fix the global-search deletion-penalty
## double-counting bug: the *other* permutation (G1->S2, G2->S1) double-counted
## the resulting orphaned sample before that fix.
toy_swap_scenario <- function() {
    sample_metadata <- data.frame(
        Sample_ID = c("sample1", "sample2", "sample3"),
        Subject_ID = c("S1", "S2", "S1"),
        Genotype_Group_ID = c("G1", "G1", "G2"),
        stringsAsFactors = FALSE
    )
    list(sample_metadata = sample_metadata)
}

## G1 has two real samples, both incorrectly reporting S2; G2 has one real
## sample correctly reporting S2 (so S2 "has a genotype" per G2). No real or
## ghost sample anywhere reports S1 (G1's true subject), so BOTH of G1's samples
## can only be resolved by inventing a placeholder ID for a duplicate that
## doesn't exist ("unknown"/LABELNOTFOUND labels).
toy_deficit_scenario <- function(shuffle = FALSE) {
    sample_metadata <- data.frame(
        Sample_ID = c("sample1", "sample2", "sample3"),
        Subject_ID = c("S2", "S2", "S2"),
        Genotype_Group_ID = c("G1", "G1", "G2"),
        stringsAsFactors = FALSE
    )
    if (shuffle) {
        # We don't actually shuffle randomly because that has a non-zero chance
        # to return the original order.
        ord <- rev(seq_len(nrow(sample_metadata)))
        sample_metadata <- sample_metadata[ord, ]
    }
    list(sample_metadata = sample_metadata)
}

## Two genotype groups (G1, G2) each with a 2-1 majority pointing at one
## subject: G1's 3 samples mostly report S1 (2 real + 1 mislabeled reporting
## S2), and G2's 3 samples mostly report S2 (2 real + 1 mislabeled reporting
## S1). Both directions of the majority vote agree (a majority of G1's samples
## report S1 *and* a majority of S1's samples report G1, and likewise for
## G2/S2), so solveMajoritySearch() alone can lock G1<->S1 and G2<->S2 without
## needing global or local search, then resolve the two mislabeled samples by
## swapping their Subject_ID labels with each other. All 6 samples share a
## single connected component with exactly 2 (as yet unlocked) genotype groups
## and 2 subjects, which doubles this scenario as a minimal case for exercising
## solveMajoritySearch()'s max_genotypes argument: max_genotypes = 1 must skip
## the whole component (2 > 1), while the default (100) must not.
toy_majority_scenario <- function() {
    sample_metadata <- data.frame(
        Sample_ID = c(
            "sample1",
            "sample2",
            "sample3",
            "sample4",
            "sample5",
            "sample6"
        ),
        Subject_ID = c("S1", "S1", "S2", "S2", "S2", "S1"),
        Genotype_Group_ID = c("G1", "G1", "G1", "G2", "G2", "G2"),
        stringsAsFactors = FALSE
    )
    list(sample_metadata = sample_metadata)
}

# Toy example with k subjects, 3 samples each, with a single cycle involving 1
# sample from each subject, such that every subject has a majority, leading to a
# single unambiguous best solution.
toy_cycle_scenario <- function(k = 10) {
    subj <- sprintf("S%02d", seq_len(k))
    geno <- sprintf("G%02d", seq_len(k))
    a_rows <- tibble(
        Sample_ID = sprintf("a%02d", seq_len(k)),
        Subject_ID = subj,
        Correct_Subject_ID = subj,
        Actually_Mislabeled = FALSE,
        Genotype_Group_ID = geno
    )
    b_rows <- tibble(
        Sample_ID = sprintf("b%02d", seq_len(k)),
        Subject_ID = subj[c(2:k, 1)],
        Correct_Subject_ID = subj,
        Actually_Mislabeled = TRUE,
        Genotype_Group_ID = geno
    )
    c_rows <- tibble(
        Sample_ID = sprintf("c%02d", seq_len(k)),
        Subject_ID = subj,
        Correct_Subject_ID = subj,
        Actually_Mislabeled = FALSE,
        Genotype_Group_ID = geno
    )
    sample_metadata <- rbind(a_rows, b_rows, c_rows)
    list(sample_metadata = sample_metadata)
}

## A "tie cycle": k genotype groups <-> k subjects, where every genotype group
## has an exact 1-1 tie between two candidate subjects (so majority search can
## never resolve any of it) and the whole cycle forms a single connected
## component. With k > global_max_genotypes (default 8), the entire component is
## permanently ineligible for global search, so only local search can make
## progress on it -- this is the situation where solveEnsemble() should keep
## skipping global search across outer-loop iterations until the component is
## fully resolved by local search alone.
big_tie_cycle_scenario <- function(k = 10) {
    subj <- sprintf("S%02d", seq_len(k))
    geno <- sprintf("G%02d", seq_len(k))
    a_rows <- data.frame(
        Sample_ID = sprintf("a%02d", seq_len(k)),
        Subject_ID = subj,
        Genotype_Group_ID = geno,
        stringsAsFactors = FALSE
    )
    b_rows <- data.frame(
        Sample_ID = sprintf("b%02d", seq_len(k)),
        Subject_ID = subj[c(2:k, 1)],
        Genotype_Group_ID = geno,
        stringsAsFactors = FALSE
    )
    sample_metadata <- rbind(a_rows, b_rows)
    label_domains <- data.frame(
        Sample_ID = sample_metadata$Sample_ID,
        Label_Domain = "omic1",
        stringsAsFactors = FALSE
    )
    list(sample_metadata = sample_metadata, label_domains = label_domains)
}
