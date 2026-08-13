## A "tie cycle": k genotype groups <-> k subjects, where every genotype
## group has an exact 1-1 tie between two candidate subjects (so majority
## search can never resolve any of it) and the whole cycle forms a single
## connected component. With k > global_max_genotypes (default 8), the
## entire component is permanently ineligible for global search, so only
## local search can make progress on it -- this is the situation where
## solveEnsemble() should keep skipping global search across outer-loop
## iterations until the component is fully resolved by local search alone.
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
    swap_cats <- data.frame(
        Sample_ID = sample_metadata$Sample_ID,
        SwapCat_ID = "omic1",
        stringsAsFactors = FALSE
    )
    list(sample_metadata = sample_metadata, swap_cats = swap_cats)
}

test_that("global search is skipped once nothing new is available to it, and solving still fully converges", {
    scenario <- big_tie_cycle_scenario(10)
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    msgs <- capture_messages(
        solved <- solveEnsemble(
            x,
            use_solvers = c("majority", "global", "local"),
            time_limit = 20
        )
    )

    expect_true(any(grepl("Skipping global search", msgs)))
    ## Skipping is never the *only* thing that happens to global search --
    ## it should still run for real at least once (the component's initial,
    ## structural too-large-for-global check plus its final resolution).
    expect_true(any(grepl("^Starting global search", msgs)))
    ## Despite all the skipping, everything should still end up solved.
    expect_equal(sum(!solved@.solve_state$relabel_data$Solved), 0)
})

test_that("global search is not skipped on an attempt where something is genuinely available", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    ## toy_swap_scenario() is immediately solvable by global search alone,
    ## so its very first attempt (nothing has run yet, so anything
    ## available at all is new) must be a real call, never a skip.
    msgs <- capture_messages(solved <- solveEnsemble(x, use_solvers = "global"))
    first_global_msg <- msgs[grepl("global search", msgs, ignore.case = TRUE)][
        1
    ]
    expect_match(first_global_msg, "^Starting global search")
})
