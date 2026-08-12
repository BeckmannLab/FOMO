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

test_that("solveEnsemble() always runs global search once more after the loop ends, when requested and not timed out", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    msgs <- capture_messages(solved <- solveEnsemble(x, use_solvers = "global"))
    ## The last "Starting ___ search" message logged should be a global
    ## search call (the unconditional post-loop safeguard), even though
    ## everything was already fully solved before it ran.
    starting_msgs <- msgs[grepl("^Starting (global|majority|local)", msgs)]
    expect_match(tail(starting_msgs, 1), "^Starting global search")
})

test_that("the post-loop safeguard still runs even when every in-loop attempt was skipped", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    ## global_max_genotypes = 0 makes every component permanently
    ## ineligible, so both in-loop call sites are skipped from the very
    ## first attempt onward -- but the safeguard call is unconditional and
    ## should still happen exactly once at the end.
    msgs <- capture_messages(
        solved <- solveEnsemble(
            x,
            use_solvers = "global",
            global_max_genotypes = 0
        )
    )
    expect_equal(sum(grepl("Skipping global search", msgs)), 2)
    expect_equal(sum(grepl("^Starting global search", msgs)), 1)
})

test_that("the post-loop safeguard does not run when global search was not requested", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    msgs <- capture_messages(
        solved <- solveEnsemble(x, use_solvers = c("majority", "local"))
    )
    expect_false(any(grepl("global search", msgs, ignore.case = TRUE)))
})

test_that("the post-loop safeguard does not run when the time limit was reached", {
    scenario <- big_tie_cycle_scenario(10)
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    ## time_limit = 0 forces the very first time-check (before any solver
    ## has run at all) to trip, so there should be no solver messages of
    ## any kind, and in particular no post-loop safeguard call.
    expect_warning(
        {
            msgs <- capture_messages(
                solved <- solveEnsemble(
                    x,
                    use_solvers = c("majority", "global", "local"),
                    time_limit = 0
                )
            )
        },
        "reached 'time_limit'"
    )
    expect_length(msgs, 0)
})
