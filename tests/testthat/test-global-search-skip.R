test_that("global search is skipped once nothing new is available to it, and solving still fully converges", {
    scenario <- big_tie_cycle_scenario(10)
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)

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
    expect_true(any(grepl("Starting global search", msgs)))
    ## Despite all the skipping, everything should still end up solved.
    expect_equal(sum(!solved@.solve_state$relabel_data$Solved), 0)
})

test_that("global search is not skipped on an attempt where something is genuinely available", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)

    ## toy_swap_scenario() is immediately solvable by global search alone,
    ## so its very first attempt (nothing has run yet, so anything
    ## available at all is new) must be a real call, never a skip.
    msgs <- capture_messages(solved <- solveEnsemble(x, use_solvers = "global"))
    first_global_msg <- msgs[grepl("global search", msgs, ignore.case = TRUE)][
        1
    ]
    expect_match(first_global_msg, "Starting global search")
})
