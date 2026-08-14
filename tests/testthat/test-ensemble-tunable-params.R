test_that("global_ghost_penalty/global_deletion_penalty are passed through to solveGlobalSearch()'s validation", {
    scenario <- toy_swap_scenario()
    x1 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    x2 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )

    suppressWarnings({
        expect_warning(
            suppressMessages(solveEnsemble(
                x1,
                use_solvers = "global",
                global_ghost_penalty = 0.5
            )),
            "ghost_penalty"
        )
        expect_warning(
            suppressMessages(solveEnsemble(
                x2,
                use_solvers = "global",
                global_deletion_penalty = 1
            )),
            "deletion_penalty"
        )
    })
})

test_that("global_max_genotypes is passed through to solveGlobalSearch(), suppressing progress when too small", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    n_unsolved_before <- nrow(x@.solve_state$unsolved_relabel_data)

    solved <- suppressMessages(solveEnsemble(
        x,
        use_solvers = "global",
        global_max_genotypes = 0
    ))

    ## With max_genotypes = 0, every component is too large for global
    ## search to ever process (it requires free_genotypes/free_subjects <=
    ## max_genotypes), so nothing gets resolved and the loop stops after
    ## one no-progress iteration.
    expect_equal(
        nrow(solved@.solve_state$unsolved_relabel_data),
        n_unsolved_before
    )
})

test_that("local_iter_per_cycle is passed through to the local solver as n_iter", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    ## toy_swap_scenario() is fully resolved by global search alone, so by
    ## the time local search runs there is nothing left for it to do -- but
    ## it still announces its first requested iteration (out of
    ## local_iter_per_cycle) before noticing that and returning early, which
    ## is enough to confirm the argument was passed through correctly.
    expect_message(
        solveEnsemble(
            x,
            use_solvers = c("global", "local"),
            local_iter_per_cycle = 3
        ),
        "iteration \\(1 of 3\\)"
    )
})
