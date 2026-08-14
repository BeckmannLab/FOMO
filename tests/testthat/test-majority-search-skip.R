test_that("majority search is skipped once nothing new is available to it", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )

    ## With majority_max_genotypes = 1, toy_swap_scenario()'s single
    ## 2-genotype component is too large for majority search from the very
    ## first check, so its available-sample set is empty both before and
    ## after -- majority search should be skipped immediately, leaving local
    ## search (global is excluded here so it can't resolve this on its own
    ## first) to resolve the scenario instead.
    msgs <- capture_messages(
        solved <- solveEnsemble(
            x,
            use_solvers = c("majority", "local"),
            majority_max_genotypes = 1
        )
    )

    expect_true(any(grepl("Skipping majority search", msgs)))
    expect_equal(sum(!solved@.solve_state$relabel_data$Solved), 0)
})

test_that("majority search is not skipped on an attempt where something is genuinely available", {
    scenario <- toy_majority_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )

    ## toy_majority_scenario() is immediately solvable by majority search
    ## alone, so its very first attempt (nothing has run yet, so anything
    ## available at all is new) must be a real call, never a skip.
    msgs <- capture_messages(solveEnsemble(x, use_solvers = "majority"))
    first_majority_msg <- msgs[
        grepl("majority search", msgs, ignore.case = TRUE)
    ][1]
    expect_match(first_majority_msg, "Starting majority search")
})
