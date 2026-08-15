test_that("solveLocalSearch() has the expected min_genotypes default", {
    expect_equal(formals(solveLocalSearch)$min_genotypes, 1)
})

test_that("solveLocalSearch() resolves toy_swap_scenario() under the default min_genotypes", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)

    solved <- suppressMessages(solveLocalSearch(x))

    expect_equal(sum(!solved@.solve_state$relabel_data$Solved), 0)
})

test_that("solveLocalSearch()'s min_genotypes skips components with too few unsolved genotypes", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)
    n_unsolved_before <- nrow(x@.solve_state$unsolved_relabel_data)

    ## toy_swap_scenario()'s single component only has 2 unsolved genotype
    ## groups (G1, G2), below a min_genotypes = 3 threshold, so it should be
    ## skipped entirely: nothing gets relabeled.
    solved <- suppressMessages(solveLocalSearch(x, min_genotypes = 3))

    expect_equal(
        nrow(solved@.solve_state$unsolved_relabel_data),
        n_unsolved_before
    )
})
