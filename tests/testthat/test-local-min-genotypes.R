test_that("solveLocalSearch() has the expected min_genotypes default", {
    expect_equal(formals(solveLocalSearch)$min_genotypes, 1)
})

test_that("solveLocalSearch() resolves toy_swap_scenario() under the default min_genotypes", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )

    solved <- suppressMessages(solveLocalSearch(x))

    expect_equal(sum(!solved@.solve_state$relabel_data$Solved), 0)
})

test_that("solveLocalSearch()'s min_genotypes skips components with too few unsolved genotypes", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
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

test_that("solveEnsemble() scales solveLocalSearch()'s min_genotypes to whichever of global/majority search are active", {
    scenario <- toy_swap_scenario()
    x1 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    x2 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    n_unsolved_before <- nrow(x1@.solve_state$unsolved_relabel_data)

    ## toy_swap_scenario()'s single component has an exact 1-1 tie (no
    ## majority possible anywhere), so majority search alone can never
    ## resolve it. With "global" excluded from use_solvers, local search is
    ## the only solver that actually could -- but the default
    ## majority_max_genotypes (100) sets its min_genotypes to 101, well
    ## above this component's 2 unsolved genotypes, so it gets skipped and
    ## the scenario is left unsolved.
    left_unsolved <- suppressMessages(solveEnsemble(
        x1,
        use_solvers = c("majority", "local")
    ))
    expect_equal(
        nrow(left_unsolved@.solve_state$unsolved_relabel_data),
        n_unsolved_before
    )

    ## Lowering majority_max_genotypes to 1 drops local search's
    ## min_genotypes to 2, which no longer exceeds this component's 2
    ## unsolved genotypes, so local search is now allowed to resolve it.
    resolved <- suppressMessages(solveEnsemble(
        x2,
        use_solvers = c("majority", "local"),
        majority_max_genotypes = 1
    ))
    expect_equal(sum(!resolved@.solve_state$relabel_data$Solved), 0)
})
