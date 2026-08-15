test_that("solveMajoritySearch() and solveEnsemble() have the expected max_genotypes default", {
    expect_equal(formals(solveMajoritySearch)$max_genotypes, 100)
    expect_equal(formals(solveEnsemble)$majority_max_genotypes, 100)
})

test_that("solveMajoritySearch() fully resolves toy_majority_scenario() under the default max_genotypes", {
    scenario <- toy_majority_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)

    solved <- suppressMessages(solveMajoritySearch(x))

    expect_equal(sum(!solved@.solve_state$relabel_data$Solved), 0)
    expect_equal(nrow(solved@.solve_state$putative_subjects), 2)
})

test_that("solveMajoritySearch()'s max_genotypes skips components that are too large", {
    scenario <- toy_majority_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)
    n_unsolved_before <- nrow(x@.solve_state$unsolved_relabel_data)

    solved <- suppressMessages(solveMajoritySearch(x, max_genotypes = 1))

    ## toy_majority_scenario()'s single component has 2 (as yet unlocked)
    ## genotype groups and 2 subjects, both over the max_genotypes = 1
    ## limit, so the component should be skipped entirely: nothing gets
    ## relabeled and no putative subjects get locked in.
    expect_equal(
        nrow(solved@.solve_state$unsolved_relabel_data),
        n_unsolved_before
    )
    expect_equal(nrow(solved@.solve_state$putative_subjects), 0)
})

test_that("majority_max_genotypes is passed through to solveMajoritySearch() by solveEnsemble()", {
    scenario <- toy_majority_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)
    n_unsolved_before <- nrow(x@.solve_state$unsolved_relabel_data)

    solved <- suppressMessages(solveEnsemble(
        x,
        use_solvers = "majority",
        majority_max_genotypes = 1
    ))

    ## With majority_max_genotypes = 1, majority search can never make
    ## progress on this scenario's single 2-genotype/2-subject component,
    ## so with only "majority" in use_solvers the ensemble loop should stop
    ## after one no-progress iteration, leaving the scenario untouched.
    expect_equal(
        nrow(solved@.solve_state$unsolved_relabel_data),
        n_unsolved_before
    )
})
