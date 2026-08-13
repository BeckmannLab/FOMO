test_that("'local_old' is an accepted use_solvers value, and rejects mixing with 'local'", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    expect_error(
        solveEnsemble(x, use_solvers = c("majority", "local", "local_old")),
        "cannot contain both"
    )
    expect_no_error(
        suppressMessages(solveEnsemble(
            x,
            use_solvers = c("majority", "global", "local_old")
        ))
    )
})

test_that("use_solvers = 'local_old' in solveEnsemble() reaches the same correct resolution as the default", {
    scenario <- toy_swap_scenario()
    x1 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
    x2 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    r1 <- suppressMessages(solveEnsemble(
        x1,
        use_solvers = c("majority", "global", "local_old")
    ))
    r2 <- suppressMessages(solveEnsemble(x2))

    expect_identical(
        r1@.solve_state$relabel_data[, c(
            "Init_Sample_ID",
            "Subject_ID",
            "Genotype_Group_ID"
        )],
        r2@.solve_state$relabel_data[, c(
            "Init_Sample_ID",
            "Subject_ID",
            "Genotype_Group_ID"
        )]
    )
})

test_that("default use_solvers runs without any deprecation warnings", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
    expect_no_warning(suppressMessages(solveEnsemble(x)))
})
