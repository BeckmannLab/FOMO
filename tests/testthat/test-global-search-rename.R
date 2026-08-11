test_that("solveComprehensiveSearch() is deprecated but matches solveGlobalSearch()", {
    scenario <- toy_swap_scenario()
    x1 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
    x2 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    expect_warning(
        r1 <- solveComprehensiveSearch(x1),
        "renamed to 'solveGlobalSearch"
    )
    r2 <- suppressMessages(solveGlobalSearch(x2))

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

test_that("solveComprehensiveSearchFast() is deprecated but matches solveGlobalSearchFast()", {
    scenario <- toy_swap_scenario()
    x1 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
    x2 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    expect_warning(
        r1 <- solveComprehensiveSearchFast(x1),
        "renamed to 'solveGlobalSearchFast"
    )
    r2 <- suppressMessages(solveGlobalSearchFast(x2))

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

test_that("solveEnsemble(use_solvers = 'comprehensive') is a deprecated alias for 'global'", {
    scenario <- toy_swap_scenario()
    x1 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
    x2 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    expect_warning(
        r1 <- solveEnsemble(
            x1,
            use_solvers = c("majority", "comprehensive", "local")
        ),
        "deprecated"
    )
    r2 <- suppressMessages(solveEnsemble(
        x2,
        use_solvers = c("majority", "global", "local")
    ))

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

test_that("use_solvers rejects mixing 'global' and 'global_fast', including via the deprecated alias", {
    scenario <- toy_swap_scenario()
    x1 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
    x2 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    expect_error(
        solveEnsemble(x1, use_solvers = c("majority", "global", "global_fast")),
        "cannot contain both"
    )
    expect_error(
        suppressWarnings(solveEnsemble(
            x2,
            use_solvers = c("majority", "comprehensive", "global_fast")
        )),
        "cannot contain both"
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
