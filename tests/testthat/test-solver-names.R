test_that("solveComprehensiveSearch() and solveComprehensiveSearchFast() no longer exist", {
    expect_false(exists("solveComprehensiveSearch", mode = "function"))
    expect_false(exists("solveComprehensiveSearchFast", mode = "function"))
    expect_false("solveComprehensiveSearch" %in% getNamespaceExports("fomo"))
    expect_false(
        "solveComprehensiveSearchFast" %in% getNamespaceExports("fomo")
    )
})

test_that("solveEnsemble(use_solvers = 'comprehensive') is no longer accepted", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    expect_error(
        solveEnsemble(x, use_solvers = c("majority", "comprehensive", "local")),
        "must only contain values from"
    )
    expect_error(
        solveEnsemble(
            x,
            use_solvers = c("majority", "comprehensive_fast", "local")
        ),
        "must only contain values from"
    )
})

test_that("solveGlobalSearchFast() and solveLocalSearchFast() no longer exist, having been renamed", {
    expect_false(exists("solveGlobalSearchFast", mode = "function"))
    expect_false(exists("solveLocalSearchFast", mode = "function"))
    expect_false("solveGlobalSearchFast" %in% getNamespaceExports("fomo"))
    expect_false("solveLocalSearchFast" %in% getNamespaceExports("fomo"))
})

test_that("solveLocalSearchOld() exists and is exported (the renamed original local search algorithm)", {
    expect_true(exists("solveLocalSearchOld", mode = "function"))
    expect_true("solveLocalSearchOld" %in% getNamespaceExports("fomo"))
})

test_that("'global_fast' and 'local_fast' are no longer accepted use_solvers values", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    expect_error(
        solveEnsemble(x, use_solvers = c("majority", "global_fast", "local")),
        "must only contain values from"
    )
    expect_error(
        solveEnsemble(x, use_solvers = c("majority", "global", "local_fast")),
        "must only contain values from"
    )
})

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
