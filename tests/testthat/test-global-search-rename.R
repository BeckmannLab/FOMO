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

test_that("use_solvers rejects mixing 'global' and 'global_fast'", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )

    expect_error(
        solveEnsemble(x, use_solvers = c("majority", "global", "global_fast")),
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
