## These functions previously had no test coverage at all (confirmed: no
## existing test file referenced plot() or plotCorrections()), which is how
## a package-wide dplyr .data$-prefix pass could touch them without any
## test catching a mistake. Basic smoke tests here just confirm each code
## path still runs without error/warning after that pass.

test_that("plot() runs without error for unsolved and solved views", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
    expect_no_error(plot(x, unsolved = TRUE))
    expect_no_error(plot(x, unsolved = FALSE))
})

test_that("plot() runs without error with collapse_samples = TRUE", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
    expect_no_error(plot(x, unsolved = FALSE, collapse_samples = TRUE))
})

test_that("plot() runs without error when queried by each supported field", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
    expect_no_error(plot(
        x,
        unsolved = FALSE,
        query_by = "Init_Component_ID",
        query_val = "Component_1"
    ))
    expect_no_error(plot(
        x,
        unsolved = FALSE,
        query_by = "Subject_ID",
        query_val = "S1"
    ))
    expect_no_error(plot(
        x,
        unsolved = FALSE,
        query_by = "Genotype_Group_ID",
        query_val = "G1"
    ))
    expect_no_error(plot(
        x,
        unsolved = FALSE,
        query_by = "Sample_ID",
        query_val = "sample1"
    ))
})

test_that("plotCorrections() runs without error, with and without a query", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
    solved <- suppressMessages(solveEnsemble(x))
    expect_no_error(plotCorrections(solved))
    expect_no_error(plotCorrections(
        solved,
        query_by = "Subject_ID",
        query_val = "S1"
    ))
})
