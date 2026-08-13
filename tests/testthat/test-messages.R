test_that("solvers emit status text via message(), not print()/cat()", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    expect_message(solveGlobalSearch(x), "Starting global search")
})

test_that("solver status messages are properly suppressible and do not leak to stdout", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    expect_no_message(suppressMessages(solveGlobalSearch(x)))

    ## Assign the result inside capture.output()'s expression so only actual
    ## cat()/print() calls made *during* solving would show up here -- not
    ## the auto-print of the returned S4 object itself (an assignment's
    ## value is invisible, unlike a bare function call's).
    x2 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    out <- capture.output(result <- suppressMessages(solveGlobalSearch(x2)))
    expect_length(out, 0)
})

test_that("solveLocalSearch() and solveLocalSearchOld() emit distinct status text", {
    scenario <- toy_swap_scenario()
    x1 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    x2 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    expect_message(solveLocalSearch(x1), "Starting local search\n$")
    expect_message(
        solveLocalSearchOld(x2),
        "Starting local search \\(old\\)\n"
    )
})
