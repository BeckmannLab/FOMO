test_that("solvers emit status text via message(), not print()/cat()", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)
    expect_message(solveGlobalSearch(x), "Starting global search")
})

test_that("solver status messages are properly suppressible and do not leak to stdout", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)
    out <- capture_output(suppressMessages(solveGlobalSearch(x)))
    expect_equal(out, "")
    expect_no_message(suppressMessages(solveGlobalSearch(x)))
})

test_that("solveLocalSearch() and solveLocalSearchOld() emit distinct status text", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)
    expect_message(solveLocalSearch(x), "Starting local search\n$")
    expect_message(
        solveLocalSearchOld(x),
        "Starting local search \\(old\\)\n"
    )
})
