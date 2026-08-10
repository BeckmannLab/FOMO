test_that("solvers emit status text via message(), not print()/cat()", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata, swap_cats = scenario$swap_cats)
    expect_message(solveGlobalSearch(x), "Starting global search")
})

test_that("solver status messages are properly suppressible and do not leak to stdout", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata, swap_cats = scenario$swap_cats)
    expect_no_message(suppressMessages(solveGlobalSearch(x)))

    ## Assign the result inside capture.output()'s expression so only actual
    ## cat()/print() calls made *during* solving would show up here -- not
    ## the auto-print of the returned S4 object itself (an assignment's
    ## value is invisible, unlike a bare function call's).
    x2 <- MislabelSolver(sample_metadata = scenario$sample_metadata, swap_cats = scenario$swap_cats)
    out <- capture.output(result <- suppressMessages(solveGlobalSearch(x2)))
    expect_length(out, 0)
})
