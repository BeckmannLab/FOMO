test_that("writeOutput() actually succeeds and produces a non-empty .xlsx file", {
    ## Regression test: writeOutput() used to unconditionally call
    ## writexl::xl_sheet() to freeze each sheet's header row, but that
    ## function has never existed in any released version of writexl, so
    ## every call to writeOutput() against a properly installed copy of the
    ## package (as opposed to one loaded via pkgload::load_all(), which
    ## tolerates the resulting NAMESPACE/import inconsistency) failed with
    ## "object 'xl_sheet' not found".
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    solved <- suppressMessages(solveEnsemble(x))

    out_file <- withr::local_tempfile(fileext = ".xlsx")
    expect_no_error(writeOutput(solved, out_file))
    expect_true(file.exists(out_file))
    expect_gt(file.info(out_file)$size, 0)
})

test_that("writeOutput() also accepts the result of collateOutput() directly", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    solved <- suppressMessages(solveEnsemble(x))
    collated <- collateOutput(solved)

    out_file <- withr::local_tempfile(fileext = ".xlsx")
    expect_no_error(writeOutput(collated, out_file))
    expect_true(file.exists(out_file))
    expect_gt(file.info(out_file)$size, 0)
})
