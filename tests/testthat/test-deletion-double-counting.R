test_that("global search does not double-count a single orphaned sample's deletion penalty", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
    solved <- suppressMessages(solveGlobalSearch(x))
    rd <- solved@.solve_state$relabel_data[, c(
        "Init_Sample_ID",
        "Subject_ID",
        "Genotype_Group_ID"
    )]
    rd <- rd[order(rd$Init_Sample_ID), ]
    rownames(rd) <- NULL

    ## The correct (lowest-cost) resolution: sample2 and sample3 swap
    ## subjects. Before the double-counting fix, this scenario's *other*
    ## permutation (leaving sample1/sample3 tangled) was artificially
    ## preferred at deletion_penalty=2 because the resulting orphaned
    ## sample's deletion was counted twice.
    expected <- data.frame(
        Init_Sample_ID = c("sample1", "sample2", "sample3"),
        Subject_ID = c("S1", "S1", "S2"),
        Genotype_Group_ID = c("G1", "G1", "G2"),
        stringsAsFactors = FALSE
    )
    expect_equal(rd, expected, ignore_attr = TRUE)
})
