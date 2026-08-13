test_that("use_solvers with only 'local' overrides local_iter_per_cycle with a message", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    expect_message(
        solveEnsemble(x, use_solvers = "local", local_iter_per_cycle = 2),
        "overriding"
    )
    expect_message(
        solveEnsemble(x, use_solvers = "local", local_iter_per_cycle = 2),
        "iteration \\(1 of 100\\)"
    )
})

test_that("use_solvers with only 'local_old' also overrides local_iter_per_cycle", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    expect_message(
        solveEnsemble(x, use_solvers = "local_old", local_iter_per_cycle = 2),
        "overriding"
    )
})

test_that("the local-only override is not applied when 'global' or 'majority' is also requested", {
    scenario <- toy_swap_scenario()
    x1 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    x2 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )

    msgs1 <- testthat::capture_messages(
        solveEnsemble(x1, use_solvers = c("global", "local"))
    )
    msgs2 <- testthat::capture_messages(
        solveEnsemble(x2, use_solvers = c("majority", "local"))
    )
    expect_false(any(grepl("overriding", msgs1)))
    expect_false(any(grepl("overriding", msgs2)))
})

test_that("use_solvers with only local search still reaches the correct resolution", {
    scenario <- toy_swap_scenario()
    x1 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )
    x2 <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        label_domains = scenario$label_domains
    )

    r1 <- suppressMessages(solveEnsemble(x1, use_solvers = "local"))
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
