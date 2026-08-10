test_that("placeholder IDs for unresolvable samples are reproducible for the same input", {
    scenario <- toy_deficit_scenario()
    x1 <- MislabelSolver(sample_metadata = scenario$sample_metadata, swap_cats = scenario$swap_cats)
    x2 <- MislabelSolver(sample_metadata = scenario$sample_metadata, swap_cats = scenario$swap_cats)
    putative_subjects <- data.frame(Genotype_Group_ID = c("G1", "G2"), Subject_ID = c("S1", "S2"))

    r1 <- .find_relabel_cycles_from_putative_subjects(
        x1@.solve_state$unsolved_relabel_data, putative_subjects,
        x1@.solve_state$unsolved_ghost_data, allow_unknowns = TRUE
    )
    r2 <- .find_relabel_cycles_from_putative_subjects(
        x2@.solve_state$unsolved_relabel_data, putative_subjects,
        x2@.solve_state$unsolved_ghost_data, allow_unknowns = TRUE
    )

    expect_true(any(grepl("LABELNOTFOUND", r1$relabel_to)))
    expect_identical(sort(r1$relabel_to), sort(r2$relabel_to))
})

test_that("placeholder IDs do not depend on input row order", {
    ordered <- toy_deficit_scenario(shuffle = FALSE)
    shuffled <- toy_deficit_scenario(shuffle = TRUE)
    x1 <- MislabelSolver(sample_metadata = ordered$sample_metadata, swap_cats = ordered$swap_cats)
    x2 <- MislabelSolver(sample_metadata = shuffled$sample_metadata, swap_cats = shuffled$swap_cats)
    putative_subjects <- data.frame(Genotype_Group_ID = c("G1", "G2"), Subject_ID = c("S1", "S2"))

    r1 <- .find_relabel_cycles_from_putative_subjects(
        x1@.solve_state$unsolved_relabel_data, putative_subjects,
        x1@.solve_state$unsolved_ghost_data, allow_unknowns = TRUE
    )
    r2 <- .find_relabel_cycles_from_putative_subjects(
        x2@.solve_state$unsolved_relabel_data, putative_subjects,
        x2@.solve_state$unsolved_ghost_data, allow_unknowns = TRUE
    )
    expect_identical(sort(r1$relabel_to), sort(r2$relabel_to))

    ## The constructed object itself stores sample_metadata sorted by
    ## Sample_ID regardless of the row order it was given in.
    expect_identical(x2@sample_metadata$Sample_ID, sort(x2@sample_metadata$Sample_ID))
})

test_that("placeholder IDs differ for genuinely different input", {
    scenario1 <- toy_deficit_scenario()
    scenario2 <- toy_deficit_scenario(ids = c("sampleA", "sampleB", "sampleC"))
    x1 <- MislabelSolver(sample_metadata = scenario1$sample_metadata, swap_cats = scenario1$swap_cats)
    x2 <- MislabelSolver(sample_metadata = scenario2$sample_metadata, swap_cats = scenario2$swap_cats)
    putative_subjects <- data.frame(Genotype_Group_ID = c("G1", "G2"), Subject_ID = c("S1", "S2"))

    r1 <- .find_relabel_cycles_from_putative_subjects(
        x1@.solve_state$unsolved_relabel_data, putative_subjects,
        x1@.solve_state$unsolved_ghost_data, allow_unknowns = TRUE
    )
    r2 <- .find_relabel_cycles_from_putative_subjects(
        x2@.solve_state$unsolved_relabel_data, putative_subjects,
        x2@.solve_state$unsolved_ghost_data, allow_unknowns = TRUE
    )
    expect_false(identical(sort(r1$relabel_to), sort(r2$relabel_to)))
})

test_that("constructing a MislabelSolver does not disturb the caller's RNG stream", {
    scenario <- toy_swap_scenario()
    set.seed(999)
    expected_after <- runif(3)

    set.seed(999)
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata, swap_cats = scenario$swap_cats)
    actual_after <- runif(3)

    expect_identical(expected_after, actual_after)
})
