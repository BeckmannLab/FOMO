test_that("placeholder IDs for unresolvable samples are reproducible for the same input", {
    scenario <- toy_deficit_scenario()
    x1 <- MislabelSolver(sample_metadata = scenario$sample_metadata)
    x2 <- MislabelSolver(sample_metadata = scenario$sample_metadata)
    putative_subjects <- data.frame(
        Genotype_Group_ID = c("G1", "G2"),
        Subject_ID = c("S1", "S2")
    )

    r1 <- .find_relabel_cycles_from_putative_subjects(
        x1@.solve_state$unsolved_relabel_data,
        putative_subjects,
        x1@.solve_state$unsolved_ghost_data,
        allow_unknowns = TRUE
    )
    r2 <- .find_relabel_cycles_from_putative_subjects(
        x2@.solve_state$unsolved_relabel_data,
        putative_subjects,
        x2@.solve_state$unsolved_ghost_data,
        allow_unknowns = TRUE
    )

    expect_true(any(grepl("LABELNOTFOUND", r1$relabel_to)))
    expect_identical(sort(r1$relabel_to), sort(r2$relabel_to))
})

test_that("placeholder IDs do not depend on input row order", {
    ordered <- toy_deficit_scenario(shuffle = FALSE)
    shuffled <- toy_deficit_scenario(shuffle = TRUE)
    x1 <- MislabelSolver(sample_metadata = ordered$sample_metadata)
    x2 <- MislabelSolver(sample_metadata = shuffled$sample_metadata)
    putative_subjects <- data.frame(
        Genotype_Group_ID = c("G1", "G2"),
        Subject_ID = c("S1", "S2")
    )

    r1 <- .find_relabel_cycles_from_putative_subjects(
        x1@.solve_state$unsolved_relabel_data,
        putative_subjects,
        x1@.solve_state$unsolved_ghost_data,
        allow_unknowns = TRUE
    )
    r2 <- .find_relabel_cycles_from_putative_subjects(
        x2@.solve_state$unsolved_relabel_data,
        putative_subjects,
        x2@.solve_state$unsolved_ghost_data,
        allow_unknowns = TRUE
    )
    expect_identical(sort(r1$relabel_to), sort(r2$relabel_to))

    ## The constructed object itself stores sample_metadata sorted by
    ## Sample_ID regardless of the row order it was given in.
    expect_identical(
        x2@sample_metadata$Sample_ID,
        sort(x2@sample_metadata$Sample_ID)
    )
})

test_that("placeholder IDs differ for genuinely different input", {
    scenario1 <- toy_deficit_scenario()
    scenario2 <- scenario1
    scenario2$sample_metadata$Sample_ID <- c("sampleA", "sampleB", "sampleC")
    x1 <- MislabelSolver(sample_metadata = scenario1$sample_metadata)
    x2 <- MislabelSolver(sample_metadata = scenario2$sample_metadata)
    putative_subjects <- data.frame(
        Genotype_Group_ID = c("G1", "G2"),
        Subject_ID = c("S1", "S2")
    )

    r1 <- .find_relabel_cycles_from_putative_subjects(
        x1@.solve_state$unsolved_relabel_data,
        putative_subjects,
        x1@.solve_state$unsolved_ghost_data,
        allow_unknowns = TRUE
    )
    r2 <- .find_relabel_cycles_from_putative_subjects(
        x2@.solve_state$unsolved_relabel_data,
        putative_subjects,
        x2@.solve_state$unsolved_ghost_data,
        allow_unknowns = TRUE
    )
    expect_false(identical(sort(r1$relabel_to), sort(r2$relabel_to)))
})

test_that("constructing a MislabelSolver does not disturb the caller's RNG stream", {
    scenario <- toy_swap_scenario()
    set.seed(999)
    expected_after <- runif(3)

    set.seed(999)
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)
    actual_after <- runif(3)

    expect_identical(expected_after, actual_after)
})

test_that("an invented placeholder label always embeds its own row's Placeholder_ID", {
    ## G1 is genetically a single subject's four samples: g1a happens to
    ## still carry a correct-looking label (S1), while g1b/g1c/g1d each
    ## carry a *different* other subject's label (S2/S3/S4) -- i.e. three
    ## other subjects each appear to have a duplicate sample that is really
    ## G1's. G2/G3/G4 each independently anchor S2/S3/S4 with one real
    ## sample. Whichever subject G1 is ultimately assigned to, no real or
    ## ghost sample anywhere else legitimately claims it once its one
    ## eligible real sample is used up, so at least two of G1's remaining
    ## samples can only be resolved by inventing placeholder duplicates.
    ## This is the scenario that originally surfaced two samples ending up
    ## with the identical (post-relabel) Sample_ID: one sample's invented
    ## label reused a *different* sample's Placeholder_ID rather than its
    ## own.
    sample_metadata <- data.frame(
        Sample_ID = c("g1a", "g1b", "g1c", "g1d", "g2", "g3", "g4"),
        Subject_ID = c("S1", "S2", "S3", "S4", "S2", "S3", "S4"),
        Genotype_Group_ID = c("G1", "G1", "G1", "G1", "G2", "G3", "G4"),
        stringsAsFactors = FALSE
    )
    x <- MislabelSolver(sample_metadata = sample_metadata)
    solved <- suppressMessages(solveGlobalSearch(x))

    rd <- solved@.solve_state$relabel_data
    unknown_rows <- grepl("LABELNOTFOUND", rd$Sample_ID)
    ## The scenario is only a meaningful test of the invariant if it
    ## actually exercises the "invented duplicate" path.
    expect_true(any(unknown_rows))

    embedded_placeholder_id <- vapply(
        strsplit(rd$Sample_ID[unknown_rows], "#"),
        \(parts) parts[4],
        character(1)
    )
    expect_identical(
        embedded_placeholder_id,
        unname(rd$Placeholder_ID[unknown_rows])
    )
})
