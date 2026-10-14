test_that("global search chooses randomly from among tied permutations", {
    sample_meta_data <- data.frame(
        Sample_ID = c(
            "Subj2_Sample1",
            "Subj1_Sample2",
            "Subj1_Sample1",
            "Subj2_Sample2"
        ),
        Subject_ID = c("Subj2", "Subj1", "Subj1", "Subj2"),
        Sex = c("M", "F", "F", "M"),
        True_Sample_ID = c(
            "Subj1_Sample1",
            "Subj1_Sample2",
            "Subj2_Sample1",
            "Subj2_Sample2"
        ),
        True_Subject_ID = c("Subj1", "Subj1", "Subj2", "Subj2"),
        Genotype_Group_ID = c(
            "Genotype_Group1",
            "Genotype_Group1",
            "Genotype_Group2",
            "Genotype_Group2"
        ),
        Mislabeled = c(TRUE, FALSE, TRUE, FALSE),
        Mislabeled_Subject = c(TRUE, FALSE, TRUE, FALSE)
    )

    correct <- map_lgl(seq_len(10), function(seed) {
        solver <- MislabelSolver(sample_meta_data)
        res <- suppressMessages(
            collateOutput(solveEnsemble(solver, seed = seed))$Sample
        )
        true_ids <- sample_meta_data %>%
            select(Sample_ID, True_Subject_ID) %>%
            tibble::deframe()
        solved_ids <- res %>%
            select(Initial_Sample_ID, Proposed_Final_Subject_ID) %>%
            tibble::deframe()
        all(true_ids == solved_ids[names(true_ids)])
    })

    # We expect to get things wrong some of the time and right some of the time.
    # All of one or the other indicates a systematic bias in selecting
    # permutations.
    expect_in(c(TRUE, FALSE), correct)
})
