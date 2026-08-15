test_that("Relabeled_By starts NA for every sample at construction", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)
    expect_true(all(is.na(x@.solve_state$relabel_data$Relabeled_By)))
})

test_that("Relabeled_By names the solver that actually relabeled a sample, and stays NA for samples that were never relabeled", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)
    solved <- suppressMessages(solveEnsemble(x))
    rd <- solved@.solve_state$relabel_data
    rd <- rd[order(rd$Init_Sample_ID), ]

    ## sample1 was already correctly labeled S1 and is never touched by the
    ## swap that resolves this scenario (sample2 <-> sample3).
    expect_true(is.na(rd$Relabeled_By[rd$Init_Sample_ID == "sample1"]))
    expect_equal(rd$Relabeled_By[rd$Init_Sample_ID == "sample2"], "global")
    expect_equal(rd$Relabeled_By[rd$Init_Sample_ID == "sample3"], "global")
})

test_that("Relabeled_By can differ from Solved_By: it tracks the last solver to actually touch a sample's own label, not just whichever solver's pass resolved its component", {
    ## A deterministic scenario built from a cyclic mislabel involving
    ## more than k genotype groups where with k = 10 exceeding
    ## global_max_genotypes's default of 8, the whole component starts
    ## out ineligible for global search, and we exclude majority
    ## search so only local search can make any progress on it at
    ## first. Local search nibbles at it, a couple of samples at a
    ## time; only once the remaining unsolved piece shrinks enough
    ## does global search take over and resolve what's left in one
    ## shot.
    k <- 10
    subj <- sprintf("S%02d", seq_len(k))
    geno <- sprintf("G%02d", seq_len(k))
    a_rows <- tibble(
        Sample_ID = sprintf("a%02d", seq_len(k)),
        Subject_ID = subj,
        Correct_Subject_ID = subj,
        Actually_Mislabeled = FALSE,
        Genotype_Group_ID = geno
    )
    b_rows <- tibble(
        Sample_ID = sprintf("b%02d", seq_len(k)),
        Subject_ID = subj[c(2:k, 1)],
        Correct_Subject_ID = subj,
        Actually_Mislabeled = TRUE,
        Genotype_Group_ID = geno
    )
    c_rows <- tibble(
        Sample_ID = sprintf("c%02d", seq_len(k)),
        Subject_ID = subj,
        Correct_Subject_ID = subj,
        Actually_Mislabeled = FALSE,
        Genotype_Group_ID = geno
    )
    sample_metadata <- rbind(a_rows, b_rows, c_rows)

    x <- MislabelSolver(sample_metadata)
    solved <- suppressMessages(solveEnsemble(
        x,
        use_solvers = c("global", "local"),
        global_max_genotypes = 8,
        time_limit = 20
    ))

    out <- collateOutput(solved)

    res <- out$Sample |>
        select(
            Initial_Subject_ID,
            Initial_Sample_ID,
            Genotype_Group_ID,
            Proposed_Final_Subject_ID,
            Proposed_Final_Sample_ID,
            Mislabeled,
            Solved_By,
            Relabeled_By
        ) |>
        inner_join(
            sample_metadata,
            dplyr::join_by(
                Initial_Subject_ID == Subject_ID,
                Initial_Sample_ID == Sample_ID,
                Genotype_Group_ID
            )
        )
    with(res, table(Mislabeled, Actually_Mislabeled, useNA = "ifany"))

    solved_table <- table(res$Solved_By, useNA = "ifany")
    relabeled_table <- table(res$Relabeled_By, useNA = "ifany")

    ## Sanity check this scenario still fully converges to the correct answer
    expect(!anyNA(solved_table), "All samples are solved")
    expect_equal(res$Mislabeled, res$Actually_Mislabeled, ignore_attr = TRUE)
    expect_equal(
        res$Proposed_Final_Subject_ID,
        res$Correct_Subject_ID,
        ignore_attr = TRUE
    )

    expect(
        all(c("local", "global") %in% res$Solved_By),
        "Both global and local solvers solved at least one component"
    )
    expect(
        all(c("local", "global") %in% res$Relabeled_By),
        "Both global and local solvers solved at least one sample"
    )
    expect_equal(relabeled_table["global"], 8, ignore_attr = TRUE)
    expect_equal(relabeled_table["local"], 2, ignore_attr = TRUE)
})

test_that("Relabeled_By is exposed in collateOutput()'s Sample sheet, immediately after Solved_By", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)
    solved <- suppressMessages(solveEnsemble(x))
    out <- collateOutput(solved)

    col_names <- colnames(out$Sample)
    expect_true("Relabeled_By" %in% col_names)
    expect_equal(
        which(col_names == "Relabeled_By"),
        which(col_names == "Solved_By") + 1
    )

    sample_row <- out$Sample[out$Sample$Initial_Sample_ID == "sample2", ]
    expect_equal(sample_row$Relabeled_By, "global")
})
