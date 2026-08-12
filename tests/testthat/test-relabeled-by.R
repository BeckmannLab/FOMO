test_that("Relabeled_By starts NA for every sample at construction", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
    expect_true(all(is.na(x@.solve_state$relabel_data$Relabeled_By)))
})

test_that("Relabeled_By names the solver that actually relabeled a sample, and stays NA for samples that were never relabeled", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
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
    ## A deterministic scenario (every solver seeds with set.seed(1)
    ## internally) built from a "tie cycle" of k genotype groups <-> k
    ## subjects, each with an exact 1-1 tie (unresolvable by majority
    ## search) -- with k = 10 exceeding global_max_genotypes's default of
    ## 8, the whole component starts out ineligible for global search, so
    ## only local search can make any progress on it at first. Local search
    ## nibbles at it, a couple of samples at a time, across several
    ## outer-loop iterations; only once the remaining unsolved piece
    ## shrinks enough does global search take over and resolve what's left
    ## in one shot. This means some samples end up in a component that
    ## global's final pass marks Solved (Solved_By = "global"), even though
    ## global didn't need to touch every sample in it to do so -- local
    ## already got some of them to their final answer earlier.
    k <- 10
    subj <- sprintf("S%02d", seq_len(k))
    geno <- sprintf("G%02d", seq_len(k))
    a_rows <- data.frame(
        Sample_ID = sprintf("a%02d", seq_len(k)),
        Subject_ID = subj,
        Genotype_Group_ID = geno,
        stringsAsFactors = FALSE
    )
    b_rows <- data.frame(
        Sample_ID = sprintf("b%02d", seq_len(k)),
        Subject_ID = subj[c(2:k, 1)],
        Genotype_Group_ID = geno,
        stringsAsFactors = FALSE
    )
    sample_metadata <- rbind(a_rows, b_rows)
    swap_cats <- data.frame(
        Sample_ID = sample_metadata$Sample_ID,
        SwapCat_ID = "omic1",
        stringsAsFactors = FALSE
    )

    x <- MislabelSolver(
        sample_metadata = sample_metadata,
        swap_cats = swap_cats
    )
    solved <- suppressMessages(solveEnsemble(
        x,
        use_solvers = c("majority", "global", "local"),
        time_limit = 20
    ))
    rd <- solved@.solve_state$relabel_data
    rd <- rd[order(rd$Init_Sample_ID), ]

    ## Sanity check this scenario still fully converges before relying on
    ## its exact, deterministic per-sample resolution below.
    expect_equal(sum(!rd$Solved), 0)

    a05 <- rd[rd$Init_Sample_ID == "a05", ]
    expect_equal(a05$Solved_By, "global")
    expect_equal(a05$Relabeled_By, "local")
    expect_true(a05$Subject_ID != a05$Init_Subject_ID)
})

test_that("Relabeled_By is exposed in collateOutput()'s Sample sheet, immediately after Solved_By", {
    scenario <- toy_swap_scenario()
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        swap_cats = scenario$swap_cats
    )
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
