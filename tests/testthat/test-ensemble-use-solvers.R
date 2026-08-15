test_that("global_ghost_penalty/global_deletion_penalty are passed through to solveGlobalSearch()'s validation", {
    # Make a cycle too big for global
    scenario <- toy_cycle_scenario(k = 10)
    n_samples <- nrow(scenario$sample_metadata)
    correct_subjects <- scenario$sample_metadata |>
        select(Sample_ID, Correct_Subject_ID) |>
        tibble::deframe()

    x <- MislabelSolver(sample_metadata = scenario$sample_metadata)

    count_wrong <- function(solver) {
        sample_out <- collateOutput(solver)$Sample
        correct <- sum(
            sample_out$Proposed_Final_Subject_ID ==
                correct_subjects[sample_out$Initial_Sample_ID]
        )
        n_samples - correct
    }

    count_wrong_with_solveEnsemble_args <- function(...) {
        count_wrong(suppressMessages(solveEnsemble(x, ...)))
    }

    # Strategies that should solve successfully (have 0 wrong):
    # Default args
    expect_equal(count_wrong_with_solveEnsemble_args(), 0)
    # Local only
    expect_equal(count_wrong_with_solveEnsemble_args(use_solvers = "local"), 0)
    # Majority only
    expect_equal(
        count_wrong_with_solveEnsemble_args(use_solvers = "majority"),
        0
    )
    # Majority with smaller limit + local
    expect_equal(
        count_wrong_with_solveEnsemble_args(
            use_solvers = c("majority", "local"),
            majority_max_genotypes = 8
        ),
        0
    )

    # Strategies that should fail (have 1 or more wrong):
    # Global only
    expect_gt(count_wrong_with_solveEnsemble_args(use_solvers = "global"), 0)
    # Majority with smaller limit
    expect_gt(
        count_wrong_with_solveEnsemble_args(
            use_solvers = "majority",
            majority_max_genotypes = 8
        ),
        0
    )
})
