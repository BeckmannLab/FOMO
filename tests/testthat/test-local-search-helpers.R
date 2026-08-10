test_that("local search's inline helpers were extracted as standalone internal functions", {
    expect_true(exists(".calc_scaled_entropy", mode = "function"))
    expect_true(exists(".calc_swapped_delta_entropy", mode = "function"))
})

test_that(".calc_scaled_entropy is 0 for a fully-concordant vote and negative for a split vote", {
    expect_equal(.calc_scaled_entropy(c(A = 3, B = 0)), 0)
    expect_lt(.calc_scaled_entropy(c(A = 2, B = 2)), 0)
})

test_that(".calc_swapped_delta_entropy is 0 for a swap that changes nothing", {
    ## A single-genotype votes matrix where the "swap" moves a vote from
    ## subject A to subject A (a no-op) on the from-side, and the to-side
    ## genotype is NA (i.e. there is no to-side component), so the overall
    ## delta must be exactly 0.
    votes <- matrix(c(3, 0), nrow = 1, dimnames = list("G1", c("A", "B")))
    base_entropies <- c(G1 = .calc_scaled_entropy(votes["G1", ]))
    d <- .calc_swapped_delta_entropy(
        votes, base_entropies,
        swap_from_subject = "A", swap_from_genotype = "G1",
        swap_to_subject = "A", swap_to_genotype = NA
    )
    ## d comes back named "G1" (inherited from indexing base_entropies
    ## internally); compare on value only.
    expect_equal(unname(d), 0)
})

test_that(".calc_swapped_delta_entropy matches a manual before/after entropy calculation", {
    ## Two genotypes, two subjects: move one vote from (G1, A) to (G1, B).
    votes <- matrix(c(3, 1, 1, 2), nrow = 2, dimnames = list(c("G1", "G2"), c("A", "B")))
    base_entropies <- c(
        G1 = .calc_scaled_entropy(votes["G1", ]),
        G2 = .calc_scaled_entropy(votes["G2", ])
    )
    d <- .calc_swapped_delta_entropy(
        votes, base_entropies,
        swap_from_subject = "A", swap_from_genotype = "G1",
        swap_to_subject = "B", swap_to_genotype = NA
    )

    ## d comes back named "G1" (inherited from indexing base_entropies
    ## internally); compare on value only.
    manual_after <- .calc_scaled_entropy(c(A = 2, B = 2))
    expect_equal(unname(d), unname(manual_after - base_entropies[["G1"]]))
})
