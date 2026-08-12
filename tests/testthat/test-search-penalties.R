test_that("ghost_penalty and deletion_penalty accept valid scalar numerics", {
    expect_silent(.validate_search_penalties(1.5, 4))
    expect_silent(.validate_search_penalties(2L, 4L)) # integers are bare numerics too
})

test_that("ghost_penalty and deletion_penalty reject non-scalar or non-numeric values", {
    expect_error(
        .validate_search_penalties(c(1.5, 2), 4),
        "single \\(length-1\\) numeric"
    )
    expect_error(
        .validate_search_penalties(1.5, c(4, 5)),
        "single \\(length-1\\) numeric"
    )
    expect_error(
        .validate_search_penalties("1.5", 4),
        "single \\(length-1\\) numeric"
    )
    expect_error(
        .validate_search_penalties(1.5, NA_real_),
        "single \\(length-1\\) numeric"
    )
})

test_that("non-positive penalties are always rejected (even if otherwise valid scalars)", {
    expect_error(.validate_search_penalties(0, 4), "positive")
    expect_error(.validate_search_penalties(-1, 4), "positive")
    expect_error(.validate_search_penalties(1.5, 0), "positive")
    expect_error(.validate_search_penalties(1.5, -2), "positive")
})

test_that("recommended-limit violations warn but do not error", {
    expect_warning(.validate_search_penalties(1, 4), "ghost_penalty")
    expect_warning(.validate_search_penalties(0.5, 4), "ghost_penalty")
    expect_warning(.validate_search_penalties(1.5, 1), "deletion_penalty")
    expect_warning(.validate_search_penalties(1.5, 1.9), "deletion_penalty")
})

test_that("deletion_penalty must be strictly greater than twice the relabel penalty, not merely equal", {
    ## At exactly 2x the relabel penalty, a plain two-sample swap (see
    ## toy_swap_scenario()) and treating the pair as one deletion score
    ## identically -- an exact tie decided only by incidental permutation
    ## order, not a real guarantee the honest fix wins. That makes
    ## deletion_penalty == 2 * relabel_penalty just as warning-worthy as
    ## deletion_penalty < 2 * relabel_penalty.
    expect_warning(.validate_search_penalties(1.5, 2), "deletion_penalty")
    expect_silent(.validate_search_penalties(1.5, 2.01))
})

test_that("deletion_penalty must also be strictly greater than ghost_penalty", {
    ## Clearing "> 2 * relabel_penalty" is not sufficient on its own: if
    ## deletion_penalty <= ghost_penalty, the algorithm can still prefer
    ## deleting a sample over relabeling it to an available ghost (a second,
    ## independent failure mode from the plain-swap-vs-deletion one above).
    expect_warning(.validate_search_penalties(3, 2.5), "deletion_penalty") # 2.5 > 2, but <= ghost_penalty (3)
    expect_warning(.validate_search_penalties(3, 3), "deletion_penalty") # exact tie with ghost_penalty
    expect_silent(.validate_search_penalties(3, 3.01))
})

test_that("solveGlobalSearch() defaults to ghost_penalty=1.5, deletion_penalty=4", {
    expect_equal(formals(solveGlobalSearch)$ghost_penalty, 1.5)
    expect_equal(formals(solveGlobalSearch)$deletion_penalty, 4)
})
