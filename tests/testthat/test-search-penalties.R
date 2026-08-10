test_that("ghost_penalty and deletion_penalty accept valid scalar numerics", {
    expect_silent(.validate_search_penalties(1.5, 4))
    expect_silent(.validate_search_penalties(2L, 4L)) # integers are bare numerics too
})

test_that("ghost_penalty and deletion_penalty reject non-scalar or non-numeric values", {
    expect_error(.validate_search_penalties(c(1.5, 2), 4), "single \\(length-1\\) numeric")
    expect_error(.validate_search_penalties(1.5, c(4, 5)), "single \\(length-1\\) numeric")
    expect_error(.validate_search_penalties("1.5", 4), "single \\(length-1\\) numeric")
    expect_error(.validate_search_penalties(1.5, NA_real_), "single \\(length-1\\) numeric")
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

test_that("solveGlobalSearch()/solveGlobalSearchFast() default to ghost_penalty=1.5, deletion_penalty=4", {
    expect_equal(formals(solveGlobalSearch)$ghost_penalty, 1.5)
    expect_equal(formals(solveGlobalSearch)$deletion_penalty, 4)
    expect_equal(formals(solveGlobalSearchFast)$ghost_penalty, 1.5)
    expect_equal(formals(solveGlobalSearchFast)$deletion_penalty, 4)
})
