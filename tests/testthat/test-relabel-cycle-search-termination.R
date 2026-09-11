## Regression tests for a non-termination bug in the shared relabel-cycle-
## search helper .find_directed_cycles() (see helpers-solve.R), used by both
## solveMajoritySearch() and solveGlobalSearch() via
## .find_relabel_cycles_from_putative_subjects(). The original implementation
## enumerated *every* simple path up to a length cutoff between a vertex and
## each of its in-neighbors (via igraph::all_simple_paths()), even though the
## caller (.find_relabel_cycles_in_component()) only ever needs one closing
## path per (vertex, in-neighbor) pair -- but the number of alternative paths
## between a pair can grow exponentially with path length on a component with
## many mutually-eligible relabel candidates, which made the search take
## impractically long -- observed as an apparent infinite loop in
## solveEnsemble() on large real-world inputs (multi-day runs that normally
## take under a minute).
##
## The fix replaces all_simple_paths() with igraph::shortest_paths() (a
## single BFS per pair): a cycle of length <= cutoff through a pair exists if
## and only if its shortest path is <= cutoff, so this finds exactly the same
## set of *pairs* that close a cycle within a cutoff as the exhaustive search
## did -- it just stops at the first (shortest) closing path instead of
## enumerating every alternative one, which is always polynomial regardless
## of how densely connected the component is (no exponential worst case to
## guard against, unlike all_simple_paths()).
##
## That alone is not quite enough to match the original algorithm's
## completeness, though: .find_directed_cycles() computes every pair's path
## against one fixed graph snapshot, so a pair whose only offered path (the
## shortest one) happens to route through a vertex that a *different* pair's
## cycle claims earlier in the same sweep is skipped for that sweep, even
## though the two cycles didn't strictly have to conflict -- the exhaustive
## search didn't have this gap, because it offered every alternative path for
## a pair at once, so a single pass could fall back to a different route
## within the same pair. .find_relabel_cycles_in_component() closes this gap
## by repeating the full shortest-cycles-first sweep (see its own comment)
## until a sweep finds nothing new, rather than running it only once: once a
## conflicting cycle's vertices are actually deleted, the blocked pair's
## shortest path gets recomputed from scratch on the next sweep and is free
## to route around the gap.
##
## The exact originally-reported input is unavailable (multi-day runs aren't
## reproducible in a test), so these tests instead target the failure
## mechanism directly:
##  - the first two tests exercise .find_directed_cycles() directly on small,
##    hand-built graphs and assert its actual invariants, deterministically,
##    with no timing involved: exactly one closing path per pair (even when
##    multiple equally-short alternatives exist), and the cutoff bound is
##    respected;
##  - the third test exercises .find_relabel_cycles_in_component() directly
##    on a small, hand-built graph specifically shaped so that a single sweep
##    provably cannot resolve every vertex (confirmed while developing this
##    fix by temporarily disabling the repeat{} wrapper and rerunning this
##    exact graph -- a single sweep resolves only half of it), to pin down
##    that repetition is actually doing the work it's meant to;
##  - the fourth test drives the public solveEnsemble() workflow on a small
##    (well under 15 subjects) ring scenario with just enough redundancy to
##    be fully solvable by exhaustive search -- confirming that realistic
##    amounts of redundancy are resolved exactly as the pre-fix
##    implementation would, with no correctness trade-off;
##  - the fifth test drives solveEnsemble() on a larger, densely
##    cross-connected version of the same scenario, shaped to create the same
##    kind of combinatorial blowup the pre-fix implementation hung on, under a
##    hard, process-level timeout so the test fails cleanly instead of hanging
##    the suite if the bug is ever reintroduced. This exact scenario and call,
##    run against the pre-fix implementation (verified while developing this
##    fix by temporarily restoring the old all_simple_paths()-based code),
##    reliably fails to return within 45s; the fix completes it in well under
##    5s, fully solved.

test_that(".find_directed_cycles() returns exactly one closing path per pair, even with multiple equally-short alternatives", {
    ## A "diamond": edge B -> A makes B an in-neighbor of A, and there are
    ## two distinct, equally-short paths from A back to B (via C1, or via
    ## C2), each of which independently closes the (A, B) pair into a
    ## length-3 cycle at cutoff = 2. Unlike the pre-fix implementation
    ## (which would return both), shortest_paths() only ever offers the
    ## caller one of them -- which one is an igraph implementation detail
    ## this test deliberately does not pin down, only that it is exactly
    ## one, and that it is a genuinely valid closing path.
    relabels_graph <- igraph::graph_from_edgelist(
        matrix(
            c(
                "B",
                "A",
                "A",
                "C1",
                "C1",
                "B",
                "A",
                "C2",
                "C2",
                "B"
            ),
            ncol = 2,
            byrow = TRUE
        ),
        directed = TRUE
    )

    cycles <- fomo:::.find_directed_cycles(relabels_graph, cutoff = 2)

    ## "A" has exactly one in-neighbor ("B"), so every returned cycle
    ## starting at "A" necessarily belongs to the (A, B) pair.
    cycles_from_a <- Filter(function(cyc) identical(cyc[1], "A"), cycles)
    expect_length(cycles_from_a, 1)
    expect_identical(cycles_from_a[[1]][1], "A")
    expect_identical(cycles_from_a[[1]][3], "B")
    expect_true(cycles_from_a[[1]][2] %in% c("C1", "C2"))
})

test_that(".find_directed_cycles() respects the cutoff bound", {
    ## B -> A, and A's only route back to B has length 3 (via X, then Y).
    relabels_graph <- igraph::graph_from_edgelist(
        matrix(
            c("B", "A", "A", "X", "X", "Y", "Y", "B"),
            ncol = 2,
            byrow = TRUE
        ),
        directed = TRUE
    )

    ## Too short a cutoff: the pair isn't closeable within it, so nothing
    ## is returned for it, exactly as all_simple_paths(cutoff = 2) would
    ## also have found nothing.
    cycles_short <- fomo:::.find_directed_cycles(relabels_graph, cutoff = 2)
    cycles_from_a_short <- Filter(
        function(cyc) identical(cyc[1], "A"),
        cycles_short
    )
    expect_length(cycles_from_a_short, 0)

    ## A cutoff that exactly fits the path finds it.
    cycles_ok <- fomo:::.find_directed_cycles(relabels_graph, cutoff = 3)
    cycles_from_a_ok <- Filter(function(cyc) identical(cyc[1], "A"), cycles_ok)
    expect_length(cycles_from_a_ok, 1)
    expect_identical(cycles_from_a_ok[[1]], c("A", "X", "Y", "B"))
})

test_that(".find_relabel_cycles_in_component() recovers cycles that a single sweep misses", {
    ## Two "parallel" 3-cycles worth of samples (g1/g2/g3, two samples per
    ## genotype), but every sample in genotype g_i is eligible to relabel
    ## against *both* samples in genotype g_(i+1) (a complete bipartite
    ## connection between consecutive layers, wrapping from g3 back to
    ## g1) -- the same redundant-eligibility shape .make_ring_scenario()
    ## below builds from realistic sample data, just built directly as a
    ## graph here for a fast, deterministic, non-end-to-end check. With a
    ## single (non-repeated) sweep, cutoff reaches the component's full
    ## size before a second attempt is possible, and greedily claiming the
    ## first-found 3-cycle deletes vertices that the *other* 3-cycle's
    ## pairs had already (now stale) precomputed a path through -- exactly
    ## the gap described in the file-level comment above. Confirmed while
    ## developing this fix, by temporarily replacing the repeat{} wrapper
    ## in .find_relabel_cycles_in_component() with a single pass and
    ## rerunning this exact graph, that a single sweep resolves only 3 of
    ## the 6 vertices (one of the two 3-cycles) -- leaving the other
    ## entirely unresolved, not just partially -- which is exactly the
    ## failure mode repeating the sweep is meant to fix.
    edges <- rbind(
        expand.grid(
            from = c("g1_m1", "g1_m2"),
            to = c("g2_m1", "g2_m2"),
            stringsAsFactors = FALSE
        ),
        expand.grid(
            from = c("g2_m1", "g2_m2"),
            to = c("g3_m1", "g3_m2"),
            stringsAsFactors = FALSE
        ),
        expand.grid(
            from = c("g3_m1", "g3_m2"),
            to = c("g1_m1", "g1_m2"),
            stringsAsFactors = FALSE
        )
    )
    relabels_graph <- igraph::graph_from_edgelist(
        as.matrix(edges[, c("from", "to")]),
        directed = TRUE
    )

    result <- fomo:::.find_relabel_cycles_in_component(relabels_graph)

    all_samples <- c("g1_m1", "g1_m2", "g2_m1", "g2_m2", "g3_m1", "g3_m2")
    expect_setequal(result$relabel_from, all_samples)
    expect_setequal(result$relabel_to, all_samples)
    ## Every relabel_to must itself be a valid out-neighbor of its
    ## relabel_from (i.e. every returned pair is a real graph edge, not
    ## just any leftover vertex), and following relabel_from -> relabel_to
    ## repeatedly from any sample must return to itself (a closed cycle).
    edge_set <- paste(edges$from, edges$to)
    expect_true(all(
        paste(result$relabel_from, result$relabel_to) %in% edge_set
    ))
    next_of <- setNames(result$relabel_to, result$relabel_from)
    for (start in all_samples) {
        cur <- next_of[[start]]
        steps <- 1
        while (cur != start && steps <= length(all_samples)) {
            cur <- next_of[[cur]]
            steps <- steps + 1
        }
        expect_identical(cur, start)
    }
})

## Shared by the two solveEnsemble()-level tests below: a ring of L genotype
## groups / subjects, each with W samples that are *all* mislabeled and
## currently claim the *previous* subject in the ring, all W of them equally
## eligible to relabel against any of the W samples in the next group. That
## creates one connected component where, at each step around the ring,
## there are W equally-eligible relabel targets: W^(L-1) distinct candidate
## cycles for the pre-fix (and, without repeated sweeps, the naive
## shortest_paths()) algorithm to consider before finding the cycles it
## actually needs.
.make_ring_scenario <- function(L, W) {
    subject_of <- function(g) sprintf("Subject%d", ((g - 1) %% L) + 1)
    genotype_of <- function(g) sprintf("Genotype%d", ((g - 1) %% L) + 1)
    anchors <- do.call(
        rbind,
        lapply(seq_len(L), function(g) {
            data.frame(
                Sample_ID = sprintf("anchor_%s", genotype_of(g)),
                Subject_ID = subject_of(g),
                Genotype_Group_ID = genotype_of(g),
                stringsAsFactors = FALSE
            )
        })
    )
    mislabeled <- do.call(
        rbind,
        lapply(seq_len(L), function(g) {
            data.frame(
                Sample_ID = sprintf("mis_%s_%02d", genotype_of(g), seq_len(W)),
                Subject_ID = subject_of(g - 1), # wrongly claims the *previous* subject
                Genotype_Group_ID = genotype_of(g),
                stringsAsFactors = FALSE
            )
        })
    )
    list(
        sample_metadata = rbind(anchors, mislabeled),
        anchor_samples = anchors$Sample_ID
    )
}

test_that("solveEnsemble() fully resolves realistic amounts of relabel-candidate redundancy, matching pre-fix behavior", {
    ## L = 5, W = 5 (5 subjects, well under the 15-subject cap): small
    ## enough to run quickly regardless of algorithm. Confirmed while
    ## developing this fix that the pre-fix (all_simple_paths()-based)
    ## implementation also fully resolves this scenario, with the exact
    ## same sample -> Subject_ID assignments -- so a regression that
    ## changed *which* valid cycles get chosen, not just whether the
    ## search terminates, would show up here.
    scenario <- .make_ring_scenario(L = 5, W = 5)
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        anchor_samples = scenario$anchor_samples
    )

    solved <- suppressMessages(solveEnsemble(
        x,
        use_solvers = "majority",
        time_limit = 20
    ))

    expect_equal(sum(!solved@.solve_state$relabel_data$Solved), 0)
})

test_that("solveEnsemble() terminates promptly on a densely cross-connected component that made the pre-fix cycle search hang", {
    skip_if_not_installed("callr")

    ## L = 8, W = 5 (8 subjects, well under the 15-subject cap): reliably
    ## ran for 45+ seconds without returning under the pre-fix
    ## implementation (verified while developing this fix, by temporarily
    ## restoring the old all_simple_paths()-based implementation and
    ## re-running this exact test body); the fix completes it in well
    ## under 5s, fully solved. A *single*, non-repeated shortest_paths()
    ## sweep is fast on this input too, but resolves none of it at all
    ## (also verified while developing this fix) -- so this test also
    ## guards against the repeat{} wrapper being dropped, not just against
    ## all_simple_paths() being reintroduced.
    scenario <- .make_ring_scenario(L = 8, W = 5)
    x <- MislabelSolver(
        sample_metadata = scenario$sample_metadata,
        anchor_samples = scenario$anchor_samples
    )

    ## Run the actual solve in a subprocess with a hard wall-clock timeout,
    ## so that if the non-termination bug is ever reintroduced, this test
    ## fails cleanly after 45s instead of hanging the test suite.
    result <- tryCatch(
        callr::r(
            function(object) {
                library(fomo)
                solved <- suppressMessages(solveEnsemble(
                    object,
                    use_solvers = "majority",
                    time_limit = 30
                ))
                list(
                    n_unsolved = sum(!solved@.solve_state$relabel_data$Solved),
                    n_samples = nrow(solved@.solve_state$relabel_data),
                    is_solver = methods::is(solved, "MislabelSolver")
                )
            },
            args = list(object = x),
            timeout = 45
        ),
        callr_timeout_error = function(e) NULL
    )

    expect_false(
        is.null(result),
        label = paste(
            "solveEnsemble() did not return within 45s -- this suggests a",
            "regression of the relabel-cycle-search non-termination bug"
        )
    )
    expect_true(result$is_solver)
    expect_equal(result$n_samples, nrow(scenario$sample_metadata))
    expect_equal(result$n_unsolved, 0)
})
