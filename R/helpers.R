tsmsg <- function(...) {
    message(date(), ": ", ...)
}

## Deterministically turn an arbitrary R object into a non-negative integer
## seed, via a hash of its serialized content. Used to seed a dedicated RNG
## stream from a MislabelSolver()'s (sorted) input -- see
## .generate_placeholder_ids() -- so that randomly-generated placeholder IDs
## are reproducible for a given dataset without being fixed across all
## datasets. Truncated to 7 hex digits (a 28-bit non-negative integer) to
## stay safely within set.seed()'s accepted range on any platform.
.hash_to_seed <- function(x) {
    strtoi(str_sub(digest(x, algo = "md5"), 1, 7), base = 16L)
}

## Generate n random, UUID-v4-*shaped* placeholder IDs using R's own seeded
## RNG (sample()). Deliberately NOT uuid::UUIDgenerate(): that function draws
## from system entropy and does not respond to set.seed()/with_seed()
## at all (confirmed empirically), so results from it can never be made
## reproducible that way. The UUID v4 shape (8-4-4-4-12 hex digits, version
## nibble fixed at "4", variant nibble drawn from {8,9,a,b}) is purely
## cosmetic, kept only so these placeholder IDs look like the
## uuid::UUIDgenerate() output they replace.
.generate_placeholder_ids <- function(n) {
    hex_digits <- c(0:9, letters[1:6])
    hex <- function(k) {
        str_c(sample(hex_digits, k, replace = TRUE), collapse = "")
    }
    vapply(
        seq_len(n),
        function(i) {
            str_c(
                hex(8),
                hex(4),
                str_c("4", hex(3)),
                str_c(sample(c("8", "9", "a", "b"), 1), hex(3)),
                hex(12),
                sep = "-"
            )
        },
        character(1)
    )
}

.genotype_matrix_to_genotype_df <- function(genotype_matrix) {
    genotype_graph <- graph_from_adjacency_matrix(genotype_matrix)
    genotype_group_ids <- components(genotype_graph)$membership
    n_genotype_groups <- length(unique(genotype_group_ids))
    n_digits <- floor(log10(n_genotype_groups)) + 1
    genotype_group_ids <- vapply(
        genotype_group_ids,
        \(x) {
            str_c(
                "Genotype_Group",
                formatC(x, width = n_digits, format = "d", flag = "0")
            )
        },
        "character"
    )
    genotype_df <- data.frame(
        Sample_ID = names(genotype_group_ids),
        Genotype_Group_ID = genotype_group_ids
    )
    return(genotype_df)
}

.generate_corrections_graph <- function(
    relabel_data
) {
    sample_corrections_df <- relabel_data |>
        filter(.data$Init_Sample_ID != .data$Sample_ID)
    corrections_edges <- sample_corrections_df |>
        select("Init_Sample_ID", "Sample_ID")

    corrections_vertices <- data.frame(
        Sample_ID = unique(c(
            sample_corrections_df[, "Init_Sample_ID"],
            sample_corrections_df[, "Sample_ID"]
        ))
    ) |>
        left_join(
            sample_corrections_df |>
                select(
                    Sample_ID = "Init_Sample_ID",
                    "Init_Component_ID",
                    "Component_ID",
                    "Subject_ID",
                    "Genotype_Group_ID",
                    "Is_Ghost",
                    "Label_Domain",
                    "Label_Domain_Shape",
                    "vertex_size_scalar"
                ),
            by = "Sample_ID"
        )
    ## For samples that don't appear in the Init_Sample_ID column (LABELNOTFOUND samples)
    ## need to manually populate fields Is_Ghost, Label_Domain, and Label_Domain_Shape
    corrections_vertices_split <- corrections_vertices |>
        filter(!is.na(.data$Is_Ghost)) |>
        mutate(Is_LABELNOTFOUND = FALSE)
    corrections_vertices_label_not_found <- corrections_vertices |>
        filter(is.na(.data$Is_Ghost)) |>
        mutate(Is_LABELNOTFOUND = TRUE) |>
        select("Sample_ID", "Is_LABELNOTFOUND") |>
        left_join(
            sample_corrections_df |>
                select(
                    "Sample_ID",
                    "Init_Component_ID",
                    "Component_ID",
                    "Subject_ID",
                    "Genotype_Group_ID",
                    "Label_Domain",
                    "Label_Domain_Shape",
                    "vertex_size_scalar"
                ),
            by = "Sample_ID"
        ) |>
        mutate(Is_Ghost = FALSE)
    corrections_vertices <- rbind(
        corrections_vertices_split,
        corrections_vertices_label_not_found
    ) |>
        mutate(
            shape = .data$Label_Domain_Shape,
            color = case_when(
                .data$Is_LABELNOTFOUND ~ PLOT_COLOR_LABEL_NOT_FOUND,
                .data$Is_Ghost ~ PLOT_COLOR_GHOST,
                TRUE ~ PLOT_COLOR_REGULAR_SAMPLE
            ),
            size = 12 * .data$vertex_size_scalar,
            label.cex = 0.5
        )

    corrections_graph <- graph_from_data_frame(
        corrections_edges,
        vertices = corrections_vertices,
        directed = TRUE
    )
    E(corrections_graph)$color <- PLOT_COLOR_DEFAULT_EDGE
    E(corrections_graph)$width <- 6

    return(corrections_graph)
}

# Internal helper to rename the components of old MislabelSolver objects so that
# they work with new package versions.
fixup_MislabelSolver <- function(object) {
    if (is.null(object@label_domains) && !is.null(object@swap_cats)) {
        object@label_domains <- object@swap_cats |>
            rename(
                Label_Domain = .data$SwapCat_ID,
                Label_Domain_Shape = .data$SwapCat_Shape
            )
        object@swap_cats <- NULL
    }
    object
}
