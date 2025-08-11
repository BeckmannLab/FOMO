.genotype_matrix_to_genotype_df <- function(genotype_matrix) {
    genotype_graph <- igraph::graph_from_adjacency_matrix(genotype_matrix)
    genotype_group_ids <- igraph::components(genotype_graph)$membership
    n_genotype_groups <- length(unique(genotype_group_ids))
    n_digits <- floor(log10(n_genotype_groups)) + 1
    genotype_group_ids <- vapply(genotype_group_ids, \(x) paste0("Genotype_Group", formatC(x, width = n_digits, format = "d", flag = "0")), "character")
    genotype_df <- data.frame(
        Sample_ID = names(genotype_group_ids),
        Genotype_Group_ID = genotype_group_ids
    )
    return(genotype_df)
}

.generate_corrections_graph <- function(
        relabel_data) {
    sample_corrections_df <- relabel_data |>
        dplyr::filter(Init_Sample_ID != Sample_ID)
    corrections_edges <- sample_corrections_df |>
        dplyr::select(Init_Sample_ID, Sample_ID)

    corrections_vertices <- data.frame(
        Sample_ID = unique(c(sample_corrections_df[, "Init_Sample_ID"],
                             sample_corrections_df[, "Sample_ID"]))
    ) |>
        dplyr::left_join(
            sample_corrections_df |> select(Sample_ID=Init_Sample_ID,
                                             Init_Component_ID,
                                             Component_ID,
                                             Subject_ID,
                                             Genotype_Group_ID,
                                             Is_Ghost,
                                             SwapCat_ID,
                                             SwapCat_Shape,
                                             vertex_size_scalar),
            by = "Sample_ID"
        )
    ## For samples that don't appear in the Init_Sample_ID column (LABELNOTFOUND samples)
    ## need to manually populate fields Is_Ghost, SwapCat_ID, and SwapCat_Shape
    corrections_vertices_split <- corrections_vertices |>
        dplyr::filter(!is.na(Is_Ghost)) |>
        dplyr::mutate(Is_LABELNOTFOUND = FALSE)
    corrections_vertices_label_not_found <- corrections_vertices |>
        dplyr::filter(is.na(Is_Ghost)) |>
        dplyr::mutate(Is_LABELNOTFOUND = TRUE) |>
        dplyr::select("Sample_ID", "Is_LABELNOTFOUND") |>
        dplyr::left_join(
            sample_corrections_df |> select(Sample_ID,
                                             Init_Component_ID,
                                             Component_ID,
                                             Subject_ID,
                                             Genotype_Group_ID,
                                             SwapCat_ID,
                                             SwapCat_Shape,
                                             vertex_size_scalar),
            by = "Sample_ID"
        ) |>
        dplyr::mutate(Is_Ghost=FALSE)
    corrections_vertices <- rbind(corrections_vertices_split,
                                  corrections_vertices_label_not_found) |>
        dplyr::mutate(
            shape = SwapCat_Shape,
            color = case_when(
                Is_LABELNOTFOUND ~ "firebrick",
                Is_Ghost ~ "lightgrey",
                TRUE ~ "orange"
            ),
            size = 12 * vertex_size_scalar,
            label.cex = 0.5
        )

    corrections_graph <- igraph::graph_from_data_frame(corrections_edges, vertices=corrections_vertices, directed=TRUE)
    igraph::E(corrections_graph)$color <- "black"
    igraph::E(corrections_graph)$width <- 6

    return(corrections_graph)
}

# .diagnose_contaminated_genotype_groups <- function(genotype_matrix) {
#
# }

