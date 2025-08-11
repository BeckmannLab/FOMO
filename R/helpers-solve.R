#' @import dplyr
#' @import igraph
#'
.update_solve_state <- function(object, initialization=FALSE) {
    if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
        return(object)
    }

    ## 1. Assign a Component_ID for each unsolved Sample_ID
    combined_graph <- .generate_graph(object@.solve_state$unsolved_relabel_data,
                                      graph_type="combined", object@.solve_state$unsolved_ghost_data)
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data %>%
        dplyr::mutate(
            Component_ID = as.character(igraph::components(combined_graph)$membership[Sample_ID])
        )
    unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data %>%
        dplyr::mutate(
            Component_ID = as.character(igraph::components(combined_graph)$membership[Sample_ID])
        )

    ## 2. Determine which Sample_ID(s) are newly solved
    ##    A Sample_ID is solved if it belongs to a Component_ID that includes only one
    ##    Genotype_Group_ID and one Subject_ID
    component_data <- rbind(unsolved_relabel_data, unsolved_ghost_data) %>%
        dplyr::group_by(Component_ID) %>%
        dplyr::summarise(
            n_Genotype_Group_ID = n_distinct(Genotype_Group_ID) - anyNA(Genotype_Group_ID),
            n_Subject_ID = dplyr::n_distinct(Subject_ID),
            n_Sample_ID = length(Sample_ID)
        ) %>%
        dplyr::mutate(Solved = n_Genotype_Group_ID <= 1 & n_Subject_ID == 1)

    solved_components <- component_data %>%
        dplyr::filter(Solved) %>%
        dplyr::pull(Component_ID)

    unsolved_relabel_data <- unsolved_relabel_data %>%
        dplyr::mutate(Solved = Component_ID %in% solved_components)
    unsolved_ghost_data <- unsolved_ghost_data %>%
        dplyr::mutate(Solved = Component_ID %in% solved_components)

    ## 3. Re-rank Component_ID(s) in order of size (so that Component1 is the largest unsolved component)
    # component_data <- component_data %>% filter(!Solved)
    n_components <- nrow(component_data)
    component_data <- component_data %>%
        dplyr::arrange(Solved, dplyr::desc(n_Sample_ID)) %>%
        dplyr::mutate(
            new_Component_ID = if (n_components != 0) seq_len(n_components) else character(0),
            new_Component_ID = paste0("Component_", formatC(new_Component_ID, width=nchar(n_components), format="d", flag="0"))
        )
    unsolved_relabel_data <- unsolved_relabel_data %>%
        dplyr::left_join(component_data[, c("Component_ID", "new_Component_ID")], by="Component_ID") %>%
        dplyr::mutate(Component_ID = new_Component_ID) %>%
        dplyr::select(-new_Component_ID)
    unsolved_ghost_data <- unsolved_ghost_data %>%
        dplyr::left_join(component_data[, c("Component_ID", "new_Component_ID")], by="Component_ID") %>%
        dplyr::mutate(Component_ID = new_Component_ID) %>%
        dplyr::select(-new_Component_ID)

    ## 4. Update putative_subjects
    ##    During initialization, lock Subject_ID/Genotype_Group_ID pairs for anchor_samples
    if (initialization) {
        anchor_putative_subjects <- object@sample_metadata %>%
            dplyr::filter(
                !is.na(Genotype_Group_ID),
                !is.na(Subject_ID),
                Sample_ID %in% object@anchor_samples) %>%
            dplyr::select(Genotype_Group_ID, Subject_ID) %>%
            dplyr::distinct()
        object <- .update_putative_subjects(object, anchor_putative_subjects)
    }
    solved_putative_subjects <- unsolved_relabel_data %>%
        dplyr::filter(Solved) %>%
        dplyr::select(Genotype_Group_ID, Subject_ID) %>%
        dplyr::distinct()
    object <- .update_putative_subjects(object, solved_putative_subjects)

    ## 5. Update relabel_data, and unsolved_relabel_data
    if (initialization) {
        col_order <- c("Init_Sample_ID", "Init_Subject_ID", "Genotype_Group_ID", "Component_ID",
                       "Sample_ID", "Subject_ID", "Solved", "Is_Ghost", "Is_Anchor", "SwapCat_ID",
                       "SwapCat_Shape", "vertex_size_scalar")
        unsolved_relabel_data <- unsolved_relabel_data %>%
            dplyr::select(all_of(col_order)) %>%
            dplyr::mutate(Init_Component_ID = Component_ID) %>%
            dplyr::relocate(Init_Component_ID, .before=Component_ID)
        unsolved_ghost_data <- unsolved_ghost_data %>%
            dplyr::select(all_of(col_order)) %>%
            dplyr::mutate(Init_Component_ID = Component_ID) %>%
            dplyr::relocate(Init_Component_ID, .before=Component_ID)
        relabel_data <- rbind(unsolved_relabel_data, unsolved_ghost_data)
    } else {
        ## Update 'relabel_data' with new sample labels in 'unsolved_relabel_data' and 'unsolved_ghost_data'
        unsolved_data <- rbind(unsolved_relabel_data, unsolved_ghost_data)
        relabel_data <- dplyr::rows_update(object@.solve_state$relabel_data, unsolved_data, by="Init_Sample_ID")
    }
    unsolved_relabel_data <- unsolved_relabel_data %>% dplyr::filter(!Solved)
    unsolved_ghost_data <- unsolved_ghost_data %>% dplyr::filter(!Solved)

    ## 6. Overwrite .solve_state
    object@.solve_state$relabel_data <- relabel_data
    object@.solve_state$unsolved_relabel_data <- unsolved_relabel_data
    object@.solve_state$unsolved_ghost_data <- unsolved_ghost_data

    return(object)
}

#' @import dplyr
#' @importFrom stats na.omit
#'
.update_putative_subjects <- function(object, proposed_putative_subjects) {
    Subject_ID <- Genotype_Group_ID <- NULL

    if (nrow(proposed_putative_subjects) == 0) return(object)
    ## Only add Genotype_Group_ID/Subject_ID combinations if neither the
    ## Genotype_Group_ID nor the Subject_ID are already in putative_subjects
    existing_genotypes <- stats::na.omit(object@.solve_state$putative_subjects$Genotype_Group_ID)
    existing_subjects <- stats::na.omit(object@.solve_state$putative_subjects$Subject_ID)
    proposed_putative_subjects <- proposed_putative_subjects %>%
        dplyr::filter(
            !(Subject_ID %in% existing_subjects),
            !(Genotype_Group_ID %in% existing_genotypes)
        )
    putative_subjects <- rbind(object@.solve_state$putative_subjects, proposed_putative_subjects)
    # .validate_putative_subjects(object@sample_genotype_data, putative_subjects)
    object@.solve_state$putative_subjects <- putative_subjects
    return(object)
}

#' @import dplyr
#' @import igraph
#' @import stringr
#' @import tidyr
#'
.generate_graph <- function(
        relabel_data,
        graph_type=c("label", "genotype", "combined"),
        ghost_data=NULL,
        genotype_matrix=NULL,
        populate_plotting_attributes=FALSE,
        collapse_samples=FALSE
) {
    graph_type_mapping <- list(
        label = "Subject_ID",
        genotype = "Genotype_Group_ID",
        combined = NA_character_
    )

    ## Collapse samples mode should only run if the graph is being plotted
    if (!populate_plotting_attributes) {
        collapse_samples <- FALSE
    }

    graph_type <- as.character(graph_type)
    graph_type <- match.arg(graph_type)

    ## Ignore ghost samples when constructing genotype graph
    if (graph_type == "genotype") {
        ghost_data <- NULL
    }

    if (collapse_samples) {
        relabel_data <- relabel_data %>%
            dplyr::group_by(Subject_ID, Genotype_Group_ID, SwapCat_ID) %>%
            dplyr::mutate(
                count = n(),
                vertex_size_scalar = sqrt(sum(vertex_size_scalar)),
                Is_Ghost = FALSE,
                Is_Anchor = any(Is_Anchor),
                Sample_ID = ifelse(count == 1, Sample_ID, paste(paste(count, "samples"), Subject_ID, SwapCat_ID, sep="\n"))
            ) %>%
            dplyr::select(Sample_ID, Subject_ID, Genotype_Group_ID, SwapCat_ID, SwapCat_Shape, count, vertex_size_scalar, Is_Ghost, Is_Anchor) %>%
            dplyr::distinct()
    }

    all_data <- relabel_data

    if (!is.null(ghost_data)) {
        if (collapse_samples) {
            ghost_data <- ghost_data %>%
                dplyr::group_by(Subject_ID, Genotype_Group_ID, SwapCat_ID) %>%
                dplyr::mutate(
                    count = n(),
                    vertex_size_scalar = sqrt(sum(vertex_size_scalar)),
                    Is_Ghost = TRUE,
                    Is_Anchor = FALSE,
                    Sample_ID = ifelse(count == 1, Sample_ID, paste(paste(count, "samples"), Subject_ID, SwapCat_ID, sep="\n"))
                ) %>%
                dplyr::select(Sample_ID, Subject_ID, Genotype_Group_ID, SwapCat_ID, SwapCat_Shape, count, vertex_size_scalar, Is_Ghost, Is_Anchor) %>%
                dplyr::distinct()
        }
        all_data <- rbind(relabel_data, ghost_data)
    }

    if (graph_type == "combined") {
        genotype_graph <- .generate_graph(relabel_data, "genotype", genotype_matrix=genotype_matrix)
        igraph::E(genotype_graph)$genotypes <- TRUE
        label_graph <- .generate_graph(relabel_data, "label", ghost_data)
        igraph::E(label_graph)$labels <- TRUE
        graph <- igraph::graph.union(genotype_graph, label_graph, byname=TRUE)
        igraph::E(graph)[is.na(igraph::E(graph)$genotypes)]$genotypes <- FALSE
        igraph::E(graph)[is.na(igraph::E(graph)$labels)]$labels <- FALSE
        igraph::E(graph)$concordant <- igraph::E(graph)$genotypes & igraph::E(graph)$labels
    } else {
        if (graph_type == "genotype" & !is.null(genotype_matrix)) {
            vertices <- all_data[, "Sample_ID", drop=FALSE]
            graph <- igraph::graph_from_adjacency_matrix(genotype_matrix, mode="undirected")
        } else {
            group_col <- graph_type_mapping[[graph_type]]
            edges <- all_data %>%
                dplyr::group_by_at(group_col) %>%
                dplyr::mutate(
                    sample_a = Sample_ID,
                    sample_b = list(Sample_ID)
                ) %>%
                dplyr::ungroup() %>%
                tidyr::unnest(sample_b) %>%
                dplyr::transmute(
                    sample1 = pmin(sample_a, sample_b),
                    sample2 = pmax(sample_a, sample_b)
                ) %>%
                dplyr::filter(sample1 != sample2) %>%
                dplyr::distinct()
            vertices <- all_data[, "Sample_ID", drop=FALSE]
            graph <- igraph::graph_from_data_frame(edges, vertices=vertices, directed=FALSE)
        }
    }

    if (!populate_plotting_attributes) {return(graph)}

    ## Specify vertex and edge attributes for plotting
    vertex_shapes <- data.frame(Sample_ID=names(igraph::V(graph))) %>%
        dplyr::left_join(all_data, by="Sample_ID") %>%
        dplyr::pull(SwapCat_Shape)
    vertex_size_scalars <- data.frame(Sample_ID=names(igraph::V(graph))) %>%
        dplyr::left_join(all_data, by="Sample_ID") %>%
        dplyr::pull(vertex_size_scalar)
    igraph::V(graph)$shape <- vertex_shapes
    igraph::V(graph)$size <- 12 * vertex_size_scalars
    igraph::V(graph)$label.cex <- 0.5

    anchor_samples <- relabel_data %>%
        dplyr::filter(Is_Anchor) %>%
        dplyr::pull(Sample_ID)
    anchor_samples <- intersect(anchor_samples, igraph::V(graph)$name)
    label_not_found_samples <- relabel_data %>%
        dplyr::filter(stringr::str_detect(Sample_ID, LABEL_NOT_FOUND)) %>%
        dplyr::pull(Sample_ID)
    label_not_found_samples <- intersect(label_not_found_samples, igraph::V(graph)$name)
    ghost_samples <- NULL
    if (!is.null(ghost_data)) {
        ghost_samples <- ghost_data %>% dplyr::pull(Sample_ID)
        ghost_samples <- intersect(ghost_samples, igraph::V(graph)$name)
    }
    igraph::V(graph)$color <- "orange"
    igraph::V(graph)[anchor_samples]$color <- "forestgreen"
    igraph::V(graph)[label_not_found_samples]$color <- "firebrick"
    igraph::V(graph)[ghost_samples]$color <- "lightgrey"

    if (graph_type == "combined") {
        igraph::E(graph)$color <- ifelse(
            igraph::E(graph)$concordant, "forestgreen",
            ifelse(igraph::E(graph)$genotypes, "orange", "cornflowerblue"))
        igraph::E(graph)[.from(ghost_samples)]$color <- "lightgrey"
    } else if (graph_type == "label") {
        igraph::E(graph)$color <- "cornflowerblue"
        igraph::E(graph)[.from(ghost_samples)]$color <- "lightgrey"
    } else {
        igraph::E(graph)$color <- "orange"
    }
    igraph::E(graph)$width <- 6

    return(graph)
}

#' @import dplyr
#'
.relabel_samples <- function(object, relabels) {
    if (nrow(relabels) == 0) return(object)

    relabels <- relabels %>%
        dplyr::rename("Sample_ID" = "relabel_to") %>%
        dplyr::left_join(
            object@sample_metadata[, c("Sample_ID", "Subject_ID")],
            by="Sample_ID"
        ) %>%
        dplyr::mutate(
            Subject_ID = ifelse(
                is.na(Subject_ID),
                vapply(Sample_ID, \(x) unlist(strsplit(x, split="#"))[2], character(1)),
                Subject_ID
            ),
            Deleted_relabel_from = grepl(LABEL_NOT_FOUND, relabel_from)
        ) %>%
        dplyr::filter(!Deleted_relabel_from) %>%
        dplyr::select(-Deleted_relabel_from)

    ## Call it relabeled sample ID instead
    unsolved_all_data <- rbind(object@.solve_state$unsolved_relabel_data,
                               object@.solve_state$unsolved_ghost_data)
    unsolved_all_data <- unsolved_all_data %>%
        dplyr::left_join(
            relabels, by=c("Sample_ID"="relabel_from"), suffix=c(".x", ".y")
        ) %>%
        dplyr::mutate(
            Sample_ID = ifelse(!is.na(Sample_ID.y), Sample_ID.y, Sample_ID),
            Subject_ID = ifelse(!is.na(Subject_ID.y), Subject_ID.y, Subject_ID.x)
        ) %>%
        dplyr::select(-dplyr::ends_with(".x"), -dplyr::ends_with(".y"))

    object@.solve_state$unsolved_relabel_data <- unsolved_all_data %>%
        dplyr::filter(!Is_Ghost)
    object@.solve_state$unsolved_ghost_data <- unsolved_all_data %>%
        dplyr::filter(Is_Ghost)
    object <- .update_solve_state(object)

    print(paste(nrow(relabels), "samples relabeled"))

    return(object)
}

#' @importFrom assertthat assert_that
#' @import igraph
#
# TODO: drop duplicate cycles
.find_directed_cycles <- function(graph, cutoff=1) {
    assertthat::assert_that(igraph::is_directed(graph),
                msg="param 'graph' must be directed")
    cycles <- list()
    for (vertex in igraph::V(graph)) {
        in_neighbors <- names(igraph::neighbors(graph, vertex, mode="in"))
        for (in_neighbor in in_neighbors) {
            simple_paths <- igraph::all_simple_paths(graph, vertex, in_neighbor, mode="out", cutoff=cutoff)
            cycles <- append(cycles, simple_paths)
        }
    }
    cycles <- lapply(cycles, names)
    return(cycles)
}

#' @import igraph
#
.find_all_relabel_cycles <- function(relabels_graph) {
    relabels <- data.frame(relabel_from=character(0), relabel_to=character(0))
    all_relabeled_samples <- NULL
    all_cycles <- list()
    cutoff <- 1
    while (igraph::vcount(relabels_graph) > 0 && cutoff < max(table(igraph::components(relabels_graph)$membership))) {
        curr_cycles <- .find_directed_cycles(relabels_graph, cutoff=cutoff)
        for (curr_cycle in curr_cycles) {
            if (!any(curr_cycle %in% all_relabeled_samples)) {
                all_cycles <- append(all_cycles, list(curr_cycle))
                all_relabeled_samples <- c(all_relabeled_samples, curr_cycle)
                relabels_graph <- igraph::delete_vertices(relabels_graph, curr_cycle)
            }
        }
        cutoff <- cutoff + 1
    }

    ## From the found cycles, construct relabels dataframe
    for (curr_cycle in all_cycles) {
        n <- length(curr_cycle)
        curr_relabels <- data.frame(
            relabel_from = curr_cycle,
            relabel_to = c(curr_cycle[2:n], curr_cycle[1])
        )
        relabels <- rbind(relabels, curr_relabels)
    }

    return(relabels)
}

#' @import dplyr
#' @import igraph
#' @importFrom uuid UUIDgenerate
#'
.find_relabel_cycles_from_putative_subjects <- function(unsolved_relabel_data,
                                                        putative_subjects,
                                                        unsolved_ghost_data=NULL,
                                                        unambiguous_only=FALSE,
                                                        allow_unknowns=FALSE) {
    allow_ghosts <- !is.null(unsolved_ghost_data)

    if (unambiguous_only) {
        allow_ghosts <- FALSE
        allow_unknowns <- FALSE
    }

    mislabel_data <- unsolved_relabel_data %>%
        dplyr::left_join(putative_subjects %>% rename(Putative_Subject_ID = Subject_ID),
                  by="Genotype_Group_ID") %>%
        dplyr::mutate(
            Inferred_Correctly_Labeled = Putative_Subject_ID == Subject_ID,
            ## Only include samples where the current label is a Subject_ID
            ## that also has a Genotype_Group_ID assigned
            Curr_Subject_ID_Genotyped = Subject_ID %in% putative_subjects$Subject_ID
        ) %>%
        dplyr::filter(!is.na(Putative_Subject_ID) & !Inferred_Correctly_Labeled & Curr_Subject_ID_Genotyped)

    relabels <- EMPTY_RELABELS

    if (nrow(mislabel_data) == 0) {
        return(relabels)
    }

    ## Every mislabeled sample with the same genotype and same swappable category
    ## must also have the same potential relabels. We search for relabels at
    ## at the SwapCat_ID/Genotype_Group_ID level
    mislabeled_genotype_swapcats <- mislabel_data %>%
        dplyr::select(SwapCat_ID, Genotype_Group_ID) %>%
        dplyr::distinct() %>%
        dplyr::left_join(putative_subjects, by="Genotype_Group_ID")

    directed_edge_mats <- vector("list", length=nrow(mislabeled_genotype_swapcats))

    all_ghost_labels <- character(0)
    all_unknown_labels <- character(0)

    ## Iterate over all SwapCat_ID/Genotype_Group_ID pairs and search for
    ## potential relabels
    for (i in seq_len(nrow(mislabeled_genotype_swapcats))) {
        swap_cat_id <- mislabeled_genotype_swapcats[i, "SwapCat_ID"]
        genotype_group_id <- mislabeled_genotype_swapcats[i, "Genotype_Group_ID"]
        subject_id <- mislabeled_genotype_swapcats[i, "Subject_ID"]

        ## These are all the mislabeled samples for a SwapCat_ID/Genotype_Group_ID pair
        mislabeled_samples <- mislabel_data %>%
            dplyr::filter(SwapCat_ID == swap_cat_id, Genotype_Group_ID == genotype_group_id) %>%
            dplyr::pull(Sample_ID)
        n_mislabeled_samples <- length(mislabeled_samples)

        ## The eligible relabels have the putative Subject_ID but are in a different Genotype_Group
        eligible_labels <- unsolved_relabel_data %>%
            dplyr::filter(SwapCat_ID == swap_cat_id, Subject_ID == subject_id, Genotype_Group_ID != genotype_group_id) %>%
            dplyr::pull(Sample_ID)
        n_eligible_labels <- length(eligible_labels)

        n_label_deficit <- n_mislabeled_samples - n_eligible_labels

        ghost_labels <- character(0)
        n_ghost_labels <- 0
        ## If there aren't enough eligible labels, try using ghosts to plug the gap
        if (allow_ghosts && n_label_deficit > 0) {
            ghost_labels <- unsolved_ghost_data %>%
                dplyr::filter(SwapCat_ID == swap_cat_id, Subject_ID == subject_id) %>%
                dplyr::pull(Sample_ID)
            n_ghost_labels <- length(ghost_labels)

            ## If we have more ghost labels than needed, select a subset
            if (n_label_deficit < n_ghost_labels) {
                ghost_labels <- ghost_labels[seq_len(n_label_deficit)]
                n_ghost_labels <- length(ghost_labels)
            }

            n_label_deficit <- n_mislabeled_samples - n_eligible_labels - n_ghost_labels
        }

        unknown_labels <- character(0)
        n_unknown_labels <- 0
        ## If there still aren't enough eligible labels, resort to plugging the gap with unknowns
        if (allow_unknowns && n_label_deficit > 0) {
            unknown_labels <- paste0(LABEL_NOT_FOUND, "#", subject_id, "#", swap_cat_id, "#", UUIDgenerate(n=n_label_deficit))
            n_unknown_labels <- length(unknown_labels)
        }

        ## We should have n_mislabeled_samples == n_eligible_labels + n_ghost_labels + n_unknown_labels

        directed_edge_mats[[i]] <- expand.grid(relabel_from=mislabeled_samples, relabel_to=c(eligible_labels, ghost_labels, unknown_labels))
        all_ghost_labels <- c(all_ghost_labels, ghost_labels)
        all_unknown_labels <- c(all_unknown_labels, unknown_labels)
    }

    ## This graph doesn't have outgoing edges from ghosts or unknowns
    ## So there are no cycles that include ghosts or unknowns yet
    directed_edge_df <- as.data.frame(do.call(rbind, directed_edge_mats))
    relabels_graph <- igraph::graph_from_data_frame(directed_edge_df, directed=TRUE)

    ## If unambiguous_only, only include samples with exactly one incoming and one outgoing edge
    if (unambiguous_only) {
        samples_with_one_incoming <- igraph::V(relabels_graph)[igraph::degree(relabels_graph, mode = "in") == 1]$name
        samples_with_one_outgoing <- igraph::V(relabels_graph)[igraph::degree(relabels_graph, mode = "out") == 1]$name
        samples_to_filter <- intersect(samples_with_one_incoming, samples_with_one_outgoing)
        relabels_graph <- igraph::subgraph(relabels_graph, samples_to_filter)
    }

    ## 1. Find relabel cycles without using ghosts or unknowns
    new_relabels <- .find_all_relabel_cycles(relabels_graph)
    relabels_graph <- relabels_graph - new_relabels$relabel_from
    relabels <- rbind(relabels, new_relabels)

    ## 2. Find relabel cycles allowing for ghosts
    if (allow_ghosts) {
        all_excess_labels <- names(igraph::V(relabels_graph)[igraph::degree(relabels_graph, mode = "in") == 0])
        igraph::V(relabels_graph)$relabel_component_id <- igraph::components(relabels_graph)$membership
        all_relabel_component_ids <- unique(igraph::V(relabels_graph)$relabel_component_id)
        ## For each component, connect ghost labels with labels that have no more incoming edges,
        ## indicating that there are no genotyped samples that can take their label
        for (curr_relabel_component_id in all_relabel_component_ids) {
            component_labels <- names(igraph::V(relabels_graph)[igraph::V(relabels_graph)$relabel_component_id == curr_relabel_component_id])
            component_ghost_labels <- intersect(all_ghost_labels, component_labels)
            component_excess_labels <- intersect(all_excess_labels, component_labels)
            new_edges <- expand.grid(relabel_from=component_ghost_labels, relabel_to=component_excess_labels)
            if (nrow(new_edges) > 0) {
                new_relabels_graph <- igraph::graph_from_data_frame(new_edges)
                relabels_graph <- igraph::union(relabels_graph, new_relabels_graph)
            }
        }
        new_relabels <- .find_all_relabel_cycles(relabels_graph)
        relabels_graph <- relabels_graph - new_relabels$relabel_from
        relabels <- rbind(relabels, new_relabels)
    }

    ## 3. Find relabel cycles allowing for unknowns
    if (allow_unknowns) {
        all_excess_labels <- names(igraph::V(relabels_graph)[igraph::degree(relabels_graph, mode = "in") == 0])
        igraph::V(relabels_graph)$relabel_component_id <- igraph::components(relabels_graph)$membership
        all_relabel_component_ids <- unique(igraph::V(relabels_graph)$relabel_component_id)
        ## For each component, connect unknown labels with labels that have no more incoming edges,
        ## indicating that there are no genotyped samples that can take their label
        for (curr_relabel_component_id in all_relabel_component_ids) {
            component_labels <- names(igraph::V(relabels_graph)[igraph::V(relabels_graph)$relabel_component_id == curr_relabel_component_id])
            component_unknown_labels <- intersect(all_unknown_labels, component_labels)
            component_excess_labels <- intersect(all_excess_labels, component_labels)
            new_edges <- expand.grid(relabel_from=component_unknown_labels, relabel_to=component_excess_labels)
            if (nrow(new_edges) > 0) {
                new_relabels_graph <- igraph::graph_from_data_frame(new_edges)
                relabels_graph <- igraph::union(relabels_graph, new_relabels_graph)
            }
        }
        new_relabels <- .find_all_relabel_cycles(relabels_graph)
        relabels_graph <- relabels_graph - new_relabels$relabel_from
        relabels <- rbind(relabels, new_relabels)
    }

    return(relabels)
}

#' @import dplyr
#' @import igraph
#' @import Matrix
#'
.find_neighbors <- function(object, include_ghost=FALSE, filter_concordant_vertices=FALSE) {
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
    unsolved_ghost_data <- NULL
    if (include_ghost) {
        unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data
    }
    combined_graph <- .generate_graph(unsolved_relabel_data, graph_type="combined", unsolved_ghost_data)
    unsolved_all_data <- rbind(unsolved_relabel_data, unsolved_ghost_data)
    putative_subjects <- object@.solve_state$putative_subjects

    ## Criteria 1:  filter only pairs of vertices that are within at exactly 2 edges of each other
    adj_matrix_sparse <- Matrix::Matrix(igraph::get.adjacency(combined_graph, sparse = TRUE))
    adj_matrix_idx <- Matrix::which(adj_matrix_sparse > 0, arr.ind = TRUE)
    dist_within_2_sparse <- adj_matrix_sparse %*% adj_matrix_sparse
    dist_within_2_idx <- Matrix::which(dist_within_2_sparse > 0, arr.ind = TRUE)
    unique_pairs_exact1 <- data.frame(
        Row = rownames(adj_matrix_sparse)[adj_matrix_idx[, 1]],
        Col = colnames(adj_matrix_sparse)[adj_matrix_idx[, 2]]) %>%
        dplyr::transmute(
            Sample_A = pmin(Row, Col),
            Sample_B = pmax(Row, Col)
        ) %>%
        dplyr::filter(Sample_A != Sample_B) %>%
        dplyr::distinct()
    unique_pairs_within2 <- data.frame(
        Row = rownames(dist_within_2_sparse)[dist_within_2_idx[, 1]],
        Col = colnames(dist_within_2_sparse)[dist_within_2_idx[, 2]]) %>%
        dplyr::transmute(
            Sample_A = pmin(Row, Col),
            Sample_B = pmax(Row, Col)
        ) %>%
        dplyr::filter(Sample_A != Sample_B) %>%
        dplyr::distinct()
    unique_pairs <- dplyr::anti_join(unique_pairs_within2, unique_pairs_exact1, by=c("Sample_A", "Sample_B"))

    ## Criteria 2: filter out pairs of vertices that include elements in anchor_samples
    unique_pairs <- unique_pairs[!(unique_pairs$Sample_A %in% object@anchor_samples |unique_pairs$Sample_B %in% object@anchor_samples), ]

    ## Criteria 3: filter only pairs of vertices that are within the same swap category
    unique_pairs <- unique_pairs %>%
        dplyr::left_join(unsolved_all_data[, c("Sample_ID", "SwapCat_ID")], by=c("Sample_A"="Sample_ID")) %>%
        dplyr::rename("SwapCat_A"="SwapCat_ID") %>%
        dplyr::left_join(unsolved_all_data[, c("Sample_ID", "SwapCat_ID")], by=c("Sample_B"="Sample_ID")) %>%
        dplyr::rename("SwapCat_B"="SwapCat_ID") %>%
        dplyr::filter(SwapCat_A == SwapCat_B) %>%
        dplyr::select(Sample_A, Sample_B)

    ## Criteria 4: filter out vertices that have at least 1 concordant edge
    if (filter_concordant_vertices) {
        concordant_edges <- igraph::E(combined_graph)[igraph::E(combined_graph)$concordant]
        concordant_vertices <- unique(c(igraph::ends(combined_graph, concordant_edges)))
        unique_pairs <- unique_pairs %>%
            dplyr::filter(
                !(Sample_A %in% concordant_vertices),
                !(Sample_B %in% concordant_vertices)
            )
    }

    ## Criteria 5: filter out pairs of vertices where either side will violate putative_subjects
    unique_pairs <- unique_pairs %>%
        dplyr::left_join(
            unsolved_all_data[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID")],
            by=c("Sample_A"="Sample_ID")) %>%
        dplyr::rename("Subject_A"="Subject_ID", "Genotype_Group_A"="Genotype_Group_ID") %>%
        dplyr::left_join(
            putative_subjects %>% dplyr::filter(!is.na(Genotype_Group_ID)),
            by=c("Genotype_Group_A"="Genotype_Group_ID")) %>%
        dplyr::rename("Putative_Subject_A"="Subject_ID") %>%
        dplyr::left_join(
            unsolved_all_data[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID")],
            by=c("Sample_B"="Sample_ID")) %>%
        dplyr::rename("Subject_B"="Subject_ID", "Genotype_Group_B"="Genotype_Group_ID") %>%
        dplyr::left_join(
            putative_subjects %>% dplyr::filter(!is.na(Genotype_Group_ID)),
            by=c("Genotype_Group_B"="Genotype_Group_ID")) %>%
        dplyr::rename("Putative_Subject_B"="Subject_ID") %>%
        dplyr::mutate(
            # The swap is invalid if
            # 1. Sample_A is already in its Putative_Subject
            # 2. Sample_B is already in its Putative_Subject
            # For swaps where the putative subjects of both samples are assigned, the swap is also invalid if
            # 3. Neither Sample_A nor Sample_B match their Putative_Subject after
            Invalid_Swap =
                (!is.na(Putative_Subject_A) & Subject_A == Putative_Subject_A) |
                (!is.na(Putative_Subject_B) & Subject_B == Putative_Subject_B) |
                (!is.na(Putative_Subject_A) & !is.na(Putative_Subject_B) & Putative_Subject_A != Subject_B & Putative_Subject_B != Subject_A)
        ) %>%
        dplyr::filter(!Invalid_Swap) %>%
        dplyr::select(Sample_A, Sample_B)

    return(unique_pairs)
}
