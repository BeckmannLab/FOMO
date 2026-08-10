if (!isGeneric("plot")) {
    setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
}

#' setGeneric("solveMajoritySearch", function(object, ...)
#'     standardGeneric("solveMajoritySearch"))
#'
#' setGeneric("solve_global_search", function(object, ...)
#'     standardGeneric("solve_global_search"))
#'
#' setGeneric("solve_local_search", function(object, ...)
#'     standardGeneric("solve_local_search"))
#'
#' setGeneric("solve_local_search_bulk", function(object, ...)
#'     standardGeneric("solve_local_search_bulk"))
#'
#' setGeneric("solve", function(object, ...)
#'     standardGeneric("solve"))
#'
#' setGeneric("write_corrections", function(object, ...)
#'     standardGeneric("write_corrections"))
