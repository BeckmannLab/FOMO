if (!methods::isGeneric("plot")) {
    methods::setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
}

#' methods::setGeneric("solveMajoritySearch", function(object, ...)
#'     standardGeneric("solveMajoritySearch"))
#'
#' methods::setGeneric("solve_global_search", function(object, ...)
#'     standardGeneric("solve_global_search"))
#'
#' methods::setGeneric("solve_local_search", function(object, ...)
#'     standardGeneric("solve_local_search"))
#'
#' methods::setGeneric("solve_local_search_bulk", function(object, ...)
#'     standardGeneric("solve_local_search_bulk"))
#'
#' methods::setGeneric("solve", function(object, ...)
#'     standardGeneric("solve"))
#'
#' methods::setGeneric("write_corrections", function(object, ...)
#'     standardGeneric("write_corrections"))
