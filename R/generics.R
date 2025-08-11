#' @import methods
if (!methods::isGeneric("plot"))
    methods::setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

#' #' @import methods
#' #'
#' methods::setGeneric("solveMajoritySearch", function(object, ...)
#'     standardGeneric("solveMajoritySearch"))
#'
#' #' @import methods
#' #'
#' methods::setGeneric("solve_comprehensive_search", function(object, ...)
#'     standardGeneric("solve_comprehensive_search"))
#'
#' #' @import methods
#' #'
#' methods::setGeneric("solve_local_search", function(object, ...)
#'     standardGeneric("solve_local_search"))
#'
#' #' @import methods
#' #'
#' methods::setGeneric("solve_local_search_bulk", function(object, ...)
#'     standardGeneric("solve_local_search_bulk"))
#'
#' #' @import methods
#' #'
#' methods::setGeneric("solve", function(object, ...)
#'     standardGeneric("solve"))
#'
#' #' @import methods
#' #'
#' methods::setGeneric("write_corrections", function(object, ...)
#'     standardGeneric("write_corrections"))
