# This is a workaround for the note "No visible binding for global variable"
# https://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
#
# ".from" is used as a symbol in igraph's edge-selector syntax, used in
# .generate_graph() (e.g. E(g)[.from(1)]), so it needs to be declared as a
# global variable so that `R CMD check` doesn't complain about it.
global_vars <- c(
    ".from"
)
globalVariables(global_vars)
