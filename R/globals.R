# This is a workaround for the note "No visible binding for global variable"
# https://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
#
# This file used to carry a much longer list here: every bare column name
# used inside a dplyr/NSE context anywhere in the package (mutate(), filter(),
# transmute(), summarize(), arrange(), group_by(), pull(), etc.), one entry
# per name, regenerated periodically by running
# codetools::checkUsagePackage("fomo", all = TRUE) (after
# pkgload::load_all()) and copying out every name reported under either of
# the two NOTE categories globalVariables() actually suppresses ("no visible
# binding for global variable" and "no visible global function definition
# for"). That approach is fragile by nature -- nothing enforces that the
# list stays in sync with the code, so it silently drifts every time a call
# site is added, renamed, or removed, and it's easy for someone to not
# notice until R CMD check is run (which is exactly what had happened to the
# previous version of this file: it had drifted enough that ~20 more
# bindings had accumulated in the package with no corresponding entries
# here, while globalVariables() itself had been commented out at some point
# and stopped running entirely, so neither the old list nor a fresh one was
# actually being applied).
#
# The more thorough fix described here previously -- following the dplyr
# non-standard-evaluation vignette's .data$-prefix pattern throughout the
# package instead --
# https://cran.r-project.org/web/packages/dplyr/vignettes/in-packages.html
# -- has now been done, across every R file in the package. Every bare
# column reference inside a data-masking context (mutate(), transmute(),
# filter(), arrange(), summarize()/summarise(), group_by(), pull(), case_when()
# nested within these, etc.) is now written as .data$column_name; every bare
# column reference inside a tidyselect-only context (select(), rename() --
# specifically the existing/old-name side -- relocate(), ungroup(), unnest())
# is written as a quoted string instead, since tidyselect deprecated .data$
# usage as of tidyselect 1.2.0 (using .data$ there produces its own
# deprecation warning instead of fixing anything). This makes each call site
# self-documenting and correct by construction, rather than depending on a
# separate list staying accurate over time.
#
# The one exception is ".from" below, which was never actually a
# globalVariables()-style NSE issue in the first place -- it's igraph's
# edge-selector syntax, used in .generate_graph() (e.g. E(g)[.from(1)]).
# `.from` doesn't exist anywhere in igraph's namespace, exported or
# internal; `[.igraph.es` recognizes the literal, unevaluated call as
# special selector syntax rather than ever actually calling a function named
# .from(). That's a well-known false-positive pattern when statically
# analyzing igraph-based code (codetools reports it as "no visible global
# function definition for '.from'"), and there is no .data$/quoted-string
# equivalent for it -- globalVariables() remains the standard, documented
# way to suppress this specific false positive.
global_vars <- c(
    ".from"
)
globalVariables(global_vars)
