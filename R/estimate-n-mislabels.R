#' Normalize a vector of sexes to a factor with 2 levels
#'
#' @param x Vector to be normalized
#' @param male_synonyms Vector of strings that should be considered equivalent
#'   to "Male". Case-insensitive.
#' @param female_synonyms Vector of strings that should be considered equivalent
#'   to "Female". Case-insensitive.
#'
#' @returns A factor the same length as `x` with levels "Female" and Male".
#' @export
#'
#' @examples
#' x <- c("m", "M", "male", "MALE", "Male", "f", "F", "Female", "FEMALE", "FeMale")
#' s <- normalize_sex(x)
#' print(s)
#' table(s)
normalize_sex <- function(
        x,
        male_synonyms = c("male", "m"),
        female_synonyms = c("female", "f")
) {
    x <- str_to_lower(x)
    male_synonyms <- str_to_lower(male_synonyms)
    female_synonyms <- str_to_lower(female_synonyms)
    assert_that(!any(male_synonyms %in% female_synonyms))
    result <- case_when(
        x %in% male_synonyms ~ "Male",
        x %in% female_synonyms ~ "Female",
        is.na(x) ~ NA_character_,
        .default = "INVALID"
    )
    assert_that(!any(na.omit(result) == "INVALID"), msg = str_c("Invalid sex value(s): ", deparse1(head(unique(x[result == "INVALID"])))))
    factor(result, levels = c("Female", "Male"))
}

#' Estimate the total number of mislabels based on known sex mismatches
#'
#' @param reported_sex Vector of reported sex for each sample.
#' @param inferred_sex Vector of inferred sex for each sample. Must be the same
#'   length and same order as `reported_sex`.
#' @param return_fraction If TRUE, return the estimated fraction of mislabeled
#'   samples instead of the estimated number.
#'
#' @returns The estimated number (or fraction) of mislabeled samples.
#'
#' @details If `reported_sex` or `inferred_sex` contains NAs, those samples will
#' not be used to estimate the fraction of mislabeled samples
#'
#' @export
#'
#' @examples
#'
#' # Simulate 100 samples with 20 mislabels
#' set.seed(1986)
#' reported <- sample(c("m", "f"), 100, replace = TRUE)
#' # Simulate 50 more samples with unknown sex
#' reported <- c(reported, rep(NA, 50))
#' inferred <- reported
#' # Simulate mislabels
#' inferred[1:20] <- sample(inferred[1:20])
#' # Should give approximately 30
#' estimate_n_mislabels_from_sex_mismatches(reported, inferred)
estimate_n_mislabels <- function(reported_sex, inferred_sex, return_fraction = FALSE) {
    # Also allow sex to be passed as logical. We arbitrarily assign each boolean
    # to a sex, since it doesn't matter which is which for the purposes of this
    # calculation.
    if (is_logical(reported_sex) && is_logical(inferred_sex)) {
        reported_sex <- if_else(reported_sex, "M", "F")
        inferred_sex <- if_else(inferred_sex, "M", "F")
    }
    reported_sex <- normalize_sex(reported_sex)
    inferred_sex <- normalize_sex(inferred_sex)
    assert_that(
        length(reported_sex) > 1,
        msg = "Not enough samples to estimate number of mislabels"
    )
    assert_that(
        length(reported_sex) == length(inferred_sex)
    )
    x <- tibble(reported_sex, inferred_sex) |>
        mutate(
            evaluable = !is.na(reported_sex) & !is.na(inferred_sex),
            matched = reported_sex == inferred_sex
        )
    if (!adjust_for_unevaluable) {
        x <- drop_na(x)
    }
    n_samples <- nrow(x)
    n_evaluable <- sum(x$evaluable)
    n_mismatch <- sum(x$evaluable & !x$matched)
    n_male <- sum(x$evaluable & (x$reported_sex == "Male"))
    frac_male <- n_male / n_evaluable
    frac_mismatch <- n_mismatch / n_evaluable

    est_frac_mislabel <- frac_mismatch / (2 * frac_male * (1 - frac_male))
    if (return_fraction) {
        return(est_frac_mislabel)
    }
    # We multiply by the total number of samples, including non-evaluable ones
    # excluded from the calculation of est_frac_mislabel
    est_n_mislabel <- n_samples * est_frac_mislabel
    return(est_n_mislabel)
}
