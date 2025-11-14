#' Count events with optional incidence rates and exact Poisson CI
#'
#' `count_bio()` is a small helper that works a bit like [dplyr::count()]:
#' it does a grouped count, and optionally also sums person-time and returns
#' incidence rates with Garwood (exact Poisson) confidence intervals.
#'
#' @param x A data frame.
#' @param ... Variables to group by.
#' @param events Optional numeric (typically 0/1 or non-negative counts).
#'   If `NULL`, the function uses `dplyr::n()` to count rows.
#' @param person_time Optional non-negative numeric column giving person-time.
#'   If supplied, the output will include person-time and incidence rates.
#' @param per Positive scalar used to scale the rate
#'   (for example, `per = 100` gives rates per 100 person-time units).
#' @param sort Logical. If `TRUE`, sort the result by the event count
#'   in descending order.
#' @param name Name of the event count column. If `NULL`, a name like
#'   `"n"`, `"nn"`, ... is chosen to avoid overwriting existing columns.
#' @param .drop Passed on to [dplyr::group_by()].
#' @param conf_level Confidence level for the Poisson confidence interval.
#'   Must be a number between 0 and 1 (default is 0.95).
#'
#' @return A data frame with the grouping columns, the event count column
#'   (`name`), and, if `person_time` is supplied, extra columns:
#'   `pt` (total person-time), `rate`, `rate_lcl`, `rate_ucl`.
#'
#' @examples
#' library(dplyr)
#'
#' # toy data
#' df <- tibble(
#'   trt = rep(c("A", "B"), each = 5),
#'   events = rbinom(10, 1, 0.3),
#'   pt = rexp(10, rate = 0.2)
#' )
#'
#' # just counts, similar to dplyr::count()
#' count_bio(df, trt)
#'
#' # counts + person-time + incidence rates per 100 person-time units
#' count_bio(df, trt,
#'           events = events,
#'           person_time = pt,
#'           per = 100)
#'
#' @export
#' @importFrom dplyr group_by summarise arrange desc n mutate
#' @importFrom rlang enquo enquos quo_is_null expr inform
#' @importFrom stats qchisq
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom data.table :=



count_bio <- function(x,
                      ...,
                      events = NULL,
                      person_time = NULL,
                      per = 1,
                      sort = FALSE,
                      name = NULL,
                      .drop = dplyr::group_by_drop_default(x),
                      conf_level = 0.95) {

  # basic checks
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame.", call. = FALSE)
  }
  if (!is.numeric(per) || length(per) != 1L || per <= 0) {
    stop("`per` must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be in (0, 1).", call. = FALSE)
  }

  # capture tidy-eval arguments
  events_quo      <- rlang::enquo(events)
  person_time_quo <- rlang::enquo(person_time)
  group_quos      <- rlang::enquos(...)

  # pick a name for the count column
  name <- check_n_name_bio(name, names(x))

  # how to compute the event count
  count_expr <- if (rlang::quo_is_null(events_quo)) {
    rlang::expr(dplyr::n())
  } else {
    rlang::expr(sum(!!events_quo, na.rm = TRUE))
  }

  # also have person-time
  has_pt <- !rlang::quo_is_null(person_time_quo)

  if (has_pt) {
    # grouped count + total person-time
    out <- x %>%
      dplyr::group_by(!!!group_quos, .drop = .drop) %>%
      dplyr::summarise(
        !!name := !!count_expr,
        pt      = sum(!!person_time_quo, na.rm = TRUE),
        .groups = "drop"
      )

    # negative person-time means bad data
    if (any(out$pt < 0, na.rm = TRUE)) {
      stop("Computed `pt` contains negative values; check `person_time`.", call. = FALSE)
    }

    # exact Poisson CI for the rate
    alpha <- 1 - conf_level
    k     <- out[[name]]
    total_pt <- out$pt

    # base rate per unit person-time
    rate <- ifelse(total_pt > 0, k / total_pt, NA_real_)

    lower_cnt <- 0.5 * stats::qchisq(alpha / 2,     2 * k)
    upper_cnt <- 0.5 * stats::qchisq(1 - alpha / 2, 2 * (k + 1))

    rate_lcl <- ifelse(total_pt > 0, lower_cnt / total_pt, NA_real_)
    rate_ucl <- ifelse(total_pt > 0, upper_cnt / total_pt, NA_real_)

    out <- out %>%
      dplyr::mutate(
        rate     = rate     * per,
        rate_lcl = rate_lcl * per,
        rate_ucl = rate_ucl * per
      )
  } else {
    out <- x %>%
      dplyr::group_by(!!!group_quos, .drop = .drop) %>%
      dplyr::summarise(
        !!name := !!count_expr,
        .groups = "drop"
      )
  }

  if (isTRUE(sort)) {
    out <- dplyr::arrange(out, dplyr::desc(.data[[name]]))
  }

  out
}

# helpers ---------------------------------------------------------------

# choose a name for the count column that does not clash with existing names
check_n_name_bio <- function(name, vars) {
  if (is.null(name)) {
    name <- n_name_bio(vars)
    if (name != "n") {
      rlang::inform(c(
        paste0("Storing counts in `", name, "`, as `n` already present in input"),
        i = "Use `name = \"new_name\"` to pick a new name."
      ))
    }
  } else {
    if (!is.character(name) || length(name) != 1L) {
      stop("`name` must be a single string.", call. = FALSE)
    }
  }
  name
}

n_name_bio <- function(vars) {
  nm <- "n"
  while (nm %in% vars) {
    nm <- paste0("n", nm)
  }
  nm
}
