# Functions related with descriptive statistics.

#' Relative Dispersion.
#'
#' Calculates the coefficient of variation (relative dispersion) of a variable. The relative dispersion
#' is defined as the standard deviation over the arithmetic mean.
#'
#' @param x A numerical variable. NA's observations are removed by default.
#' @return The coefficient of variation (relative dispersion).
#' @examples
#' height <- rnorm(100, 170, 8)
#' rel_dis(height)
#' @export
rel_dis <- function(x) {
  rd <- sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  rd
}

#' Reference range (reference interval).
#'
#' \code{reference_range} estimates the reference range (reference interval) of a numerical variable.
#'
#' The reference range assumes normality and represents the limits that would include 95% of the expected
#' observations.
#'
#' @param avg The arithmetic mean (a scalar numerical value).
#' @param std The standard deviation (a scalar numerical value).
#' @return A data frame with the reference range limits.
#' @examples
#' x <- rnorm(100, 170, 8)
#' round(mean(x), 2)
#' round(sd(x), 2)
#'
#' round(reference_range(mean(x), sd(x)), 2)
#' @export
reference_range <- function(avg, std) {
  lower.ri <- avg - qnorm(0.975) * std
  upper.ri <- avg + qnorm(0.975) * std
  ri <- data.frame(lower.ri, upper.ri)
  ri
}

#' Geometric mean.
#'
#' @param x A numeric variable with no negative values.
#' @return A scalar, the calculated geometric mean.
#' @examples
#' data(IgM, package = "ISwR")
#' Ab <- data.frame(IgM)
#' estat(~IgM, data = Ab)
#' geo_mean(IgM)
#' @export
geo_mean <- function(x) {
  positive.check <- which(x <= 0)
  if (length(positive.check) >= 1) {
    print("Some observations are not positive")
  } else {
    lx <- length(x)
    geo <- exp(sum(log(x), na.rm = TRUE) / lx)
    geo
  }
}

#' Harmonic mean.
#'
#' @param x A numeric variable with no zero values.
#' @return A scalar, the calculated harmonic mean.
#' @examples
#' data(IgM, package = "ISwR")
#' Ab <- data.frame(IgM)
#' estat(~IgM, data = Ab)
#' harm_mean(IgM)
#' @export
harm_mean <- function(x) {
  zero.check <- which(x == 0)
  if (length(zero.check) >= 1) {
    print("Some observations are equal to zero")
  } else {
    lx <- length(x)
    inv <- sum(1 / x, na.rm = TRUE)
    harm <- lx / inv
    harm
  }
}

#' Bootstrap Confidence Intervals.
#'
#' \code{bst} estimates confidence intervals around the \code{\link[base]{mean}}, \code{\link[stats]{median}} or \code{\link{geo_mean}}.
#'
#' @param x A numerical variable. Missing observations are removed by default.
#' @param stat Statistic, either "mean" (default), "median" or "gmean" (geometric mean).
#' @param n Number of replicates for the bootstrap (n=1000 by default).
#' @param CI Confidence intervals (CI=95 by default).
#' @param digits Number of digits for rounding (default = 2).
#' @return A data frame with the estimate and confidence intervals.
#' @examples
#' data(IgM, package = "ISwR")
#' bst(IgM, "median")
#'
#' bst(IgM, "gmean")
#' @export
bst <- function(x, stat = "mean", n = 1000, CI = 95, digits = 2) {
  xmeans <- numeric(n)
  if (stat == "median") {
    for (i in 1:n) {
      xmeans[i] <- median(sample(x, replace = TRUE),
        na.rm = TRUE
      )
    }
    estimate <- median(x, na.rm = TRUE)
  } else if (stat == "gmean") {
    for (i in 1:n) xmeans[i] <- geo_mean(sample(x, replace = TRUE))
    estimate <- geo_mean(x)
  } else if (stat == "mean") {
    for (i in 1:n) {
      xmeans[i] <- mean(sample(x, replace = TRUE),
        na.rm = TRUE
      )
    }
    estimate <- mean(x, na.rm = TRUE)
  } else {
    stop("Statistic not available")
  }
  alpha <- 1 - (CI / 100)
  tail <- alpha / 2
  lower <- round(quantile(xmeans, tail, na.rm = TRUE), digits = digits)
  upper <- round(quantile(xmeans, 1 - tail, na.rm = TRUE), digits = digits)
  estimate <- round(estimate, digits = digits)
  cis <- data.frame(stat, estimate, CI, lower, upper)
  rownames(cis) <- ""
  colnames(cis)[3] <- "%CI"
  cis
}

#' Internal function to calculate descriptive statistics.
#'
#' \code{stats_quotes} is an internal function called by \code{estat}.
#'
#' @param x a numeric variable
#' @param data2 A data frame where \code{x} can be found.
#' @param digits Number of digits for rounding.
stats_quotes <- function(x, data2, digits = 2) {
  n <- sum(!is.na(data2[[x]]))
  x.min <- min(data2[[x]], na.rm = TRUE)
  x.max <- max(data2[[x]], na.rm = TRUE)
  avg <- mean(data2[[x]], na.rm = TRUE)
  med <- median(data2[[x]], na.rm = TRUE)
  std <- sd(data2[[x]], na.rm = TRUE)
  cv <- std / avg
  res <- data.frame(n, x.min, x.max, avg, med, std, cv)
  names(res) <- c("N", "Min", "Max", "Mean", "Median", "SD", "CV")
  res <- round(res, digits = digits)
  res
}

#' Descriptive statistics for continuous variables.
#'
#' \code{estat} calculates descriptives of numerical variables.
#'
#' @param object When chaining, this holds an object produced in the earlier portions of the chain. Most users can safely ignore this argument. See details and examples.
#' @param formula A formula with shape: ~ x or ~ x|z (for groups).
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param digits Number of digits for rounding (default = 2).
#' @param label Label used to display the name of the variable (see examples).
#' @return A data frame with descriptive statistics.
#' @seealso \code{\link{summary}}.
#' @examples
#' data(kfm, package = "ISwR")
#' require(sjlabelled, quietly = TRUE)
#' kfm <- kfm |>
#'   var_labels(
#'     dl.milk = "Breast-milk intake (dl/day)",
#'     sex = "Sex",
#'     weight = "Child weight (kg)",
#'     ml.suppl = "Milk substitute (ml/day)",
#'     mat.weight = "Maternal weight (kg)",
#'     mat.height = "Maternal height (cm)"
#'   )
#'
#' kfm |>
#'   estat(~dl.milk)
#'
#' estat(~ dl.milk | sex, data = kfm)
#'
#' kfm |>
#'   estat(~ weight | sex)
#' @export
estat <- function(object = NULL, formula = NULL, data = NULL, digits = 2, label = NULL) {
  if (inherits(object, "formula")) {
    formula <- object
    object <- NULL
  }
  if (inherits(object, "data.frame")) {
    data <- object
    object <- NULL
  }
  vars <- all.vars(formula)
  y <- vars[1]
  outcome <- data[[y]]
  if (!is.numeric(outcome)) {
    stop(y, " must be a numerical variable")
  }
  if (is.null(sjlabelled::get_label(outcome)) == FALSE & is.null(label)) {
    lab <- sjlabelled::get_label(outcome)
  } else {
    lab <- ifelse(is.null(label), y, label)
  }
  nv <- length(vars)
  if (nv == 1) {
    res <- stats_quotes(y, data2 = data)
    res <- round(res, digits = digits)
    res <- data.frame(lab, res)
    names(res)[1] <- ""
    res
  } else {
    x <- vars[2]
    exposure <- data[[x]]
    if (!is.factor(exposure)) {
      stop(x, " must be a factor")
    }
    lev <- levels(exposure)
    nl <- length(lev)
    res <- numeric(nl * 7)
    dim(res) <- c(nl, 7)
    for (i in 1:nl) {
      res[i, ] <- as.numeric(stats_quotes(y,
        data2 = subset(data, exposure == lev[i])
      ))
    }
    res <- as.data.frame(res)
    names(res) <- c(
      "N", "Min", "Max", "Mean", "Median",
      "SD", "CV"
    )
    res <- round(res, digits = digits)
    res <- data.frame(var = lev, res)
    names(res)[1] <- x
    var <- character(nl)
    var[1] <- lab
    res <- data.frame(var, res)
    names(res)[1] <- ""
    names(res)[2] <- ifelse(is.null(sjlabelled::get_label(exposure)) == FALSE,
      sjlabelled::get_label(exposure), x
    )
    res
  }
}

#' Table of descriptive statistics by categorical variable.
#'
#' \code{cross_tbl} is a wrapper to function  from package \code{\link[gtsummary]{tbl_summary}} that constructs tables of descriptive statistics stratified by levels of a categorical outcome.
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param by The quoted name of the  categorical variable (factor) used for the stratification.
#' @param head_label Character, label to be used as head for the variable's column.
#' @param bold Display labels in bold?
#' @param show_total Logical, show column with totals?
#' @param p_val Logical, show p-values?
#' @param pad Numerical, padding above and bellow rows.
#' @param method An integer indicating methods for continuous variables.
#' 1 Reports means and standard deviations.
#' 2 Reports medians and interquartile ranges.
#' @param ... Additional arguments passed to \code{\link[gtsummary]{tbl_summary}}.
#' @details Function \code{cross_tbl} is a relatively simple wrapper to function \code{\link[gtsummary]{tbl_summary}}. It constructs contingency tables and can also be used to report a table with descriptives for all variables stratified by one of the variables. Please see examples to see how to list variables. If data is labelled, the label of the stratifying variable is used as part of the header.
#' @return A huxtable with descriptive statistics stratified by levels of the outcome.
#' @seealso \code{\link[gtsummary]{tbl_summary}}
#' @examples
#' require(dplyr, quietly = TRUE)
#'
#' #' data(Oncho)
#'
#' ## A two by two contingency table:
#' Oncho |>
#'   select(mf, area) |>
#'   cross_tbl(by = "mf", bold = TRUE) |>
#'   theme_pubh(2)
#'
#' ## Reporting prevalence:
#' Oncho |>
#'   select(mf, area) |>
#'   cross_tbl(by = "area", bold = TRUE) |>
#'   theme_pubh(2)
#'
#' ## Descriptive statistics for all variables in the \code{Oncho} data set except \code{id}.
#' Oncho |>
#'   select(-id) |>
#'   cross_tbl(by = "mf", bold = TRUE) |>
#'   theme_pubh(2)
#' @export
cross_tbl <- function(data, by, head_label = " ",
                      bold = TRUE, show_total = TRUE, p_val = FALSE,
                      pad = 3, method = 2, ...) {
  if (method == 1) {
    tbl <- data |>
      gtsummary::tbl_summary(
        by = tidyselect::all_of(by),
        statistic = list(all_continuous() ~ "{mean} ({sd})"),
        ...
      ) |>
      gtsummary::modify_spanning_header(all_stat_cols() ~ paste0(
        "**",
        sjlabelled::get_label(data[[by]]),
        "**"
      ))
  } else {
    tbl <- data |>
      gtsummary::tbl_summary(
        by = tidyselect::all_of(by),
        ...
      ) |>
      gtsummary::modify_spanning_header(all_stat_cols() ~ paste0(
        "**",
        sjlabelled::get_label(data[[by]]),
        "**"
      ))
  }

  if (show_total == FALSE) {
    tbl <- tbl
  } else {
    tbl <- tbl |>
      gtsummary::add_overall(last = TRUE)
  }

  tbl <- tbl |>
    gtsummary::modify_header(label ~ head_label) |>
    gtsummary::modify_footnote(everything() ~ NA)

  if (bold == FALSE) {
    tbl <- tbl
  } else {
    tbl <- tbl |> gtsummary::bold_labels()
  }

  if (p_val == TRUE) {
    tbl <- tbl |> gtsummary::add_p()
  } else {
    tbl <- tbl
  }

  tbl <- tbl |>
    gtsummary::as_hux_table() |>
    huxtable::set_align(huxtable::everywhere, -1, "right") |>
    huxtable::set_top_padding(pad) |>
    huxtable::set_bottom_padding(pad) |>
    huxtable::set_align(1, huxtable::everywhere, "center")

  tbl
}

