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
rel_dis <- function(x)
{
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
reference_range <- function(avg, std)
{
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
#' estat(~ IgM, data = Ab)
#' geo_mean(IgM)
geo_mean <- function(x)
{
  positive.check <- which(x <= 0)
  if (length(positive.check) >= 1) {
    print('Some observations are not positive')
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
#' estat(~ IgM, data = Ab)
#' harm_mean(IgM)
harm_mean <- function(x)
{
  zero.check <- which(x == 0)
  if (length(zero.check) >= 1) {
    print('Some observations are equal to zero')
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
bst <- function (x, stat = "mean", n = 1000, CI = 95, digits = 2)
{
  xmeans <- numeric(n)
  if (stat == "median") {
    for (i in 1:n) xmeans[i] <- median(sample(x, replace = TRUE),
                                       na.rm = TRUE)
    estimate <- median(x, na.rm = TRUE)
  }
  else if (stat == "gmean") {
    for (i in 1:n) xmeans[i] <- geo_mean(sample(x, replace = TRUE))
    estimate <- geo_mean(x)
  }
  else if (stat == "mean") {
    for (i in 1:n) xmeans[i] <- mean(sample(x, replace = TRUE),
                                     na.rm = TRUE)
    estimate <- mean(x, na.rm = TRUE)
  }
  else {
    stop("Statistic not available")
  }
  alpha <- 1 - (CI/100)
  tail <- alpha/2
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
stats_quotes <- function(x, data2, digits = 2)
{
  n <- sum(!is.na(data2[[x]]))
  x.min <- min(data2[[x]], na.rm = TRUE)
  x.max <- max(data2[[x]], na.rm = TRUE)
  avg <- mean(data2[[x]], na.rm = TRUE)
  med <- median(data2[[x]], na.rm = TRUE)
  std <- sd(data2[[x]], na.rm = TRUE)
  cv <- std / avg
  res <- data.frame(n, x.min, x.max, avg, med, std, cv)
  names(res) <- c("N", "Min.", "Max.", "Mean", "Median", "SD", "CV")
  res <- round(res, digits=digits)
  res
}

#' Cross-tabulation.
#'
#' \code{cross_tab} is a wrapper to functions from package \code{moonBook} to construct tables of descriptive statistics stratified by levels of a categorical outcome.
#'
#' @param object When chaining, this holds an object produced in the earlier portions of the chain. Most users can safely ignore this argument. See details and examples.
#' @param formula A formula with shape: \code{y ~ x}, where \code{y} is a categorical outcome and \code{x} is the explanatory variable or a set of explanatory variables (see Details and Examples).
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param label A character, label to be used for the outcome (for non-labelled data).
#' @param show.total Logical, show column with totals?
#' @param p_val Logical, show p-values?
#' @param ... Additional arguments passed to \code{\link[moonBook]{mytable_sub}}.
#' @details Function \code{cross_tab} is a relatively simple wrapper to function \code{mytable} of package \code{moonBook}. Its main purpose is to construct contingency tables but it can also be used to report a table with descriptives for all variables as long as they are still stratified by the outcome. Please see examples to see how to list explanatory variables. For categorical explanatory variables, the function reports column percentages. If data is labelled with \code{sjlabelled}, the label of the outcome (dependent) variable is used to name the outcome; this name can be changed with argument \code{label}.
#' @return A huxtable with descriptive statistics stratified by levels of the outcome.
#' @seealso \code{\link[moonBook]{mytable}}
#' @examples
#' data(Oncho)
#'
#' ## A two by two contingency table:
#' Oncho %>%
#'   cross_tab(mf ~ area)
#'
#' ## Reporting prevalence:
#' Oncho %>%
#'   cross_tab(area ~ mf)
#'
#' ## Contingency table for both sex and area of residence:
#' Oncho %>%
#'   cross_tab(mf ~ sex + area, p_val = TRUE)
#'
#' ## Descriptive statistics for all variables in the \code{Oncho} data set except \code{id}.
#' require(dplyr, quietly = TRUE)
#' Oncho %>%
#'   select(- id) %>%
#'   cross_tab(mf ~ .)
cross_tab <- function(object = NULL, formula = NULL, data = NULL,
                      label = NULL, show.total = TRUE, p_val = FALSE, ...)
{
  if (inherits(object, "formula")) {
    formula <- object
    object <- NULL
  }
  if (inherits(object, "data.frame")) {
    data <- object
    object <- NULL
  }
  vars <- all.vars(formula)
  dependent <- vars[1]
  explanatory <- vars[-1]
  nv <- length(explanatory)
  outcome <- data[[dependent]]
  k <- length(levels(outcome))
  nl <- k + 2
  nm <- k + 1
  nl2 <- nl + 1
  if (is.null(sjlabelled::get_label(outcome)) == FALSE & is.null(label)) {
    lab <- sjlabelled::get_label(outcome)
  }
  else {
    lab <- ifelse(is.null(label), dependent, label)
  }
  tbl <- moonBook::mytable(formula, data = data, show.total = show.total, ...) %>%
    moonBook::mytable2df()
  if(p_val == TRUE) {
    tvl <- tbl
  } else{
    tvl <- tbl %>% dplyr::select(-ncol(tbl))
  }
  names(tvl)[1] <- " "
  res <- tvl %>%
    huxtable::as_hux() %>% huxtable::set_align(huxtable::everywhere, 2:nl, "right") %>%
    huxtable::insert_row(c(" ", lab), fill = " ") %>%
    huxtable::merge_cells(1, 2:nm) %>%
    huxtable::set_header_rows(1, TRUE) %>%
    huxtable::set_align(1:3, huxtable::everywhere, "center")
  if(nv > 1){
    res = res
  } else {
    res = res %>% huxtable::set_bold(4, 1)
  }
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
#' @seealso \code{\link{summary}}, \code{\link[moonBook]{mytable}}.
#' @examples
#' data(kfm, package = "ISwR")
#' require(sjlabelled, quietly = TRUE)
#' kfm = kfm %>%
#' var_labels(
#'   dl.milk = 'Breast-milk intake (dl/day)',
#'   sex = 'Sex',
#'   weight = 'Child weight (kg)',
#'   ml.suppl = 'Milk substitute (ml/day)',
#'   mat.weight = 'Maternal weight (kg)',
#'   mat.height = 'Maternal height (cm)'
#'  )
#'
#' kfm %>%
#'   estat(~ dl.milk)
#'
#' estat(~ dl.milk|sex, data = kfm)
#'
#' kfm %>%
#'   estat(~ weight|sex)
estat <- function (object = NULL, formula = NULL, data = NULL, digits = 2, label = NULL)
{
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
  if (!is.numeric(outcome))
    stop(y, " must be a numerical variable")
  if (is.null(sjlabelled::get_label(outcome)) == FALSE & is.null(label)){
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
  }
  else {
    x <- vars[2]
    exposure <- data[[x]]
    if (!is.factor(exposure))
      stop(x, " must be a factor")
    lev <- levels(exposure)
    nl <- length(lev)
    res <- numeric(nl * 7)
    dim(res) <- c(nl, 7)
    for (i in 1:nl) {
      res[i, ] <- as.numeric(stats_quotes(y,
                                          data2 = subset(data, exposure == lev[i])))
    }
    res <- as.data.frame(res)
    names(res) <- c("N", "Min.", "Max.", "Mean", "Median",
                    "SD", "CV")
    res <- round(res, digits = digits)
    res <- data.frame(var = lev, res)
    names(res)[1] <- x
    var <- character(nl)
    var[1] <- lab
    res <- data.frame(var, res)
    names(res)[1] <- ""
    names(res)[2] <- ifelse(is.null(sjlabelled::get_label(exposure)) == FALSE,
                            sjlabelled::get_label(exposure), x)
    res
  }
}
