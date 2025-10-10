# Epidemiology-related functions

#' Expand a data frame.
#'
#' \code{expand_df} expands a data frame by a vector of frequencies.
#'
#' This is a generic function that resembles weighted frequencies in other statistical packages
#' (for example, Stata). \code{expand.df} was adapted from a function developed by deprecated package
#' \code{epicalc} (now package \code{epiDisplay}).
#'
#' @param aggregate.data A data frame.
#' @param index.var A numerical variable with the frequencies (counts).
#' @param retain.freq Logical expression indicating if frequencies should be kept.
#' @return An expanded data frame with replicates given by the frequencies.
#' @examples
#' Freq <- c(5032, 5095, 41, 204)
#' Mortality <- gl(2, 2, labels = c("No", "Yes"))
#' Calcium <- gl(2, 1, 4, labels = c("No", "Yes"))
#' anyca <- data.frame(Freq, Mortality, Calcium)
#' anyca
#' anyca.exp <- expand_df(anyca)
#' with(anyca.exp, table(Calcium, Mortality))
#' @export
expand_df <- function(aggregate.data, index.var = "Freq", retain.freq = FALSE) {
  output <- NULL
  for (i in 1:nrow(aggregate.data)) {
    if (retain.freq) {
      output <- rbind(output, aggregate.data[rep(i, aggregate.data[
        ,
        which(names(aggregate.data) == index.var)
      ][i]), ])
    } else {
      output <- rbind(output, aggregate.data[rep(i, aggregate.data[
        ,
        which(names(aggregate.data) == index.var)
      ][i]), ][, -which(names(aggregate.data) == index.var)])
    }
  }
  data.frame(output, row.names = 1:nrow(output))
}

#' Mantel-Haenszel odds ratio.
#'
#' \code{mhor} computes odds ratios by levels of the stratum variable as well as the Mantel-Haenszel
#' pooled odds ratio. The test for effect modification (test for interaction) is also displayed.
#'
#' @param object When chaining, this holds an object produced in the earlier portions of the chain. Most users can safely ignore this argument. See details and examples.
#' @param formula A formula with shape: outcome ~ stratum/exposure.
#' @param data A data frame containing the variables used in \code{formula}.
#' @return Odds ratios with 95% confidence intervals on the effect of \code{exposure} on
#'   \code{outcome} by levels of \code{stratum}. The Mantel-Haenszel pooled OR and the test
#'   for effect modification is also reported.
#' @examples
#' data(oswego, package = "epitools")
#' require(dplyr, quietly = TRUE)
#' require(sjlabelled, quietly = TRUE)
#' oswego <- oswego |>
#'   mutate(
#'     ill = factor(ill, labels = c("No", "Yes")),
#'     sex = factor(sex, labels = c("Female", "Male")),
#'     chocolate.ice.cream = factor(chocolate.ice.cream, labels = c("No", "Yes"))
#'   ) |>
#'   var_labels(
#'     ill = "Developed illness",
#'     sex = "Sex",
#'     chocolate.ice.cream = "Consumed chocolate ice cream"
#'   )
#'
#' oswego |>
#'   mhor(ill ~ sex / chocolate.ice.cream)
#' @export
mhor <- function(object = NULL, formula = NULL, data = NULL) {
  if (inherits(object, "formula")) {
    formula <- object
    object <- NULL
  }
  if (inherits(object, "data.frame")) {
    data <- object
    object <- NULL
  }
  vars <- all.vars(formula)
  outcome <- vars[1]
  exposure <- vars[3]
  stratum <- vars[2]
  model <- glm(formula, data = data, family = binomial)
  model2 <- glm(data[[outcome]] ~ data[[stratum]] * data[[exposure]], family = binomial)
  model.ci <- Epi::ci.lin(model, Exp = TRUE)
  model2.aov <- car::Anova(model2)
  nr2 <- nrow(model.ci)
  n1 <- length(levels(data[[stratum]])) + 1
  modelci <- data.frame(
    "OR" = round(model.ci[n1:nr2, 5], 2),
    "Lower CI" = round(model.ci[n1:nr2, 6], 2),
    "Upper CI" = round(model.ci[n1:nr2, 7], 2),
    "p" = round_pval(model.ci[n1:nr2, 4])
  )
  names(modelci)[4] <- "Pr(>|z|)"
  mht <- mantelhaen.test(data[[outcome]], data[[exposure]], data[[stratum]])
  res <- data.frame(round(mht$estimate[[1]], 2), round(mht$conf.int[1], 2), round(mht$conf.int[2], 2), round_pval(mht$p.value))
  names(res) <- c("Common OR", "Lower CI", "Upper CI", "Pr(>|z|)")
  rownames(res) <- "Cochran-Mantel-Haenszel: "
  cat("\n")
  print(modelci)
  cat("\n")
  print(res)
  cat("\n")
  if (model2.aov[[3, 3]] < 0.001) {
    cat("Test for effect modification (interaction): p ", round_pval(model2.aov[[3, 3]]), "\n", "\n")
  } else {
    cat("Test for effect modification (interaction): p = ", round(model2.aov[[3, 3]], 4), "\n", "\n")
  }
}

#' Proportion, p1 from proportion p2 and OR.
#'
#' \code{prop_or} is a simple function to calculate a proportion, from another proportion and the odds
#' ratio between them.
#'
#' @param p2 The value of a proportion in the unexposed group (p2).
#' @param or The odds ratio of p1/p2.
#' @return \code{p1}, the proportion in the exposed group (p1).
#' @examples
#' flu <- matrix(c(20, 80, 220, 140), nrow = 2)
#' colnames(flu) <- c("Yes", "No")
#' rownames(flu) <- c("Vaccine", "Placebo")
#' flu
#'
#' or <- (20 * 140) / (80 * 220)
#' p2 <- 80 / 220
#' prop_or(p2 = p2, or = or)
#' 20 / 240
#' @export
prop_or <- function(p2, or) {
  p1 <- 1 - p2
  por <- (or * p2) / (p1 + or * p2)
  por
}

#' Function to calculate OR using Wald CI, and plot trend.
#'
#' \code{odds_trend} calculates the odds ratio with confidence intervals (Wald) for different levels.
#' (three or more) of the exposure variable.
#' @details \code{odds_trend} is a wrap function that calls \code{oddsratio} from package \code{epitools}.
#' @param data A data frame.
#' @param x A variable in the data frame.
#' @param y A variable in the data frame.
#' @param method Method for calculating confidence interval around odds ratio.
#' @param ... Passes optional arguments to \code{oddsratio}.
#' @details Additional methods for confidence intervals include: \code{"midp"}, \code{"fisher"}, and \code{"small"}.
#' @return A data frame.
#' @seealso \code{\link[epitools]{oddsratio}}.
#' @examples
#' ## A cross-sectional study looked at the association between obesity and a biopsy resulting
#' ## from mammography screening.
#'
#' Freq <- c(3441, 34, 39137, 519, 20509, 280, 12149, 196, 11882, 199)
#' Biopsy <- gl(2, 1, 10, labels = c("No", "Yes"))
#' Weight <- gl(5, 2, 10, labels = c(
#'   "Underweight", "Normal", "Over (11-24%)",
#'   "Over (25-39%)", "Over (> 39%)"
#' ))
#' breast <- data.frame(Freq, Biopsy, Weight)
#' breast
#'
#' breast <- expand_df(breast)
#' require(sjlabelled, quietly = TRUE)
#'
#' breast <- var_labels(breast,
#'   Weight = "Weight group"
#' )
#'
#' breast |>
#'   odds_trend(x = Weight, y = Biopsy)
#'
#' breast |>
#'   odds_trend(x = Weight, y = Biopsy) |>
#'   effect_plot(x = OR, y = Exposure)
#' @export
odds_trend <- function(
    data, x, y,
    method = "wald", ...
) {
  exposure <- data[[deparse(substitute(x))]]
  outcome <- data[[deparse(substitute(y))]]
  orwald <- epitools::oddsratio(exposure, outcome, method = method, ...)
  n <- nrow(orwald$measure)
  or.df <- data.frame(
    x = 1:n, round(orwald$measure, 2),
    round_pval(orwald$p.value[, 3]),
    round_pval(orwald$p.value[, 2])
  )
  names(or.df)[5:6] <- c("chi.square", "fisher.exact")
  names(or.df)[1:2] <- c("Exposure", "OR")
  names(or.df)[3:4] <- c("CI_low", "CI_high")
  nam <- row.names(or.df)
  or.df$Exposure <- factor(or.df$Exposure, labels = nam)
  or.df <- data.frame(or.df, row.names = NULL)
  or.df
}
