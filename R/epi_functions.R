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
#' Mortality <- gl(2, 2, labels=c("No", "Yes"))
#' Calcium <- gl(2, 1, 4, labels=c("No", "Yes"))
#' anyca <- data.frame(Freq, Mortality, Calcium)
#' anyca
#' anyca.exp <- expand_df(anyca)
#' with(anyca.exp, table(Calcium, Mortality))
expand_df <- function(aggregate.data, index.var = "Freq", retain.freq = FALSE)
{
  output <- NULL
  for (i in 1:nrow(aggregate.data)) {
    if (retain.freq) {
      output <- rbind(output, aggregate.data[rep(i, aggregate.data[,
                                                                   which(names(aggregate.data) == index.var)][i]),
                                             ])
    } else {
      output <- rbind(output, aggregate.data[rep(i, aggregate.data[,
                                                                   which(names(aggregate.data) == index.var)][i]),
                                             ][, -which(names(aggregate.data) == index.var)])
    }
  }
  data.frame(output, row.names = 1:nrow(output))
}

#' Mantel-Haenszel odds ratio.
#'
#' \code{mhor} computes odds ratios by levels of the stratum variable as well as the Mantel-Haenszel
#' pooled odds ratio. The test for effect modification (test for interaction) is also displayed.
#'
#' @seealso \link{mh}
#' @param object When chaining, this holds an object produced in the earlier portions of the chain. Most users can safely ignore this argument. See details and examples.
#' @param formula A formula with shape: outcome ~ stratum/exposure.
#' @param data A data frame containing the variables used in \code{formula}.
#' @return Odds ratios with 95% confidence intervals on the effect of \code{exposure} on
#'   \code{outcome} by levels of \code{stratum}. The Mantel-Haenszel pooled OR and the test
#'   for effect modification is also reported.
#' @examples
#' data(oswego, package = "epitools")
#' require(dplyr)
#' require(sjlabelled)
#' oswego <- oswego %>%
#'   mutate(
#'     ill = factor(ill, labels = c("No", "Yes")),
#'     sex = factor(sex, labels = c("Female", "Male")),
#'     chocolate.ice.cream = factor(chocolate.ice.cream, labels = c("No", "Yes"))
#'   ) %>%
#'   var_labels(
#'     ill = "Developed illness",
#'     sex = "Sex",
#'     chocolate.ice.cream = "Consumed chocolate ice cream"
#'   )
#'
#' require(moonBook)
#' mytable(ill ~ sex + chocolate.ice.cream, data = oswego, show.total = TRUE)
#'
#' oswego %>%
#'   mhor(ill ~ sex/chocolate.ice.cream)
mhor <- function(object = NULL, formula = NULL, data = NULL)
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
    "OR" = round(model.ci[n1:nr2,5], 2),
    "Lower CI" = round(model.ci[n1:nr2,6], 2),
    "Upper CI" = round(model.ci[n1:nr2,7], 2),
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
  if(model2.aov[[3,3]] < 0.001) {
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
prop_or <- function(p2, or)
{
	p1 <- 1-p2
	por <- (or*p2)/(p1 + or*p2)
	por
}

#' Function to calculate OR using Wald CI, and plot trend.
#'
#' \code{odds_trend} calculates the odds ratio with confidence intervals (Wald) for different levels
#' (three or more) of the exposure variable, constructs the corresponding plot and calculates if the trend is
#' significant or not.
#'
#' @details \code{odds_trend} is a wrap function that calls \code{oddsratio} from package \code{epitools}.
#' @param formula A formula with shape: outcome ~ exposure.
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param angle Angle of for the x labels (default = 45).
#' @param hjust Horizontal adjustment for x labels (default = 1).
#' @param ... Passes optional arguments to \code{oddsratio}.
#' @return A list with components \code{df} a data frame with the results and \code{fig} corresponding plot.
#' @seealso \code{\link[epitools]{oddsratio}}.
#' @examples
#' ## A cross-sectional study looked at the association between obesity and a biopsy resulting
#' ## from mammography screening.
#'
#' Freq <- c(3441, 34, 39137, 519, 20509, 280, 12149, 196, 11882, 199)
#' Biopsy <- gl(2, 1, 10, labels = c("No", "Yes"))
#' Weight <- gl(5, 2, 10, labels = c("Underweight", "Normal", "Over (11-24%)",
#'              "Over (25-39%)", "Over (> 39%)"))
#' breast <- data.frame(Freq, Biopsy, Weight)
#' breast
#'
#' breast <- expand_df(breast)
#' require(sjlabelled)
#'
#' breast = var_labels(breast,
#'   Weight = 'Weight group'
#'   )
#'
#' odds_trend(Biopsy ~ Weight, data = breast)$df
#'
#' odds_trend(Biopsy ~ Weight, data = breast)$fig
odds_trend <- function(formula, data, angle = 45,
                       hjust = 1, ...)
{
  vars <- all.vars(formula)
  x <- vars[2]
  outcome <- data[[vars[1]]]
  exposure <- data[[vars[2]]]
  orwald <- epitools::oddsratio(exposure, outcome, ...)
  n <- nrow(orwald$measure)
  or.df <- data.frame(x = 1:n, round(orwald$measure, 2),
                      round_pval(orwald$p.value[, 3]),
                      round_pval(orwald$p.value[, 2]))
  names(or.df)[5:6] <- c("chi.square", "fisher.exact")
  names(or.df)[1:2] <- c("Exposure", "OR")
  nam <- row.names(or.df)
  or.df$Exposure <- factor(or.df$Exposure, labels = nam)
  or.df <- data.frame(or.df, row.names = NULL)
  df <- or.df
  fit = glm(formula, data = data, binomial)
  fig = sjPlot::plot_model(fit, type = 'pred', terms = x,
                           title = '', dot.size = 1) %>%
    gf_theme(axis.text.x = ggplot2::element_text(angle = angle, hjust = hjust))
  res <- list(df = df, fig = fig)
  res
}

#' Diagnostic tests from variables.
#'
#' \code{diag_test} is a wrap function that calls \code{epi.tests} from package \code{epiR}.
#' It computes sensitivity, specificity and other statistics related with screening tests.
#'
#' @param object When chaining, this holds an object produced in the earlier portions of the chain. Most users can safely ignore this argument. See details and examples.
#' @param formula A formula with shape: outcome ~ predictor (see details).
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param ... Further arguments passed to \code{epi.tests}.
#' @details For the \code{formula}, the outcome is the gold standard and the explanatory variable is the new (screening) test. See examples.
#' @seealso \code{\link[epiR]{epi.tests}}.
#' @examples
#' ## We compare the use of lung’s X-rays on the screening of TB against the gold standard test.
#' Freq <- c(1739, 8, 51, 22)
#' BCG <- gl(2, 1, 4, labels=c("Negative", "Positive"))
#' Xray <- gl(2, 2, labels=c("Negative", "Positive"))
#' tb <- data.frame(Freq, BCG, Xray)
#' tb <- expand_df(tb)
#'
#' tb %>%
#'   diag_test(BCG ~ Xray)
diag_test <- function(object = NULL, formula = NULL, data = NULL, ...)
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
  outcome <- vars[1]
  exposure <- vars[2]
  dat  <- epitools::epitable(data[[exposure]], data[[outcome]], rev='both')
  print(epiR::epi.tests(dat, ...))
  cat("\n")
}

#' Diagnostic tests from direct input.
#'
#' \code{diag_test2} is a wrap that calls \code{epi.tests} from package \code{epiR}.
#' It computes sensitivity, specificity and other statistics related with screening tests.
#'
#' @details \code{diag.test} uses direct input variables.
#'
#' @param aa Number of cases where both screening test and the gold standard are positive.
#' @param bb Number of cases where screening test is positive but gold standard is negative.
#' @param cc Number of cases where screening test is negative but gold standard is positive.
#' @param dd Number of cases where both screening test and the gold standard are negative.
#' @seealso \code{\link[epiR]{epi.tests}}.
#' @examples
#' ## We compare the use of lung’s X-rays on the screening of TB against the gold standard test.
#' diag_test2(22, 51, 8, 1739)
diag_test2 <- function(aa, bb, cc, dd)
{
	cat("\n")
	dat  <- as.table(matrix(c(aa, bb, cc, dd), nrow = 2, byrow = TRUE))
	print(epiR::epi.tests(dat))
	cat("\n")
}

#' Measures of association from two by two contingency tables (formula).
#'
#' \code{contingency} is a wrap that calls \code{epi.2by2} from package \code{epiR}.
#'
#' @param object When chaining, this holds an object produced in the earlier portions of the chain. Most users can safely ignore this argument. See details and examples.
#' @param formula A formula with shape: outcome ~ exposure.
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param method A character string with options: "cohort.count", "cohort.time", "case.control", or "cross.sectional".
#' @param ... Further arguments passed to \code{epi.2by2}.
#' @details \code{contingency} uses a formula as a way to input variables.
#' @details \code{contingency} displays the contingency table as a way for the user to check that the reference levels
#' in the categorical variables (outcome and exposure) are correct. Then displays measures of association
#' (table from \code{epi.2by2}). It also reports either chi-squared
#' test or exact Fisher's test;
#' \code{contingency} checks which one of the tests two is appropriate.
#' @seealso \code{\link[epiR]{epi.2by2}}.
#' @examples
#' ## A case-control study on the effect of alcohol on oesophageal cancer.
#' Freq <- c(386, 29, 389, 171)
#' status <- gl(2, 1, 4, labels=c("Control", "Case"))
#' alcohol <- gl(2, 2, labels=c("0-39", "40+"))
#' cancer <- data.frame(Freq, status, alcohol)
#' cancer <- expand_df(cancer)
#' contingency(status ~ alcohol, data = cancer, method = "case.control")
#'
#' data(Oncho)
#' require(moonBook)
#'
#' mytable(mf ~ area, data = Oncho, show.total = TRUE)
#'
#' Oncho %>%
#'   contingency(mf ~ area)
contingency <- function(object = NULL, formula = NULL, data = NULL, method="cohort.count", ...)
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
  outcome <- vars[1]
  exposure <- vars[2]
  dat  <- epitools::epitable(data[[exposure]], data[[outcome]], rev='both')
  print(dat)
  cat("\n")
  dat2  <- as.table(matrix(c(as.numeric(dat[1]), as.numeric(dat[3]), as.numeric(dat[2]), as.numeric(dat[4])), nrow = 2, byrow = TRUE))
  dat.epi <- epiR::epi.2by2(dat2, method = method, ...)
  print(dat.epi)
  cat("\n")
  if(chisq.fisher(dat) == 1) print(fisher.test(dat))
  else print(chisq.test(dat))
}

#' Measures of association from two by two contingency tables (direct input).
#'
#' \code{contingency2} is a wrap that calls \code{epi.2by2} from package \code{epiR}.
#'
#' @param aa Number of cases where both exposure and outcome are present.
#' @param bb Number of cases where exposure is present but outcome is absent.
#' @param cc Number of cases where exposure is absent but outcome is present.
#' @param dd Number of cases where both exposure and outcome are absent.
#' @param ... Further arguments passed to \code{\link[epiR]{epi.2by2}}.
#' @seealso \code{\link[epiR]{epi.2by2}}.
#' @examples
#' ## A case-control study on the effect of alcohol on oesophageal cancer.
#' Freq <- c(386, 29, 389, 171)
#' status <- gl(2, 1, 4, labels=c("Control", "Case"))
#' alcohol <- gl(2, 2, labels=c("0-39", "40+"))
#' cancer <- data.frame(Freq, status, alcohol)
#' cancer <- expand_df(cancer)
#'
#' contingency2(171, 389, 29, 386, method = "case.control")
contingency2 <- function(aa, bb, cc, dd, ...)
{
	dat  <- as.table(matrix(c(aa, bb, cc, dd), nrow = 2, byrow = TRUE))
	colnames(dat) <- c("Yes", "No")
	rownames(dat) <- c("Yes", "No")
	cat("\n")
	print(dat)
	cat("\n")
	dat.epi <- epiR::epi.2by2(dat, ...)
	print(dat.epi)
	cat("\n")
	if(chisq.fisher(dat) == 1) print(fisher.test(dat))
		else print(chisq.test(dat))
}


#' Internal test for chi-squared assumption.Fisher (2 by 2). If results = T, it fails
#'
#' \code{chisq.fisher} is an internal function called by \code{contingency} and \code{contingency2} that uses the Fisher exact test if results from the assumptions for the chi-squared test fail.
#' @param tab A numeric two by two table.
chisq.fisher <- function(tab)
{
	n <- sum(tab)
	sr <- rowSums(tab)
	sc <- colSums(tab)
	elem <- outer(sr, sc, "*")/n
	res <- any(elem < 5)*1
	res
}
