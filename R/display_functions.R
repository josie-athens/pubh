# Functions to display coefficients

#' Rounding p-values.
#'
#' \code{round_pval} is an internal function called by \code{glm_coef} to round p-values from model coefficients.
#' @param pval vector of p-values, numeric.
round_pval <- function(pval)
{
	res <- ifelse(pval < 0.001, "< 0.001", round(pval, 3))
}

#' Table of coefficients from generalised linear models.
#'
#' \code{glm_coef} displays estimates with confidence intervals and p-values from generalised linear models (see Details).
#'
#' @details \code{glm_coef} recognises objects (models) from the following classes: \code{clm}, \code{clogit},
#' \code{coxph}, \code{gee}, \code{glm}, \code{glmerMod}, \code{lm}, \code{lme}, \code{multinom}, \code{negbin},
#' \code{polr} and \code{surveg}
#' @details For models from logistic regression (including conditional logistic, ordinal and multinomial),
#' Poisson or survival analysis, coefficient estimates and corresponding confidence intervals are
#' automatically exponentiated (back-transformed).
#' @details By default, \code{glm_coef} uses robust standard errors for calculating confidence intervals.
#' @details \code{glm_coef} can display two different data frames depending on the option of \code{type},
#' for type \code{type = "cond"} (the default), \code{glm_coef} displays the standard table of coefficients
#' with confidence intervals and p-values; for \code{type = "ext"}, \code{glm_coef} displays each number
#' in a different column and includes standard errors.
#' @details Please read the Vignette on Regression for more details.
#' @param model A model from any of the classes listed in the details section.
#' @param digits A scalar, number of digits for rounding the results (default = 2).
#' @param alpha Significant level (default = 0.05) used to calculate confidence intervals.
#' @param labels An optional character vector with the names of the coefficients (including intercept).
#' @param se_rob Logical, should robust errors be used to calculate confidence intervals? (default = TRUE).
#' @param type Character, either "cond" (condensed) or "ext" (extended). See details.
#' @param exp_norm Logical, should estimates and confidence intervals should be exponentiated? (for family == "gaussian").
#' @return A data frame with estimates, confidence intervals and p-values from \code{glm} objects.
#' @examples
#' ## Continuous outcome.
#' data(birthwt, package = "MASS")
#' require(dplyr)
#' birthwt <- mutate(birthwt,
#'   smoke = factor(smoke, labels = c("Non-smoker", "Smoker")),
#'   Race = factor(race > 1, labels = c("White", "Non-white")))
#'
#' model_norm <- glm(bwt ~ smoke + race, data = birthwt)
#' glm_coef(model_norm)
#'
#' model_norm %>%
#'   glm_coef(labels=c("Constant", "Smoker vs Non-smoker", "Non-white vs White"))
#'
#' ## Logistic regression.
#' data(diet, package = "Epi")
#' model_binom <- glm(chd ~ fibre, data = diet, family = binomial)
#' model_binom %>%
#'   glm_coef(labels = c("Constant", "Fibre intake (g/day)"))
#'
#' model_binom %>%
#' glm_coef(labels = c("Constant", "Fibre intake (g/day)"), type = "ext")
#'
#' ## Poisson regression.
#' library(MASS)
#' data(quine)
#' levels(quine$Eth) <- list(White = "N", Aboriginal = "A")
#' levels(quine$Sex) <- list(Male = "M", Female = "F")
#' model_pois <- glm(Days ~ Eth + Sex + Age, family = poisson, data = quine)
#'
#' model_pois %>%
#'   glm_coef()
#'
#' deviance(model_pois) / df.residual(model_pois) # to check for overdispersion
#'
#' model_negbin <- glm.nb(Days ~ Eth + Sex + Age, data = quine)
#' unadj <- glm_coef(model_negbin,
#'                   labels=c("Constant",
#'                   "Race: Aboriginal/White",
#'                   "Sex: Female/Male",
#'                   "F1/Primary",
#'                   "F2/Primary",
#'                   "F3/Primary"))
#' unadj # Not-adjusted for multiple comparisons
#'
#' ## For more examples, please read the Vignette on Regression.
glm_coef <- function(model, digits = 2, alpha = 0.05, labels = NULL, se_rob = TRUE,
                      type = "cond", exp_norm = FALSE) {
  if (type != "cond" & type != "ext")
    stop("Option type not supported")
  mod <- summary(model)
  n <- nrow(mod)
  conf.int <- 1 - alpha
  zcrit <- qnorm((1 + conf.int)/2)
  if (class(model)[1] == "clm" | class(model)[1] == "polr") {
    int.low <- exp(mod$coefficients[, 1] - mod$coefficients[, 2] * zcrit)
    int.up <- exp(mod$coefficients[, 1] + mod$coefficients[, 2] * zcrit)
    estim <- mod$coefficients[, 1]
    coeff <- exp(estim)
    coeff.error <- mod$coefficients[, 2]
    zval <- mod$coefficients[, 1]/mod$coefficients[, 2]
    pval <- round_pval((1 - pnorm(abs(zval))) * 2)
    if (type == "cond") {
      coeff_ci <- paste0(round(coeff, digits), " (", round(int.low, digits),
                         ", ", round(int.up, digits), ")")
      coeff.df <- data.frame(coeff_ci, pval)
      names(coeff.df) <- c("Ordinal OR", "Pr(>|Z|)")
    } else {
      coeff_ci <- round(data.frame(coeff, int.low, int.up, coeff.error), digits)
      coeff.df <- data.frame(coeff_ci, pval)
      names(coeff.df) <- c("Ordinal OR", "Lower CI", "Upper CI", "Std. Error",
                           "Pr(>|Z|)")
    }
    if (is.null(labels)) {
      rownames(coeff.df) <- rownames(coeff.df)
    } else {
      rownames(coeff.df) <- labels
    }
    coeff.df
  } else if (class(model)[1] == "multinom") {
    n.mod <- length(rownames(mod$coefficients))
    if (type == "cond") {
      for (i in 1:n.mod) {
        est <- mod$coefficients[i, ]
        back <- exp(mod$coefficients[i, ])[-1]
        modse <- mod$standard.errors[i, ]
        low <- exp(mod$coefficients[i, ] - modse * zcrit)[-1]
        up <- exp(mod$coefficients[i, ] + modse * zcrit)[-1]
        zval <- est/modse
        pval <- round_pval((1 - pnorm(abs(zval))) * 2)[-1]
        back_ci <- paste0(round(back, digits), " (", round(low, digits),
                          ", ", round(up, digits), ")")
        df <- data.frame(back_ci, pval)
        colnames(df) <- c("Multinomial OR", "Pr(>|z|)")
        if (is.null(labels)) {
          rownames(df) <- colnames(mod$coefficients)[-1]
        } else {
          rownames(df) <- labels
        }
        rownames(df) <- colnames(mod$coefficients)[-1]
        cat("\n")
        print(paste(rownames(mod$coefficients)[i], " vs ", mod$lev[1]))
        print(df)
        cat("\n")
      }
    } else {
      for (i in 1:n.mod) {
        est <- mod$coefficients[i, ]
        back <- exp(mod$coefficients[i, ])
        back[1] <- NA
        modse <- mod$standard.errors[i, ]
        low <- exp(mod$coefficients[i, ] - modse * zcrit)
        low[1] <- mod$coefficients[i, 1] - modse[1] * zcrit
        up <- exp(mod$coefficients[i, ] + modse * zcrit)
        up[1] <- mod$coefficients[i, 1] + modse[1] * zcrit
        zval <- est/modse
        pval <- round_pval((1 - pnorm(abs(zval))) * 2)
        df <- round(cbind(back, low, up, zval), digits)
        df <- data.frame(df, pval)
        colnames(df) <- c("Multinomial OR", paste("lower", 100 - 100 * alpha,
                                                  "ci", sep = ""), paste("upper", 100 - 100 * alpha, "ci", sep = ""),
                          "z value", "Pr(>|z|)")
        if (is.null(labels)) {
          rownames(df) <- colnames(mod$coefficients)
        } else {
          rownames(df) <- labels
        }
        rownames(df) <- colnames(mod$coefficients)
        cat("\n")
        print(paste(rownames(mod$coefficients)[i], " vs ", mod$lev[1]))
        print(df)
        cat("\n")
      }
    }
  } else if (class(model)[1] == "lme") {
    mod <- summary(model)$tTable
    estim <- mod[, 1]
    low <- estim - mod[, 2] * zcrit
    up <- estim + mod[, 2] * zcrit
    stats <- mod[, 2:4]
    if (type == "cond") {
      estim_ci <- paste0(round(estim, digits), " (", round(low, digits), ", ",
                         round(up, digits), ")")
      out.df <- data.frame(estim_ci, round_pval(mod[, 5]))
      names(out.df) <- c("Coefficient", "Pr(>|t|)")
    } else {
      out.df <- round(data.frame(estim, low, up, stats), digits)
      out.df <- data.frame(out.df, round_pval(mod[, 5]))
      names(out.df) <- c("Coeff", "Lower CI", "Upper CI", "SE", "DF", "t value",
                         "Pr(>|t|)")
    }
    if (is.null(labels)) {
      rownames(out.df) <- rownames(out.df)
    } else {
      rownames(out.df) <- labels
    }
    out.df
  } else if (class(model)[1] == "clogit") {
    mod.coef <- mod$coefficients
    or.low <- exp(mod.coef[, 1] - mod.coef[, 3] * zcrit)
    or.up <- exp(mod.coef[, 1] + mod.coef[, 3] * zcrit)
    OR <- mod.coef[, 2]
    if (type == "cond") {
      or_ci <- paste0(round(OR, digits), " (", round(or.low, digits), ", ",
                      round(or.up, digits), ")")
      out.df <- data.frame(or_ci, round_pval(mod.coef[, 5]))
      colnames(out.df) <- c("Odds ratio", "Pr(>|z|)")
    } else {
      cis <- data.frame(OR, mod.coef[, 2], or.low, or.up)
      cis <- round(cis, digits)
      out.df <- data.frame(cis, round_pval(mod.coef[, 5]))
      colnames(out.df) <- c("OR", "Std. Error", "Lower CI", "Upper CI", "Pr(>|z|)")
    }
    if (is.null(labels)) {
      rownames(out.df) <- rownames(out.df)
    } else {
      rownames(out.df) <- labels
    }
    out.df
  } else if (class(model)[1] == "survreg") {
    mod.coef <- mod$table
    hc <- sandwich::vcovHAC(model)
    ct <- lmtest::coeftest(model, vcov = hc)
    np <- nrow(mod.coef)
    if (se_rob == TRUE) {
      mod.coef <- ct[, ]
    } else {
      mod.coef <- mod$table
    }
    Estimate <- mod.coef[, 1]
    OR <- exp(Estimate)
    or.low <- exp(Estimate - mod.coef[, 2] * zcrit)
    or.up <- exp(Estimate + mod.coef[, 2] * zcrit)
    if (type == "cond") {
      or_ci <- paste0(round(OR, digits), " (", round(or.low, digits), ", ",
                      round(or.up, digits), ")")
      out.df <- data.frame(or_ci, round_pval(mod.coef[, 4]))
      out.df <- out.df[-1, ]
      colnames(out.df) <- c("Survival time ratio", "Pr(>|z|)")
    } else {
      cis <- data.frame(OR, mod.coef[, 2], or.low, or.up)
      cis <- round(cis, digits)
      out.df <- data.frame(cis, round_pval(mod.coef[, 4]))
      out.df <- out.df[-1, ]
      colnames(out.df) <- c("Survival time ratio", "Std. Error", "Lower CI",
                            "Upper CI", "Pr(>|z|)")
    }
    if (is.null(labels)) {
      labs <- rownames(out.df)
      rownames(out.df) <- gsub("Log(scale)", replacement = "Scale", labs, fixed = TRUE)
    } else {
      rownames(out.df) <- labels
    }
    out.df
  } else if (class(model)[1] == "coxph") {
    mod.coef <- mod$coefficients
    or.low <- exp(mod.coef[, 1] - mod.coef[, 3] * zcrit)
    or.up <- exp(mod.coef[, 1] + mod.coef[, 3] * zcrit)
    OR <- mod.coef[, 2]
    if (type == "cond") {
      or_ci <- paste0(round(OR, digits), " (", round(or.low, digits), ", ",
                      round(or.up, digits), ")")
      out.df <- data.frame(or_ci, round_pval(mod.coef[, 5]))
      colnames(out.df) <- c("Hazard ratio", "Pr(>|z|)")
    } else {
      cis <- data.frame(OR, mod.coef[, 3], or.low, or.up)
      cis <- round(cis, digits)
      out.df <- data.frame(cis, round_pval(mod.coef[, 5]))
      colnames(out.df) <- c("Hazard ratio", "Std. Error", "Lower CI", "Upper CI",
                            "Pr(>|z|)")
    }
    if (is.null(labels)) {
      rownames(out.df) <- rownames(out.df)
    } else {
      rownames(out.df) <- labels
    }
    out.df
  } else if (class(model)[1] == "gee" & (family(model)$family == "binomial" | family(model)$family ==
                                         "poisson") | family(model)$family == "quasi") {
    mod <- mod$coefficients
    n <- nrow(mod)
    estim <- mod[, 1]
    coeff <- exp(mod[, 1])
    coeff[1] <- NA
    if (se_rob == TRUE) {
      segee <- mod[, 4]
      my.z <- mod[, 5]
    } else {
      segee <- mod[, 2]
      my.z <- mod[, 3]
    }
    low <- estim - segee * zcrit
    low[2:n] <- exp(low[2:n])
    up <- estim + segee * zcrit
    up[2:n] <- exp(up[2:n])
    pvals <- (1 - pnorm(abs(my.z))) * 2
    if (type == "cond") {
      coeff_ci <- paste0(round(coeff, digits), " (", round(low, digits), ", ",
                         round(up, digits), ")")
      out.df <- data.frame(coeff_ci, round_pval(pvals))
      names(out.df) <- c("Exp(Coeff)", "Pr(>|z|)")
    } else {
      out.df <- round(data.frame(estim, coeff, low, up, segee), digits)
      out.df <- data.frame(out.df, round_pval(pvals))
      names(out.df) <- c("Coeff", "Exp(Coeff)", "Lower CI", "Upper CI", "SE",
                         "Pr(>|z|)")
    }
    if (is.null(labels)) {
      rownames(out.df) <- rownames(out.df)
    } else {
      rownames(out.df) <- labels
    }
    out.df
  } else if (class(model)[1] == "gee" & family(model)$family == "gaussian") {
    mod <- mod$coefficients
    n <- nrow(mod)
    estim <- mod[, 1]
    if (se_rob == TRUE) {
      segee <- mod[, 4]
      my.z <- mod[, 5]
    } else {
      segee <- mod[, 2]
      my.z <- mod[, 3]
    }
    low <- estim - segee * zcrit
    up <- estim + segee * zcrit
    pvals <- (1 - pnorm(abs(my.z))) * 2
    if (type == "cond") {
      estim_ci <- paste0(round(estim, digits), " (", round(low, digits), ", ",
                         round(up, digits), ")")
      out.df <- data.frame(estim_ci, round_pval(pvals))
      names(out.df) <- c("Coefficient", "Pr(>|z|)")
    } else {
      out.df <- round(data.frame(estim, low, up, segee), digits)
      out.df <- data.frame(out.df, round_pval(pvals))
      names(out.df) <- c("Coeff", "Lower CI", "Upper CI", "SE", "Pr(>|z|)")
    }
    if (is.null(labels)) {
      rownames(out.df) <- rownames(out.df)
    } else {
      rownames(out.df) <- labels
    }
    out.df
  } else if (class(model)[1] == "glmerMod" & (family(model)$family == "binomial" |
                                              family(model)$family == "poisson")) {
    mod <- mod$coefficients
    n <- nrow(mod)
    estim <- mod[, 1]
    coeff <- exp(mod[, 1])
    coeff[1] <- NA
    low <- estim - mod[, 2] * zcrit
    low[2:n] <- exp(low[2:n])
    up <- estim + mod[, 2] * zcrit
    up[2:n] <- exp(up[2:n])
    if (type == "cond") {
      coeff_ci <- paste0(round(coeff, digits), " (", round(low, digits), ", ",
                         round(up, digits), ")")
      out.df <- data.frame(coeff_ci, round_pval(mod[, 4]))[-1, ]
      names(out.df) <- c("Exp(Coeff)", "Pr(>|z|)")
    } else {
      stats <- mod[, 2:3]
      out.df <- round(data.frame(estim, coeff, low, up, stats), digits)
      out.df <- data.frame(out.df, round_pval(mod[, 4]))
      names(out.df) <- c("Coeff", "Exp(Coeff)", "Lower CI", "Upper CI", "SE",
                         "z value", "Pr(>|z|)")
    }
    if (is.null(labels)) {
      rownames(out.df) <- rownames(out.df)
    } else {
      rownames(out.df) <- labels
    }
    out.df
  } else if (class(model)[1] == "glmerMod" & family(model)$family == "gaussian") {
    mod <- mod$coefficients
    n <- nrow(mod)
    estim <- mod[, 1]
    low <- estim - mod[, 2] * zcrit
    up <- estim + mod[, 2] * zcrit
    if (type == "cond") {
      estim_ci <- paste0(round(estim, digits), " (", round(low, digits), ", ",
                         round(up, digits), ")")
      out.df <- data.frame(estim_ci, round_pval(mod[, 4]))
      names(out.df) <- c("Coefficient", "Pr(>|z|)")
    } else {
      stats <- mod[, 2:3]
      out.df <- round(data.frame(estim, low, up, stats), digits)
      out.df <- data.frame(out.df, round_pval(mod[, 4]))
      names(out.df) <- c("Coeff", "Lower CI", "Upper CI", "SE", "z value",
                         "Pr(>|z|)")
    }
    if (is.null(labels)) {
      rownames(out.df) <- rownames(out.df)
    } else {
      rownames(out.df) <- labels
    }
    out.df
  } else {
    if (class(model)[1] != "glm" & class(model)[1] != "negbin" & class(model)[1] !=
        "lm")
      stop("Object's class or family error not supported.")
    if (class(model)[1] == "glm" | class(model)[1] == "negbin" |
        class(model)[1] == "lm") {
      tcrit <- qt((1 + conf.int)/2, df.residual(model))
      hc <- sandwich::vcovHAC(model)
      ct <- lmtest::coeftest(model, vcov = hc)
      if (se_rob == TRUE) {
        mod.coef <- ct[, ]
      } else {
        mod.coef <- mod$coefficients
      }
      or.low <- exp(mod.coef[, 1] - mod.coef[, 2] * zcrit)
      or.up <- exp(mod.coef[, 1] + mod.coef[, 2] * zcrit)
      ci.low <- mod.coef[, 1] - mod.coef[, 2] * tcrit
      ci.up <- mod.coef[, 1] + mod.coef[, 2] * tcrit
      Estimate <- mod.coef[, 1]
      OR <- exp(Estimate)
      est_ci <- paste0(round(Estimate, digits), " (", round(ci.low, digits),
                       ", ", round(ci.up, digits), ")")
      or_ci <- paste0(round(OR, digits), " (", round(or.low, digits), ", ",
                      round(or.up, digits), ")")
      if (family(model)$family == "binomial" & type == "cond") {
        cis <- data.frame(or_ci)
        out.df <- data.frame(cis, round_pval(mod.coef[, 4]))
        colnames(out.df) <- c("Odds ratio", "Pr(>|z|)")
      }
      if (family(model)$family == "binomial" & type == "ext") {
        cis <- round(data.frame(OR, mod.coef[, 2], or.low, or.up), digits)
        out.df <- data.frame(cis, round_pval(mod.coef[, 4]))
        colnames(out.df) <- c("Odds ratio", "Std. Error", "Lower CI", "Upper CI",
                              "Pr(>|z|)")
      }
      if ((family(model)$family == "poisson" |
           class(model)[1] == "negbin" |
           family(model)$family == "quasipoisson") & type == "cond") {
        cis <- data.frame(or_ci)
        out.df <- data.frame(cis, round_pval(mod.coef[, 4]))
        colnames(out.df) <- c("Rate ratio", "Pr(>|z|)")
      }
      if ((family(model)$family == "poisson" |
           class(model)[1] == "negbin" |
           family(model)$family == "quasipoisson") & type == "ext") {
        cis <- round(data.frame(OR, mod.coef[, 2], or.low, or.up), digits)
        out.df <- data.frame(cis, round_pval(mod.coef[, 4]))
        colnames(out.df) <- c("Rate ratio", "Std. Error", "Lower CI", "Upper CI",
                              "Pr(>|z|)")
      }
      if (family(model)$family == "gaussian" & type == "cond") {
        cis <- data.frame(est_ci)
        out.df <- data.frame(cis, round_pval(mod.coef[, 4]))
        colnames(out.df) <- c("Coefficient", "Pr(>|t|)")
      }
      if (family(model)$family == "gaussian" & type == "ext") {
        cis <- round(data.frame(Estimate, mod.coef[, 2], ci.low, ci.up),
                     digits)
        out.df <- data.frame(cis, round_pval(mod.coef[, 4]))
        colnames(out.df) <- c("Coefficient", "Std. Error", "Lower CI", "Upper CI",
                              "Pr(>|t|)")
      }
      if (family(model)$family == "gaussian" & type == "cond" & exp_norm == TRUE) {
        cis <- data.frame(or_ci)
        out.df <- data.frame(cis, round_pval(mod.coef[, 4]))
        colnames(out.df) <- c("Ratio", "Pr(>|t|)")
      }
      if (family(model)$family == "gaussian" & type == "ext" & exp_norm == TRUE) {
        cis <- round(data.frame(OR, mod.coef[, 2], or.low, or.up),
                     digits)
        out.df <- data.frame(cis, round_pval(mod.coef[, 4]))
        colnames(out.df) <- c("Ratio", "Std. Error", "Lower CI", "Upper CI",
                              "Pr(>|t|)")
      }
    }
    if (is.null(labels)) {
      rownames(out.df) <- rownames(out.df)
    } else {
      rownames(out.df) <- labels
    }
    out.df
  }
}
