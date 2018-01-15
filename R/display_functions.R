# Functions to display coefficients

#' Rounding p-values.
#'
#' \code{round_pval} is an internal function called by \code{glm_coef} to round p-values from model coefficients.
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @param pval vector of p-values, numeric.
round_pval <- function(pval)
{
	res <- ifelse(pval < 0.001, "< 0.001", round(pval, 3))
}

#' Table of coefficients from generalised linear models.
#'
#' \code{glm_coef} displays estimates with confidence intervals and p-values from generalised linear models (see Details).
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @details \code{glm_coef} recognises objects (models) from the following classes: \code{clm}, \code{clogit},
#' \code{coxph}, \code{gee}, \code{glm}, \code{glmerMod}, \code{lm}, \code{lme}, \code{multinom}, \code{negbin},
#' \code{polr} and \code{surveg}
#' @details For models from logistic regression (including conditional logistic, ordinal and multinomial),
#' Poisson or survival analysis, coefficient estimates and corresponding confidence intervals are
#' automatically exponentiated (back-transformed).
#' @details By default, \code{glm_coef} uses robust standard errors for calculating confidence intervals.
#' @details Please read the Vignette on Regression for more details.
#' @param model A model from any of the classes listed in the details section.
#' @param digits A scalar, number of digits for rounding the results (default = 2).
#' @param alpha Significant level (default = 0.05) used to calculate confidence intervals.
#' @param labels An optional character vector with the names of the coefficients (including intercept).
#' @param se.rob Logical, should robust errors be used to calculate confidence intervals? (default = TRUE).
#' @return A data frame with estimates, standard errors, confidence intervals and p-values from \code{glm} objects.
#' @examples
#' ## Continuous outcome.
#' data(birthwt, package = "MASS")
#' birthwt$smoke <- factor(birthwt$smoke, labels=c("Non-smoker", "Smoker"))
#' birthwt$race <- factor(birthwt$race > 1, labels=c("White", "Non-white"))
#' model_norm <- glm(bwt ~ smoke + race, data = birthwt)
#' glm_coef(model_norm)
#' glm_coef(model_norm, labels=c("Constant", "Smoker vs Non-smoker", "Non-white vs White"))
#'
#' ## Logistic regression.
#' data(diet, package = "Epi")
#' model_binom <- glm(chd ~ fibre, data = diet, family = binomial)
#' glm_coef(model_binom, labels = c("Constant", "Fibre intake (g/day)"))
#'
#' ## Poisson regression.
#' library(MASS)
#' data(quine)
#' levels(quine$Eth) <- list(White = "N", Aboriginal = "A")
#' levels(quine$Sex) <- list(Male = "M", Female = "F")
#' model_pois <- glm(Days ~ Eth + Sex + Age, family = poisson, data = quine)
#' glm_coef(model_pois)
#' deviance(model_pois) / df.residual(model_pois) # to check for overdispersion
#'
#' model_negbin <- glm.nb(Days ~ Eth + Sex + Age, data = quine)
#' unadj <- glm_coef(model_negbin, labels=c("Constant",
#'                                    "Race: Aboriginal/White",
#'                                    "Sex: Female/Male",
#'                                    "F1/Primary",
#'                                    "F2/Primary",
#'                                    "F3/Primary"))
#' unadj # Not-adjusted for multiple comparisons
#'
#' library(multcomp)
#' model_glht <- glht(model_negbin, linfct  = mcp(Age = "Tukey"))
#' age_glht <- xymultiple(model_glht, Exp = TRUE, plot = FALSE)
#' age_glht
#'
#' final <- unadj
#' final[, 5] <- as.character(final[, 5])
#' age_glht[, 5] <- as.character(age_glht[, 5])
#' final[4:6, 3:5] <- age_glht[1:3, 3:5]
#' final # Final table, with CIs and p-values adjusted for multiple comparisons (Westfall).
#'
#' ## For more examples, please read the Vignette on Regression.
glm_coef <- function(model, digits = 2, alpha = 0.05, labels = NULL, se.rob = TRUE)
{
  mod <- summary(model)
  n <- nrow(mod)
  conf.int <- 1 - alpha
  zcrit <- qnorm((1 + conf.int)/2)
  if (class(model)[1] == "clm" | class(model)[1] == "polr") {
    # Confidence intervals:
    int.low <- exp(mod$coefficients[, 1] - mod$coefficients[, 2] * zcrit)
    int.up <- exp(mod$coefficients[, 1] + mod$coefficients[, 2] * zcrit)
    # Coefficients:
    estim <- mod$coefficients[, 1]
    coeff <- exp(estim)
    coeff.error <- mod$coefficients[, 2]
    coeff.df <- round(data.frame(coeff, int.low, int.up, coeff.error), digits)
    zval <- mod$coefficients[, 1] / mod$coefficients[, 2]
    pval <- round_pval((1 - pnorm(abs(zval))) * 2)
    coeff.df <- data.frame(coeff.df, pval)
    names(coeff.df) <- c("Ordinal OR", "Lower CI", "Upper CI", "Std. Error", "Pr(>|Z|)")
    if (is.null(labels)) {
      rownames(coeff.df) <- rownames(coeff.df)
    } else {
      rownames(coeff.df) <- labels
    }
    coeff.df
  } else
    if (class(model)[1] == "multinom") {
    n.mod <- length(rownames(mod$coefficients))
    for(i in 1:n.mod) {
      est <- mod$coefficients[i, ]
      back <- exp(mod$coefficients[i, ])
      back[1] <- NA
      modse <- mod$standard.errors[i, ]
      low <- exp(mod$coefficients[i, ] - modse * zcrit)
      low[1] <- mod$coefficients[i, 1] - modse[1] * zcrit
      up <- exp(mod$coefficients[i, ] + modse * zcrit)
      up[1] <- mod$coefficients[i, 1] + modse[1] * zcrit
      zval <- est / modse
      pval <- round_pval((1 - pnorm(abs(zval))) * 2)
      df <- round(cbind(back, low, up, zval), digits)
      df <- data.frame(df, pval)
      colnames(df) <- c("Multinomial OR", paste("lower", 100 - 100 * alpha, "ci", sep = ""), paste("upper", 100 - 100 * alpha, "ci", sep = ""), "z value", "Pr(>|z|)")
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
  } else
    if (class(model)[1] == "lme") {
      mod <- summary(model)$tTable
      estim <- mod[, 1]
      low <- estim - mod[, 2] * zcrit
      up <- estim - mod[, 2] * zcrit
      stats <- mod[, 2:4]
      out.df <- round(data.frame(estim, low, up, stats), digits)
      out.df <- data.frame(out.df, round_pval(mod[, 5]))
      names(out.df) <- c("Coeff", "Lower CI", "Upper CI", "SE", "DF", "t value", "Pr(>|t|)")
      if (is.null(labels)) {
        rownames(out.df) <- rownames(out.df)
      } else {
        rownames(out.df) <- labels
      }
      out.df
    } else
      if (class(model)[1] == "clogit") {
        mod.coef <- mod$coefficients
        or.low <- exp(mod.coef[, 1] - mod.coef[, 3] * zcrit)
        or.up <- exp(mod.coef[, 1] + mod.coef[, 3] * zcrit)
        OR <- mod.coef[, 2]
        Estimate <- mod.coef[, 1]
        cis <- data.frame(OR, mod.coef[, 2], or.low, or.up)
        cis <- round(cis, digits)
        out.df <- data.frame(cis, round_pval(mod.coef[, 5]))
        colnames(out.df) <- c("OR", "Std. Error", "Lower CI", "Upper CI", "Pr(>|z|)")
        if (is.null(labels)) {
          rownames(out.df) <- rownames(out.df)
        } else {
          rownames(out.df) <- labels
        }
        out.df
      } else
        if (class(model)[1] == "survreg") {
          mod.coef <- mod$table
          hc <- sandwich::vcovHAC(model)
          ct <- lmtest::coeftest(model, vcov=hc)
          np <- nrow(mod.coef)
          if (se.rob == TRUE) {
            mod.coef <- ct[, ]
          } else {
            mod.coef <- mod$table
          }
          Estimate <- mod.coef[, 1]
          OR <- exp(Estimate)
          or.low <- exp(Estimate - mod.coef[, 2] * zcrit)
          or.up <- exp(Estimate + mod.coef[, 2] * zcrit)
          cis <- data.frame(OR, mod.coef[, 2], or.low, or.up)
          cis <- round(cis, digits)
          out.df <- data.frame(cis, round_pval(mod.coef[, 4]))
          out.df <- out.df[-1, ]
          colnames(out.df) <- c("Survival time ratio", "Std. Error", "Lower CI", "Upper CI", "Pr(>|z|)")
          if (is.null(labels)) {
            labs <- rownames(out.df)
            rownames(out.df) <- gsub("Log(scale)", replacement = "Scale", labs, fixed = TRUE)
          } else {
            rownames(out.df) <- labels
          }
          out.df
        } else
          if (class(model)[1] == "coxph") {
            mod.coef <- mod$coefficients
            or.low <- exp(mod.coef[, 1] - mod.coef[, 3] * zcrit)
            or.up <- exp(mod.coef[, 1] + mod.coef[, 3] * zcrit)
            OR <- mod.coef[, 2]
            Estimate <- mod.coef[, 1]
            cis <- data.frame(OR, mod.coef[, 3], or.low, or.up)
            cis <- round(cis, digits)
            out.df <- data.frame(cis, round_pval(mod.coef[, 5]))
            colnames(out.df) <- c("Hazard ratio", "Std. Error", "Lower CI", "Upper CI", "Pr(>|z|)")
            if (is.null(labels)) {
              rownames(out.df) <- rownames(out.df)
            } else {
              rownames(out.df) <- labels
            }
            out.df
          } else
            if (class(model)[1] == "gee" & (family(model)$family == "binomial" |
                                            family(model)$family == "poisson") |
                                            family(model)$family == "quasi") {
              mod <- mod$coefficients
              n <- nrow(mod)
              estim <- mod[, 1]
              coeff <- exp(mod[, 1])
              coeff[1] <- NA
              if (se.rob == TRUE) {
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
              out.df <- round(data.frame(estim, coeff, low, up, segee), digits)
              out.df <- data.frame(out.df, round_pval(pvals))
              names(out.df) <- c("Coeff", "Exp(Coeff)", "Lower CI", "Upper CI", "SE", "Pr(>|z|)")
              if (is.null(labels)) {
                rownames(out.df) <- rownames(out.df)
              } else {
                rownames(out.df) <- labels
              }
              out.df
              } else
                if (class(model)[1] == "gee" & family(model)$family == "gaussian") {
                  mod <- mod$coefficients
                  n <- nrow(mod)
                  estim <- mod[, 1]
                  if (se.rob == TRUE) {
                    segee <- mod[, 4]
                    my.z <- mod[, 5]
                  } else {
                    segee <- mod[, 2]
                    my.z <- mod[, 3]
                  }
                  low <- estim - segee * zcrit
                  up <- estim + segee * zcrit
                  pvals <- (1 - pnorm(abs(my.z))) * 2
                  out.df <- round(data.frame(estim, low, up, segee), digits)
                  out.df <- data.frame(out.df, round_pval(pvals))
                  names(out.df) <- c("Coeff", "Lower CI", "Upper CI", "SE", "Pr(>|z|)")
                  if (is.null(labels)) {
                    rownames(out.df) <- rownames(out.df)
                  } else {
                    rownames(out.df) <- labels
                  }
                  out.df
                  } else
                    if (class(model)[1] == "glmerMod" & (family(model)$family == "binomial" | family(model)$family == "poisson")) {
                      mod <- mod$coefficients
                      n <- nrow(mod)
                      estim <- mod[, 1]
                      coeff <- exp(mod[, 1])
                      coeff[1] <- NA
                      low <- estim - mod[, 2] * zcrit
                      low[2:n] <- exp(low[2:n])
                      up <- estim + mod[, 2] * zcrit
                      up[2:n] <- exp(up[2:n])
                      stats <- mod[, 2:3]
                      out.df <- round(data.frame(estim, coeff, low, up, stats), digits)
                      out.df <- data.frame(out.df, round_pval(mod[, 4]))
                      names(out.df) <- c("Coeff", "Exp(Coeff)", "Lower CI", "Upper CI", "SE", "z value", "Pr(>|z|)")
                      if (is.null(labels)) {
                        rownames(out.df) <- rownames(out.df)
                      } else {
                        rownames(out.df) <- labels
                      }
                      out.df
                      } else
                        if (class(model)[1] == "glmerMod" & family(model)$family == "gaussian") {
                          mod <- mod$coefficients
                          n <- nrow(mod)
                          estim <- mod[, 1]
                          low <- estim - mod[, 2] * zcrit
                          up <- estim + mod[, 2] * zcrit
                          stats <- mod[, 2:3]
                          out.df <- round(data.frame(estim, low, up, stats), digits)
                          out.df <- data.frame(out.df, round_pval(mod[, 4]))
                          names(out.df) <- c("Coeff", "Lower CI", "Upper CI", "SE", "z value", "Pr(>|z|)")
                          if (is.null(labels)) {
                            rownames(out.df) <- rownames(out.df)
                          } else {
                            rownames(out.df) <- labels
                          }
                          out.df
                          } else  {
                            if (class(model)[1] != "glm" & class(model)[1] != "negbin" & class(model)[1] != "lm")
                              stop("Object's class or family error not supported.")
                            tcrit <- qt((1 + conf.int)/2, df.residual(model))
                            hc <- sandwich::vcovHAC(model)
                            ct <- lmtest::coeftest(model, vcov=hc)
                            if (se.rob == TRUE) {
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
                            if (family(model)$family=="binomial" |
        family(model)$family=="poisson" |
        substr(family(model)$family,1,12)=="Negative Bin" |
        family(model)$family=="quasipoisson") {
                              cis <- data.frame(OR, mod.coef[,2], or.low, or.up)
                              } else {
                                cis <- data.frame(Estimate, mod.coef[,2], ci.low, ci.up)
                                }
                            cis <- round(cis, digits)
                            out.df <- data.frame(cis, round_pval(mod.coef[,4]))
                            if (family(model)$family=="binomial"){
                              colnames(out.df) <- c("OR", "Std. Error", "Lower CI", "Upper CI", "Pr(>|z|)")
                              } else {
                                if (family(model)$family=="poisson" |
          substr(family(model)$family,1,12)=="Negative Bin" |
          family(model)$family=="quasipoisson") {
                                  colnames(out.df) <- c("IRR", "Std. Error", "Lower CI", "Upper CI", "Pr(>|z|)")
                                  } else  {
                                    colnames(out.df) <- c("Estimate", "Std. Error", "Lower CI", "Upper CI", "Pr(>|t|)")
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


