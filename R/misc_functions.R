# Miscellaneous functions

#' Sum of squares for Jackknife.
#'
#' \code{ss_jk} is an internal function called by \link{jack_knife}. It calculates the squared
#' difference of a numerical variable around a given value (for example, the mean).
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @param obs A numerical vector with no missing values (NA's).
#' @param stat The value of the statistic that is used as a reference.
#' @return The squared difference between a variable and a given value.
#' @examples
#' x <- rnorm(10, 170, 8)
#' x
#' mean(x)
#' ss_jk(x, mean(x))
#' jack_knife(x)
ss_jk <- function(obs, stat)
{
	n <- length(obs)
	ss <- numeric(n)
	for(i in 1:n) ss[i] <- ((obs[i] - stat)^2)
	ss
}

#' Ranks leverage observations from Jackknife method.
#'
#' \code{jack_knife} Ranks the squared differences between mean values from Jackknife analysis
#' (arithmetic mean estimated by removing one observation at a time) and the original mean value.
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @param x A numeric variable. Missing values are removed by default.
#' @return Data frame with the ranked squared differences.
#' @seealso \link{rank_leverage}
#' @examples
#' x <- rnorm(10, 170, 8)
#' x
#' mean(x)
#' jack_knife(x)
#'
#' x <- rnorm(100, 170, 8)
#' mean(x)
#' head(jack_knife(x))
jack_knife <- function(x)
{
	jk.means <- knife_mean(x)
	mu <- mean(x, na.rm = TRUE)
	jk.ss <- ss_jk(jk.means, mu)
	jk.ord <- order(jk.ss, decreasing = TRUE)
	Observation <- jk.ord
	SS <- jk.ss[jk.ord]
	Mean <- jk.means[jk.ord]
	res <- data.frame(Observation, Value = x[jk.ord], Mean, SS)
	rownames(res) <- NULL
	res
}

#' Leverage.
#'
#' \code{leverage} is an internal function called by \link{rank_leverage}.
#'
#' Estimates the leverage of each observation around the arithmetic mean.
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @param x A numeric variable. Missing values are removed by default.
#' @return Variable with corresponding leverage estimations
#' @examples
#' x <- rnorm(10, 170, 8)
#' x
#' mean(x)
#' leverage(x)
#' rank_leverage(x)
leverage <- function(x)
{
	inv <- 1 / length(x)
	ssx <- (x - mean(x, na.rm = TRUE))^2
	ssx2 <- sum(ssx)
	lev <- inv + ssx / ssx2
	lev
}

#' Ranks observations by leverage.
#'
#' \code{rank_leverage} ranks observations by their leverage (influence) on the arithmetic mean.
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @param x A numeric variable. Missing values are removed by default.
#' @return A data frame ranking observations by their leverage around the mean.
#' @seealso \link{jack_knife}
#' @examples
#' x <- rnorm(10, 170, 8)
#' x
#' mean(x)
#' rank_leverage(x)
#'
#' x <- rnorm(100, 170, 8)
#' mean(x)
#' head(rank_leverage(x))
rank_leverage <- function(x)
{
	lev <- leverage(x)
	lev.ord <- order(lev, decreasing = TRUE)
	Observation <- lev.ord
	Leverage <- lev[lev.ord]
	res <- data.frame(Observation, Leverage)
	res
}


#' Given y solve for x in a simple linear model.
#'
#' \code{predict_inv} Calculates the value the predictor x that generates value y with a simple linear model.
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @param model A simple linear model object (class lm).
#' @param y A numerical scalar, the value of the outcome for which we want to calculate the predictor x.
#' @return The estimated value of the predictor.
#' @examples
#' ## Spectrophotometry example. Titration curve for riboflavin (nmol/ml). The sample has an absorbance
#' ## of 1.15. Aim is to estimate the concentration of riboflavin in the sample.
#'
#' Riboflavin <- seq(0, 80, 10)
#' OD <- 0.0125 * Riboflavin + rnorm(9, 0.6, 0.03)
#' titration <- data.frame(Riboflavin, OD)
#'
#' xyplot(OD ~ Riboflavin, data = titration, pch = 16, col = 1, aspect = 3/4) +
#'   layer(panel.smoother(lwd = 1.5, col = 2, method = "lm", ...))
#'
#' ## Model with intercept different from zero:
#' model <- lm(OD ~ Riboflavin, data = titration)
#' glm_coef(model)
#' predict_inv(model, 1.15)
predict_inv <- function(model, y)
{
  param <- length(coef(model))
  m <- ifelse(param == 1, as.numeric(coef(model)), coef(model)[2])
  b <- ifelse(param == 1, 0, coef(model)[1])
	x <- (y - b) / m
	names(x) <- names(m)
	x
}


#' Jackknife for means.
#'
#' \code{knife_mean} is an internal function. Calculates arithmetic means by removing one observation
#' at a time.
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @param x A numerical variable. Missing values are removed for the mean calculation.
#' @return A vector with the mean calculations.
#' @examples
#' x <- rnorm(10, 170, 8)
#' x
#' mean(x)
#' knife_mean(x)
knife_mean <- function(x)
{
	n <- length(x)
	jk <- numeric(n)
	for(i in 1:n) jk[i] <- mean(x[-i], na.rm = TRUE)
	jk
}

#' Relative and Cumulative Frequency.
#'
#' \code{freq_cont} tabulates a continuous variable by given classes.
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @param x A numerical (continuous) variable. Ideally, relatively long (greater than 100 observations).
#' @param bks Breaks defining the classes (see example).
#' @param dg Number of digits for rounding (default = 2).
#' @return A data frame with the classes, the mid-point, the frequencies, the relative and cumulative frequencies.
#' @seealso \link{hist}
#' @examples
#' data(IgM, package="ISwR")
#' Ab <- data.frame(IgM)
#' estat(~ IgM, data = Ab)
#' freq_cont(IgM, seq(0, 4.5, 0.5))
freq_cont <- function(x, bks, dg = 2)
{
	cx <- cut(x, breaks = bks)
	n <- length(x)
	x.tab <- as.data.frame(table(cx))
	rel <- x.tab$Freq / n
	n2 <- length(rel)
	mids <- numeric(n2)
	cum <- cumsum(rel)
	for(i in 1:n2) mids[i] <- mean(c(bks[i], bks[i + 1]))
	rel <- round(rel, digits = dg)
	cum <- round(cum, digits = dg)
	res <- data.frame(x.tab[, 1], mids, x.tab[, 2], rel, cum)
	names(res) <- c('Class', 'Mids', 'Counts', 'rel.freq', 'cum.freq')
	res
}

#' Goodness of fit for Logistic Regression.
#'
#' \code{logistic_gof} performs the Hosmer and Lemeshow test to test the goodness of fit of a logistic
#' regression model. This function is part of \code{residuals.lrm} from package \code{rms}.
#'
#' @author Frank Harell, Vanderbilt University <f.harrell@vanderbilt.edu>
#' @param model A logistic regression model object.
#' @references Hosmer DW, Hosmer T, Lemeshow S, le Cessie S, Lemeshow S. A comparison of goodness-of-fit
#' tests for the logistic regression model. Stat in Med 16:965–980, 1997.
#' @examples
#' data(diet, package = "Epi")
#' model <- glm(chd ~ fibre, data = diet, family = binomial)
#' glm_coef(model, labels = c("Constant", "Fibre intake (g/day)"))
#' logistic_gof(model)
logistic_gof <- function(model)
{
	L <- model$linear.predictors
	if (length(L) == 0)
      stop("you did not use linear.predictors=TRUE for the fit")
	cof <- model$coefficients
	P <- 1 / (1 + exp(-L))
	Y <- model$y
	rnam <- names(Y)
	cnam <- names(cof)
	if (!is.factor(Y))
      Y <- factor(Y)
  lev <- levels(Y)
  Y <- unclass(Y) - 1
	X <- glm(model$formula, family = binomial, x = TRUE, data = model$data)$x
  y <- Y >= 1
  p <- P
  sse <- sum((y - p)^2)
  wt <- p * (1 - p)
  d <- 1 - 2 * p
  z <- lm.wfit(X, d, wt)
  res <- z$residuals * sqrt(z$weights)
  sd <- sqrt(sum(res^2))
  ev <- sum(wt)
  z <- (sse - ev) / sd
  P <- 2 * (1 - pnorm(abs(z)))
  stats <- data.frame(sse, ev, sd, z, P)
	names(stats) <- c("Sum of squared errors", "Expected value|H0",
  "SD", "Z", "P")
  return(drop(stats))
}



#' Coefficient of determination.
#'
#' \code{coef_det} estimates the coefficient of determination (r-squared) from fitted (predicted) and
#' observed values. Outcome from the model is assumed to be numerical.
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @param obs Vector with observed values (numerical outcome).
#' @param fit Vector with fitted (predicted) values.
#' @return A scalar, the coefficient of determination (r-squared).
#' @examples
#' ## Linear regression:
#' Riboflavin <- seq(0, 80, 10)
#' OD <- 0.0125*Riboflavin + rnorm(9, 0.6, 0.03)
#' titration <- data.frame(Riboflavin, OD)
#' model1 <- lm(OD ~ Riboflavin, data=titration)
#' summary(model1)
#' coef_det(titration$OD, fitted(model1))
#'
#' ## Non-linear regression:
#' library(nlme)
#' data(Puromycin)
#' mm.tx <- gnls(rate ~ SSmicmen(conc, Vm, K), data=Puromycin, subset=state=="treated")
#' summary(mm.tx)
#' coef_det(Puromycin$rate[1:12], mm.tx$fitted)
coef_det <- function(obs, fit)
{
	obm <- mean(obs, na.rm=TRUE)
	SSres <- sum((obs-fit)^2, na.rm = TRUE)
	SStot <- sum((obs-obm)^2, na.rm = TRUE)
	res <- 1 - (SSres/SStot)
	res
}

#' Ranks observations based upon influence measures on models.
#'
#' \code{rank_influence} calculates influence measures of each data observation on models and then ranks them.
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @author Several R core team members and John Fox, originally in his \code{car} package.
#' @param model A generalised linear model object.
#' @details \code{rank_influence} is a wrap function that calls \link{influence.measures}, ranks observations on
#' their significance influence on the model and displays the 10 most influential observations
#' (if they are significant).
#' @examples
#' data(diet, package = "Epi")
#' model <- glm(chd ~ fibre, data = diet, family = binomial)
#' rank_influence(model)
rank_influence <- function(model)
{
	 infl <- influence.measures(model)
	 mat <- (infl$is.inf)*1
	 signifa <- sort(apply(mat, 1, sum), decreasing=TRUE)
	 mat.sum <- signifa[signifa>0]
	 mat.names <- as.numeric(names(mat.sum))
	 pos <- match(mat.names, row.names(infl$infmat), nomatch = 0)
	 if (length(pos) > 10) {
		 pos2 <- pos[1:10]
		 signifb <- as.vector(mat.sum)[1:10]
		 res <- data.frame(infl$infmat[pos2,], signif=signifb)
	 }
	 if (length(pos) < 10 & length(pos) > 0) {
		 signifb <- as.vector(mat.sum)
		 res <- data.frame(infl$infmat[pos,], signif=signifb)
	 }
	 if (length(pos)==0) {
		 res <- 'No significant influence measures'
	 }
	 res
}

#' Inverse of the logit
#'
#' \code{inv_logit} Calculates the inverse of the logit (probability in logistic regression)
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @param x Numerical value used to compute the inverse of the logit.
inv_logit <- function(x)
{
  if (any(omit <- is.na(x))) {
    lv <- x
    lv[omit] <- NA
    if (any(!omit))
      lv[!omit] <- Recall(x[!omit])
    return(lv)}
  exp(x)/(1 + exp(x))
}

#' Pseudo R2 (logistic regression)
#' \code{pseudo_r2} Calculates R2 analogues (pseudo R2) of logistic regression.
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @param model A logistic regression model.
#' @details \code{pseudo_r2} calculates three pseudo R2 of logistic regression models: 1) Nagelkerke, @0 Cox and Snell, 3) Hosmer and Lemeshow.
#' @return A data frame with the calculated pseudo R2 values.
#' @examples
#' model_oncho <- glm(mf ~ area, data = Oncho, binomial)
#' glm_coef(model_oncho, labels = c("Constant", "Area (rainforest/savannah)"))
#' pseudo_r2(model_oncho)
pseudo_r2 <- function(model)
{
  if (!(model$family$family == "binomial" && (model$family$link ==
                                              "logit" || model$family$link == "probit"))) {
    stop("No logistic regression model, no pseudo R^2 computed.")
  }
  dev <- model$deviance
  null_dev <- model$null.deviance
  model_n <- length(model$fitted)
  r_hl <- 1 - dev / null_dev
  r_cs <- 1 - exp(-(null_dev - dev) / model_n)
  r_n <- r_cs / (1 - (exp(-(null_dev / model_n))))
  Index <- c("Nagelkerke", "Hosmer and Lemeshow", "Cox and Snell")
  Estimate <- round(c(r_n, r_hl, r_cs), 3)
  res <- data.frame(Index, Estimate)
  res
}
