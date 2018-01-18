# Functions for Plotting

#' Multiple comparisons with plot.
#'
#' \code{xymultiple} displays results from post-doc analysis and constructs corresponding plot.
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @author Deepayan Sarkar, R-Core.
#' @param model.glht An object of class "glht" (from multiple comparisons).
#' @param method Method passed to \link{summary.glht} (see details).
#' @param Exp Logical, should results be exponentiated? (default = FALSE).
#' @param dg Number of digits for rounding (default = 2).
#' @param plot Logical, should a plot be constructed? (default = TRUE).
#' @param ... Passes additional information to \code{xyplot}.
#' @details
#' The default adjusting method is "Westfall". Other options are: "single-step", "Shaffer",
#' "free", "holm", "hochberg", "hommel", "bonferroni".
#' @seealso \code{glht}, \code{glht-methods}.
#' @return A data frame with CIs and p-values adjusted for multiple comparisons.
#' @examples
#' library(multcomp)
#' data(birthwt)
#' birthwt$race <- factor(birthwt$race, labels = c("White", "Black", "Other"))
#' model1 <- aov(bwt ~ race, data = birthwt)
#' model1_glht <- glht(model1, linfct = mcp(race = "Tukey"))
#' xymultiple(model1_glht)
#'
#' model2 <- glm(low ~ race, data = birthwt, family = binomial)
#' model2_glht <- glht(model2, linfct = mcp(race="Tukey"))
#' xymultiple(model2_glht, Exp = TRUE)
xymultiple <- function(model.glht, method = "Westfall", Exp = FALSE, dg = 2, plot = TRUE, ...)
{
  model.summary <- summary(model.glht, adjusted(type = method))
  cisnorm <- as.data.frame(confint(model.summary, adjusted(type=method))$confint)
  if (Exp == TRUE) {
    cisnorm <- exp(cisnorm)
    names(cisnorm)[1] <- "Ratio"
  } else {
    cisnorm <- cisnorm
    names(cisnorm)[1] <- "Difference"
  }
  cis <- round(cisnorm, dg)
  pval <- round_pval(as.vector(model.summary$test$pvalues))
  model.ci <- data.frame(Comparison = row.names(cis), cis, pval)
  names(model.ci)[5] <- "Pr(>|Z|)"
  rownames(model.ci) <- NULL
  if (plot == FALSE) {
    model.ci
  } else {
    if (Exp == FALSE) {
      fig <- xyplot(cbind(Difference, lwr, upr) ~ Comparison, data = model.ci, ylab = "Difference",
                    xlab = NULL, panel = panel.errbars, pch = 20, col = 1, scales = list(x = list(rot = 45)), ...) +
        layer(panel.abline(h = 0, lty = 2, col = 2, lwd = 0.5))
      print(model.ci)
      fig
    } else {
      fig <- xyplot(cbind(Ratio, lwr, upr) ~ Comparison, data=model.ci, ylab = "Ratio",
                    xlab = NULL, panel = panel.errbars, pch = 20, col = 1, scales = list(x = list(rot = 45)), ...) +
        layer(panel.abline(h = 1, lty = 2, col = 2, lwd = 0.5))
      print(model.ci)
      fig
    }
  }
}

#' Bland-Altman agreement plots.
#'
#' @details \code{bland_altman} constructs Bland-Altman agreement plots.
#' @details Variables in \code{formula} are continuous paired observations. When the distribution of the outcome
#' is not normal, but becomes normal with a log-transformation, \code{bland_altman} can plot the ratio between
#' outcomes (difference in the log scale) by using option \code{transform = TRUE}.
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @author Deepayan Sarkar, R-Core.
#' @param formula A formula of the form y ~ x (see details).
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param transform Logical, should ratios instead of difference be used to construct the plot?
#' @param aspect Physical aspect ratio passed to \link{xyplot}.
#' @param ... Further arguments passed to \link{xyplot}.
#' @examples
#' data(wright, package = "ISwR")
#' bland_altman(mini.wright ~ std.wright, data = wright, pch = 16)
#' bland_altman(mini.wright ~ std.wright, data = wright, pch = 16,
#'              ylab = "Large-mini expiratory flow rate (l/min)",
#'              xlab = "Mean expiratory flow rate (l/min)")
bland_altman <- function(formula, data, transform = FALSE, aspect = 3/4, ...)
{
  y <- NULL; rm(y)
  vars <- all.vars(formula)
  resp <- vars[1]
  expl <- vars[2]
  diff <- data[[expl]] - data[[resp]]
  ratio <- data[[expl]] / data[[resp]]
  avg <- apply(cbind(data[[expl]], data[[resp]]), 1, mean, na.rm = TRUE)
  df <- data.frame(Difference = diff, Ratio = ratio, Average = avg)
  if (transform == TRUE) {
    xyplot(Ratio ~ Average, data = df, aspect=aspect, ...) +
      layer(panel.abline(h = mean(y, na.rm = TRUE), col = 2, lwd = 1.5)) +
      layer(panel.abline(h = as.numeric(reference_range(mean(y, na.rm = TRUE), sd(y, na.rm = TRUE)))[1], col = 2, lwd = 1.5, lty = 2)) +
      layer(panel.abline(h = as.numeric(reference_range(mean(y, na.rm = TRUE), sd(y, na.rm = TRUE)))[2], col = 2, lwd = 1.5, lty = 2))
  } else {
    xyplot(Difference ~ Average, data = df, ...) +
      layer(panel.abline(h = mean(y, na.rm = TRUE), col = 2, lwd = 1.5)) +
      layer(panel.abline(h = as.numeric(reference_range(mean(y, na.rm = TRUE), sd(y, na.rm = TRUE)))[1], col = 2, lwd = 1.5, lty = 2)) +
      layer(panel.abline(h = as.numeric(reference_range(mean(y, na.rm = TRUE), sd(y, na.rm = TRUE)))[2], col = 2, lwd = 1.5, lty = 2))
  }
}

#' Residual vs Fitted plot.
#'
#' \code{rvf_plot} plots studentized residuals against fitted values from \link{glm} objects using package \code{lattice}.
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @author Deepayan Sarkar, R-Core.
#' @param model A \link{glm} or \link{lm} object with a numerical outcome.
#' @param pch Point character passed to \link{xyplot}.
#' @param col Colour passed to \link{xyplot}.
#' @param aspect Physical aspect ratio passed to \link{xyplot}.
#' @param ... Further arguments passed to \link{xyplot}.
#' @examples
#' data(thuesen, package = "ISwR")
#' model <- lm(short.velocity ~ blood.glucose, data = thuesen)
#' plot(model, which = 1)
#' rvf_plot(model)
rvf_plot <- function(model, pch=20, col=1, aspect=3/4, ...)
{
  res = MASS::studres(model); fit = fitted(model)
  rvf <- data.frame(res, fit)
  xyplot(res ~ fit, data = rvf, xlab = "Fitted values", pch = pch, col = col, aspect = aspect,
         ylab="Studentized residuals", ...) +
    layer(panel.abline(h = 0, lty = 2, col = "gray", lwd = 0.8)) +
    layer(panel.abline(h = 2, lty = 2, col = "gray", lwd = 0.8)) +
    layer(panel.abline(h = -2, lty = 2, col = "gray", lwd = 0.8)) +
    layer(panel.smoother(lwd = 1.2, col = 2, se = FALSE, ...))
}

#' Construct "pretty" box plots in lattice.
#'
#' \code{box_plot} is a wrap function that calls \link{bwplot} to construct more aesthetic box plots.
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @author Deepayan Sarkar, R-Core.
#' @param formula A formula of the form \code{y ~ x} where \code{y} is a numerical variable and \code{x} is a factor.
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param box.fill Colour used for the box passed to \link{bwplot}.
#' @param box.ratio Ratio of box passed to \link{bwplot}.
#' @param aspect Physical aspect ratio passed to \link{bwplot}.
#' @param ... Further arguments passed to \link{bwplot}.
#' @examples
#' data(kfm, package = "ISwR")
#' box_plot(dl.milk ~ sex, data = kfm, ylab = "Breast-milk intake (dl/day)")
box_plot <- function(formula, data, box.fill = "gray70", box.ratio = 0.7, aspect = 3/4, ...)
{
  bwplot(formula, data, pch="|", box.ratio = box.ratio, aspect = aspect,
         par.settings = list(plot.symbol = list(pch = 20, col = 1),
                             box.rectangle = list(col = 1, lwd = 2, fill = box.fill),
                             box.umbrella = list(col = 1, lty = 1)), ...)
}

#' Quantile-quantile plots against the standard Normal distribution.
#'
#' \code{qq_plot} constructs quantile-quantile plots against the standard normal distribution
#' (also known as quantile-normal plots).
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @author Deepayan Sarkar, R-Core.
#' @param formula A formula of the form \code{ ~ x} or \code{ ~ x|z} where \code{x} is a numerical variable and \code{z} is a factor.
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param pch Point character passed to \link{qqmath}.
#' @param col Colour passed to \link{qqmath}.
#' @param aspect Physical aspect ratio passed to \link{qqmath}.
#' @param ... Further arguments passed to \link{qqmath}.
#' @examples
#' data(kfm, package = "ISwR")
#' qq_plot(~ dl.milk, data = kfm, ylab = "Breast-milk intake (dl/day)")
#' qq_plot(~ dl.milk|sex, data = kfm, ylab = "Breast-milk intake (dl/day)", aspect = 1)
qq_plot <- function(formula, data = NULL, pch = 20, col = 1, aspect = 3/4, ...)
{
  x <- NULL; rm(x)
  qqmath(formula, data = data, aspect = aspect, pch = pch, col = col, ...) +
    layer(panel.qqmathline(x, col = 2, lwd = 1.5))
}

#' Histogram with Normal density curve.
#'
#' \code{hist_norm} constructs histograms in lattice and adds corresponding Normal density curve.
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @author Deepayan Sarkar, R-Core.
#' @param formula A formula of the form \code{ ~ x} where \code{x} is a numerical variable.
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param col Colour passed to \link{histogram}.
#' @param aspect Physical aspect ratio passed to \link{histogram}.
#' @param ... Further arguments passed to \link{histogram}.
#' @examples
#' data(birthwt, package = "MASS")
#' hist_norm(~ bwt, data = birthwt, nint = 15, xlab = "Birth weight (g)")
hist_norm <- function(formula, data = NULL, col = "gray70", aspect = 3/4, ...)
{
	x <- all.vars(formula)
	m <- mean(data[[x]], na.rm = TRUE)
	s <- sd(data[[x]], na.rm = TRUE)
	histogram(formula, data=data, type='density', col=col, aspect=aspect, ...,
            panel = function(x, ...) {
            panel.histogram(x, ...)
            panel.mathdensity(dmath = dnorm, col = "black", lwd=1.5,
                              args = list(mean=m,sd=s))})
}

#' Strip plots with error bars.
#'
#' \code{strip_error} constructs strip plots in \code{lattice} with error bars showing 95% bootstrapped
#' confidence intervals around mean values.
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @author Deepayan Sarkar, R-Core.
#' @param formula A formula of the form \code{y ~ x} or \code{y ~ x|z} where \code{y} is a
#' numerical variable and both \code{x} and \code{z} are factors.
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param pch Point character passed to \link{stripplot}.
#' @param aspect Physical aspect ratio passed to \link{stripplot}.
#' @param ... Further arguments passed to \link{stripplot}.
#' @examples
#' data(energy, package="ISwR")
#' strip_error(expend ~ stature, data = energy, xlab = "Stature", ylab = "Energy expenditure (MJ)")
#'
#' ## Adding an horizontal line to show significant difference:
#' fig <- strip_error(expend~stature, data=energy, xlab="Stature",
#'                    ylab="Energy expenditure (MJ)", ylim=c(5.5,14))
#' fig + layer(panel.segments(1, 13.3, 2, 13.3, lwd=1.5)) + layer(panel.text(1.5, 13.5, "*"))
#'
#' data(birthwt, package = "MASS")
#' birthwt$smoke <- factor(birthwt$smoke, labels = c("Non-smoker", "Smoker"))
#' birthwt$Race <- factor(birthwt$race > 1, labels = c("White", "Non-white"))
#' strip_error(bwt ~ Race|smoke, data = birthwt, ylab = "Birth weight (g)")
strip_error <- function(formula, data, pch = 20, aspect = 3/4, ...)
{
  varss <- all.vars(formula)
  nvs <- length(varss)
  CI <- gen_bst_df(formula, data = data)
  names(CI)[1] <- "resp"
  names(CI)[4] <- "expl"
  if (nvs == 2) {
    stripplot(formula, data = data, jitter = TRUE, pch = pch, aspect = aspect, ...)  +
  	xyplot(cbind(resp, LowerCI, UpperCI) ~ expl, data = CI,
  		panel=panel.errbars, pch = 16, ewidth = 0.2, col = 1, make.grid = "none")
  } else {
    names(CI)[5] <- "strat"
    stripplot(formula, data = data, jitter = TRUE, pch = pch, ...)  +
  	xyplot(cbind(resp, LowerCI, UpperCI) ~ expl|strat, data = CI,
  		panel=panel.errbars, pch = 16, ewidth = 0.2, col = 1, make.grid = "none")
  }
}

#' Bar charts with error bars.
#'
#' \code{bar_error} constructs bar charts in \code{lattice} with error bars showing 95% bootstrapped
#' confidence intervals around mean values. High of bars represent mean values.
#' @details Limits for the y-axis have to be estimated; lower limit should be zero and upper limit higher than
#' the maximum upper confidence interval.
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @author Deepayan Sarkar, R-Core.
#' @param formula A formula of the form \code{y ~ x} or \code{y ~ x|z} where \code{y} is a
#' numerical variable and both \code{x} and \code{z} are factors.
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param col Colour passed to \link{barchart}.
#' @param aspect Physical aspect ratio passed to \link{barchart}.
#' @param ... Further arguments passed to \link{barchart}.
#' @examples
#' data(birthwt, package = "MASS")
#' birthwt$smoke <- factor(birthwt$smoke, labels = c("Non-smoker", "Smoker"))
#' gen_bst_df(bwt ~ smoke, data = birthwt) # To estimate limits of y-axis.
#' bar_error(bwt ~ smoke, data = birthwt, ylab = "Birth weight (g)", ylim = c(0, 3500))
#'
#' birthwt$race <- factor(birthwt$race, labels = c("White", "African American", "Other"))
#' gen_bst_df(bwt ~ race|smoke, data = birthwt) # To estimate limits of y-axis.
#' bar_error(bwt ~ race|smoke, data = birthwt, ylab = "Birth weight (g)", ylim = c(0, 3800))
#'
#' bar_error(bwt ~ race|smoke, data = birthwt, ylab = "Birth weight (g)", ylim = c(0, 3800),
#'           col = c("gray95", "gray20", "gray50"))
bar_error <- function(formula, data, col = "gray70", aspect = 3/4, ...)
{
  varss <- all.vars(formula)
  nvs <- length(varss)
  CI <- gen_bst_df(formula, data=data)
  names(CI)[1] <- "resp"
  names(CI)[4] <- "expl"
  if (nvs == 2) {
    barchart(resp~expl, data=CI, col=col, aspect=aspect, ...)  +
  	xyplot(cbind(resp,LowerCI,UpperCI) ~ expl, data=CI,
  		panel=panel.errbars, pch='-', ewidth=0.2, col=1, make.grid="none")
  } else {
    names(CI)[5] <- "strat"
    barchart(resp~expl|strat, data=CI, col=col, ...) +
    xyplot(cbind(resp,LowerCI,UpperCI) ~ expl|strat, data=CI,
      panel=panel.errbars, pch='-', ewidth=0.2, col=1, make.grid="none")
  }
}

#' Plot of model coefficients.
#'
#' \code{coef_plot} Constructs plot displaying estimates of parameters with bars representing confidence intervals.
#' @details \code{coef_plot} does not show estimate for the constant (intercept). Estimates and confidence intervals
#' can be optionally exponentiated, in which case estimates would represent ratios instead of differences.
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @author Deepayan Sarkar, R-Core.
#' @param model A regression object (like \link{glm}).
#' @param Exp Logical, should estimates and confidence intervals be exponentiated?
#' @param CI Proportion representing the confidence intervals.
#' @param pch Point character passed to \link{xyplot}.
#' @param col Colour passed to \link{xyplot}.
#' @param ... Further arguments passed to \link{xyplot}.
#' @examples
#' data(birthwt, package = "MASS")
#' birthwt$smoke <- factor(birthwt$smoke, labels=c("Non-smoker", "Smoker"))
#' birthwt$race <- factor(birthwt$race > 1, labels=c("White", "Non-white"))
#' model1 <- glm(bwt ~ smoke + race, data = birthwt)
#' glm_coef(model1, labels=c("Constant", "Smoker vs Non-smoker", "Non-white vs White"))
#' coef_plot(model1)
coef_plot <- function(model, Exp = FALSE, CI = 0.95, pch = 20, col = 1, ...)
{
	alpha <- 1 - CI
	if (Exp == FALSE) {
		cis <- glm_coef(model)
		cis <- cis[-1, c(1, 3:4)]
	 	param <- row.names(car::Anova(model))
	 	cis2 <- data.frame(Parameter = param, Estimate = cis[, 1], low = cis[, 2], high = cis[, 3])
	 	xyplot(cbind(Estimate, low, high) ~ Parameter, data = cis2, ylab = "Coefficient (difference)",
	 	panel = panel.errbars, pch = pch, col = col, ...) +
	 	layer(panel.abline(h = 0, lty = 2, col = 1, lwd = 0.5))
	} else {
		cis <- glm_coef(model)
		cis <- cis[-1, c(2, 4:5)]
	 	param <- row.names(car::Anova(model))
	 	cis2 <- data.frame(Parameter = param, Estimate = cis[, 1], low = cis[, 2], high = cis[, 3])
	 	xyplot(cbind(Estimate, low, high) ~ Parameter, data = cis2, ylab = "Coefficient (ratio)",
	 	panel = panel.errbars, pch = pch, col = col, ...) +
	 	layer(panel.abline(h = 1, lty = 2, col = 1, lwd = 0.5))}
}

#' Internal function for displaying error bars in lattice plots.
#'
#' \code{panel.errbars1} is an internal function called by \code{panel.errbars1}.
#' @param x A numeric vector with the x positions for the bars.
#' @param y0 A numeric vector with the lower limits of the bars.
#' @param y1 A numeric vector with the upper limits of the bars.
#' @param ewidth An integer.
panel.errbars1 <- function(x, y0, y1, ewidth = 0)
{
  x <- as.numeric(x)
  offs <- ewidth / 2
  panel.segments(x0 = x, x1 = x, y0 = y0, y1 = y1)
  panel.segments(x0 = x - offs, x1 = x + offs, y0 = y0, y1 = y0)
  panel.segments(x0 = x - offs, x1 = x + offs, y0 = y1, y1 = y1)
}

#' Internal function for displaying error bars in lattice plots.
#'
#' \code{panel.errbars} is an internal function called by \code{bar_error} and \code{strip_error}.
#' @param x A numeric vector with the x positions for the bars.
#' @param y A numeric matrix with the mid values, lower values and upper values of the bars.
#' @param panel.xy A lattice panel.
#' @param make.grid Type of grid, a character with options: "horizontal", "vertical", "both" and "none".
#' @param ewidth An integer.
#' @param ... Further arguments passed to lattice.
panel.errbars <- function(x, y, ..., panel.xy = panel.xyplot, make.grid = c("horizontal", "vertical", "both", "none"), ewidth = 0)
{
   Y <- matrix(y, nrow = length(y) / 3, ncol = 3)
   y <- Y[, 1]
   y0 <- Y[, 2]
   y1 <- Y[, 3]
   make.grid <- match.arg(make.grid)
   if (make.grid == "horizontal") {
     panel.grid(h = -1, v = 0)
   } else if (make.grid == "vertical"){
     panel.grid(v = -1, h = 0)
   } else if (make.grid == "both"){
     panel.grid(v = -1, h = -1)
   }
   panel.errbars1(x, y0, y1, ewidth = ewidth)
   panel.xy(x, y, ...)
}

#' Generate a data frame with estimate and bootstrap CIs.
#'
#' \code{gen_bst_df} is an internal function called by \link{bar_error} that generates a data frame with
#' confidence intervals of a continuous variable by levels of one or two categorical ones (factors).
#'
#' @author Josie Athens, Department of Preventive and Social Medicine, University of Otago, New Zealand.
#' @param formula A formula of the form \code{y ~ x} or \code{y ~ x|z} where \code{y} is a
#' numerical variable and both \code{x} and \code{z} are factors.
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param stat Statistic used for \link{bst}.
#' @param ... Passes optional arguments to \link{bst}.
#' @return A data frame with the confidence intervals by level.
#' @examples
#' data(kfm, package = "ISwR")
#' gen_bst_df(dl.milk ~ sex, data = kfm)
#' bar_error(dl.milk ~ sex, data = kfm, ylim = c(0,9), ylab  = "Breast-milk intake (dl/day)")
#'
#' data(birthwt, package = "MASS")
#' birthwt$smoke <- factor(birthwt$smoke, labels=c("Non-smoker", "Smoker"))
#' birthwt$Race <- 0
#' birthwt$Race[birthwt$race>1] <- 1
#' birthwt$Race <- factor(birthwt$Race, labels=c("White", "Non-white"))
#' gen_bst_df(bwt ~ smoke|Race, data = birthwt)
#' bar_error(bwt ~ smoke|Race, data = birthwt, ylim = c(0, 3800), ylab  = "Birth weight (g)")
gen_bst_df <- function(formula, data, stat = "mean", ...)
{
	vars <- all.vars(formula)
	resp <- vars[1]
	expl <- vars[2]
  nv <- length(vars)
  if (nv == 2) {
    n <- length(levels(data[[expl]]))
  	first <- tapply(data[[resp]], data[[expl]], bst, stat, ...)
  	second <- numeric(n * 3)
  	dim(second) <- c(n, 3)
  	pos <- c(2, 4, 5)
  	for (i in 1:n) second[i,] <- as.numeric(first[[i]][pos])
  	df <- as.data.frame(second)
  	names(df) <- c(resp, "LowerCI", "UpperCI")
  	final  <- data.frame(df, factor(levels(data[[expl]]), levels = levels(data[[expl]])))
  	names(final)[4] <- expl
  	final
  } else {
    strat <- vars[3]
  	n1 <- length(levels(data[[expl]]))
  	n2 <- length(levels(data[[strat]]))
  	n <- n1*n2
  	first <- tapply(data[[resp]], survival::strata(data[[expl]], data[[strat]]), bst, stat, ...)
  	second <- numeric(n * 3)
  	dim(second) <- c(n, 3)
  	pos <- c(2, 4, 5)
  	for (i in 1:n) second[i,] <- as.numeric(first[[i]][pos])
  	df <- as.data.frame(second)
  	names(df) <- c(resp, "LowerCI", "UpperCI")
  	gp1 <- factor(rep(levels(data[[expl]]), each = n2), levels = levels(data[[expl]]))
  	gp2 <- factor(rep(levels(data[[strat]]), n1), levels = levels(data[[strat]]))
  	final  <- data.frame(df, gp1, gp2)
  	names(final)[4] <- expl
  	names(final)[5] <- strat
  	final
  }
}
