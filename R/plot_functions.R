# Functions for Plotting

#' Palette inspired by the British Medical Journal.
bmj = c(
  "#2A6EBBFF", "#F0AB00FF", "#C50084FF", "#7D5CC6FF",
  "#E37222FF", "#69BE28FF", "#00B2A9FF", "#CD202CFF",
  "#747678FF"
)

#' Annotating a plot to display differences between groups.
#'
#' \code{gg_star} is a function used to display differences between groups (see details).
#'
#' @param fig A gformula object.
#' @param x1 Position in \code{x} for the start of the horizontal line.
#' @param y1 Position in \code{y} for the start of the vertical line, below to the horizontal line.
#' @param x2 Position in \code{x} for the end of the horizontal line.
#' @param y2 Position in \code{y} where the horizontal line is drawn.
#' @param y3 Position in \code{y} where the text is added.
#' @param legend Character text used for annotating the plot.
#' @param ... Additional information passed to \code{\link[ggformula]{gf_text}}.
#' @details This function draws an horizontal line from coordinate (\code{x1}, \code{y2})
#' to coordinate (\code{x2}, \code{y2}). Draws vertical lines below the horizontal line,
#' towards data, from (\code{x1}, \code{y1}) to (\code{x1}, \code{y2}) and from
#' (\code{x2}, \code{y1}) to (\code{x2}, \code{y2}). Finally, adds
#' text above the horizontal line, at the mid point between \code{x1} and \code{x2}.
#' See examples.
#' @examples
#' data(kfm, package = "ISwR")
#' require(sjlabelled, quietly = TRUE)
#' data(energy, package = "ISwR")
#' energy <- energy |>
#'   var_labels(
#'     expend = "Energy expenditure (MJ/day)",
#'     stature = "Stature"
#'   )
#'
#' p <- energy |>
#'   strip_error(x = stature, y = expend)
#'
#' gg_star(p, x1 = 1, y1 = 13, x2 = 2, y2 = 13.2, y3 = 13.4, "**") +
#' ylim(0, 13.6)
#' @export
gg_star <- function(fig, x1, y1, x2, y2, y3, legend = "*", ...) {
  fig + ggplot2::geom_segment(
      x = x1, y = y2, xend = x2, yend = y2,
      linewidth = 0.3, col = 'black'
    ) +
    ggplot2::geom_segment(
      x = x1, y = y1, xend = x1, yend = y2,
      linewidth = 0.3, col = 'black'
    ) +
    ggplot2::geom_segment(
      x = x2, y = y1, xend = x2, yend = y2,
      linewidth = 0.3, col = 'black'
    ) +
    ggplot2::geom_text(
      x = mean(c(x1, x2)), y = y3, label = legend, ...
    )
}

#' Multiple comparisons with plot.
#'
#' \code{multiple} displays results from post-doc analysis and constructs corresponding plot.
#'
#' @param model A fitted model supported by \code{emmeans}, such as the result of a call to \code{aov}, \code{lm}, \code{glm}, etc.
#' @param formula A formula with shape: \code{~ y} or \code{~ y|x} (for interactions). Where \code{y} is the term of the model on which comparisons are made and \code{x} is a term interacting with \code{y}.
#' @param adjust Method to adjust CIs and p-values (see details).
#' @param type Type of prediction  (matching "linear.predictor", "link", or "response").
#' @param reverse Logical argument. Determines the direction of comparisons.
#' @param level Confidence interval significance level.
#' @param digits Number of digits for rounding (default = 2).
#' @param ... Further arguments passed to \code{\link[emmeans]{emmeans}}.
#' @details
#' The default adjusting method is "mvt" which uses the multivariate t distribution.
#' Other options are: "bonferroni", "holm", "hochberg", "tukey" and "none". The default option for argument \code{reverse} is to make reverse comparisons, i.e., against the reference level matching comparisons from \code{lm} and \code{glm}.
#' @return A list with objects: \code{df} A data frame with adjusted p-values, \code{fig_ci} a plot with estimates and adjusted confidence intervals and \code{fig_pval} a plot comparing adjusted p-values.
#' @seealso \code{\link[emmeans]{emmeans}}, \code{\link[emmeans]{pwpp}}.
#' @examples
#' data(birthwt, package = "MASS")
#' birthwt$race <- factor(birthwt$race, labels = c("White", "African American", "Other"))
#'
#' model_1 <- aov(bwt ~ race, data = birthwt)
#' multiple(model_1, ~race)
#'
#' multiple(model_1, ~race) |>
#'   effect_plot(x = estimate, y = contrast, orientation = "y") +
#'   xlab("Difference in birth weights (g)") +
#'   ylab("") + geom_vline(xintercept = 0, col = bmj[2], lty = 2)
#' @export
multiple <- function(
  model, formula, adjust = "hochberg", type = "response", reverse = TRUE, level = 0.95, digits = 2, ...
) {
  term_emm <- emmeans::emmeans(model, formula, type = type, ...)
  emm_df <- as.data.frame(pairs(term_emm, adjust = adjust, reverse = reverse))
  emm_ci <- as.data.frame(confint(pairs(term_emm, adjust = adjust, reverse = reverse), level = level))
  log10_pval <- log10(emm_df[["p.value"]])
  emm_df$p.value <- round_pval(emm_df[["p.value"]])
  emm_df <- emm_df |>
    dplyr::select(-df)
  emm_ci <- emm_ci |>
    dplyr::select(-df)
  vars <- all.vars(formula)
  n <- length(vars)
  if (n == 1) {
    names(emm_ci)[4:5] <- c("CI_low", "CI_high")
  } else {
    names(emm_ci)[5:6] <- c("CI_low", "CI_high")
  }
  emm_df <- cbind(emm_df, CI_low = emm_ci[["CI_low"]], CI_high = emm_ci[["CI_high"]])
  emm_df <- sjmisc::round_num(emm_df, digits = digits)
  emm_df
}

#' Bland-Altman agreement plots.
#'
#' @details \code{bland_altman} constructs Bland-Altman agreement plots.
#' @details Variables in \code{formula} are continuous paired observations. When the distribution of the outcome
#' is not normal, but becomes normal with a log-transformation, \code{bland_altman} can plot the ratio between
#' outcomes (difference in the log scale) by using option \code{transform = TRUE}.
#'
#' @param data A data frame.
#' @param x A numerical variable.
#' @param y A numerical variable.
#' @param size Size of the symbol using to plot data.
#' @param transform Logical, should ratios instead of differences be used to construct the plot?
#' @examples
#' data(Sharples)
#'
#' Sharples |>
#'   bland_altman(
#'     x = weight, y = srweight,
#'     transform = TRUE, size = 0.7
#'   ) +
#'   xlab("Mean of weights (kg)") +
#'   ylab("Measured / Self-reported weight")
#' @export
bland_altman <- function(
  data, x, y, size = 1, transform = FALSE
) {
  expl <- deparse(substitute(x))
  resp <- deparse(substitute(y))

  diff <- data[[expl]] - data[[resp]]
  ratio <- data[[expl]] / data[[resp]]
  avg <- apply(cbind(data[[expl]], data[[resp]]), 1, mean, na.rm = TRUE)
  df <- data.frame(Difference = diff, Ratio = ratio, Mean = avg)
  lw1 <- as.numeric(reference_range(
    mean(diff, na.rm = TRUE),
    sd(diff, na.rm = TRUE)
  ))[1]
  up1 <- as.numeric(reference_range(
    mean(diff, na.rm = TRUE),
    sd(diff, na.rm = TRUE)
  ))[2]
  lw2 <- as.numeric(reference_range(
    mean(ratio, na.rm = TRUE),
    sd(ratio, na.rm = TRUE)
  ))[1]
  up2 <- as.numeric(reference_range(
    mean(ratio, na.rm = TRUE),
    sd(ratio, na.rm = TRUE)
  ))[2]
  if (transform == TRUE) {
    ggplot2::ggplot(
      data = df, ggplot2::aes(x = df[["Mean"]], y = df[["Ratio"]])
    ) +
    ggplot2::geom_point(size = size) +
    ggplot2::geom_hline(yintercept = mean(ratio), col = bmj[2]) +
    ggplot2::geom_hline(yintercept = lw2, col = bmj[2], lty = 2) +
    ggplot2::geom_hline(yintercept = up2, col = bmj[2], lty = 2)
  } else {
    ggplot2::ggplot(
      data = df, ggplot2::aes(x = df[["Mean"]], y = df[["Difference"]])
    ) +
    ggplot2::geom_point(size = size) +
    ggplot2::geom_hline(yintercept = mean(diff), col = bmj[2]) +
    ggplot2::geom_hline(yintercept = lw1, col = bmj[2], lty = 2) +
    ggplot2::geom_hline(yintercept = up1, col = bmj[2], lty = 2)
  }
}

#' Quantile-quantile plots against the standard Normal distribution.
#'
#' \code{qq_plot} constructs quantile-quantile plots against the standard normal distribution
#' (also known as quantile-normal plots).
#' @param data A data frame.
#' @param y A numerical variable in the data frame.
#' @param size Size of the marker.
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
#'   qq_plot(dl.milk)
#' @export
qq_plot <- function(
  data, y, size = 0.7
) {
    yl <- deparse(substitute(y))
    ggplot2::ggplot(
      data = data, ggplot2::aes(sample = {{y}})
    ) +
    ggplot2::geom_qq(size = size) +
    ggplot2::geom_qq_line(col = bmj[2]) +
    ggplot2::xlab("Theoretical quantiles") +
    ggplot2::ylab(get_label(data[[yl]]))
}

#' Strip plots with error bars.
#'
#' \code{strip_error} constructs strip plots with error bars showing 95% bootstrapped
#' confidence intervals around mean values.
#' @param data A data frame.
#' @param x A categorical variable (the exposure).
#' @param y A numerical variable (the outcome or response variable).
#' @param size Size of the marker.
#' @examples
#' data(energy, package = "ISwR")
#' require(sjlabelled, quietly = TRUE)
#' energy <- energy |>
#'   var_labels(
#'     expend = "Energy expenditure (MJ/day)",
#'     stature = "Stature"
#'   )
#'
#' energy |>
#'   strip_error(x = stature, y = expend)
#' @export
strip_error <- function(data, x, y, size = 0.5){
  ggplot2::ggplot(
    data = data, ggplot2::aes(x = {{x}}, y = {{y}})
  ) +
    ggplot2::geom_jitter(height = 0, width = 0.1, size = size) +
    ggplot2::stat_summary(
      fun.data = "mean_cl_boot", geom = "errorbar",
      color = bmj[2], width = 0.2
    ) +
    ggplot2::stat_summary(
      fun = "mean", geom = "point", color = bmj[2],
      size = 1
    )
}

#' Bar charts with error bars.
#'
#' \code{bar_error} constructs bar charts in with error bars showing 95% bootstrapped
#' confidence intervals around mean values. High of bars represent mean values.
#' @param data A data frame.
#' @param x A variable in the data frame.
#' @param y A variable in the data frame
#' @param alpha Opacity of the colour fill (0 = invisible, 1 = opaque).
#' @examples
#' require(dplyr, quietly = TRUE)
#' require(sjlabelled, quietly = TRUE)
#' data(birthwt, package = "MASS")
#'
#' birthwt <- birthwt |>
#'   mutate(
#'     smoke = factor(smoke, labels = c("Non-smoker", "Smoker"))
#'   ) |>
#'   var_labels(
#'     bwt = "Birth weight (g)",
#'     smoke = "Smoking status"
#'   )
#'
#' birthwt |>
#'   bar_error(x = smoke, y = bwt)
#' @export
bar_error <- function(
  data, x, y, alpha = 0.7
) {
  ggplot2::ggplot(
    data = data, ggplot2::aes(x = {{x}}, y = {{y}})
  ) +
    ggplot2::stat_summary(
      geom = "bar", fun = "mean", alpha = alpha
    ) +
    ggplot2::stat_summary(
      fun.data = "mean_cl_boot", geom = "errorbar",
      width = 0.2, color = "black"
    )
}

#' Generate a data frame with estimate and bootstrap CIs.
#'
#' \code{gen_bst_df} is a function called that generates a data frame with
#' confidence intervals of a continuous variable by levels of one or two categorical ones (factors).
#' @param object When chaining, this holds an object produced in the earlier portions of the chain. Most users can safely ignore this argument. See details and examples.
#' @param formula A formula with shape: \code{y ~ x} or \code{y ~ x|z} where \code{y} is a
#' numerical variable and both \code{x} and \code{z} are factors.
#' @param data A data frame where the variables in the \code{formula} can be found.
#' @param stat Statistic used for \link{bst}.
#' @param ... Passes optional arguments to \link{bst}.
#' @return A data frame with the confidence intervals by level.
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
#'   gen_bst_df(dl.milk ~ sex)
#'
#' data(Fentress)
#'
#' Fentress |> gen_bst_df(pain ~ group)
#' @export
gen_bst_df <- function (object = NULL, formula = NULL, data = NULL, stat = "mean", ...)
{
  if (inherits(object, "formula"))
  {
    formula <- object
    object <- NULL
  }
  if (inherits(object, "data.frame")) {
    data <- object
    object <- NULL
  }
  vars <- all.vars(formula)
  resp <- vars[1]
  expl <- vars[2]
  outcome <- data[[resp]]
  exposure <- data[[expl]]
  nv <- length(vars)
  if (nv == 2) {
    n <- length(levels(exposure))
    first <- tapply(outcome, exposure, bst, stat, ...)
    second <- numeric(n * 3)
    dim(second) <- c(n, 3)
    pos <- c(2, 4, 5)
    for (i in 1:n) second[i, ] <- as.numeric(first[[i]][pos])
    df <- as.data.frame(second)
    names(df) <- c(resp, "CI_low", "CI_high")
    final <- data.frame(df, factor(levels(exposure), levels = levels(exposure)))
    resp <- ifelse(is.null(sjlabelled::get_label(outcome)),
                   resp, sjlabelled::get_label(outcome))
    expl <- ifelse(is.null(sjlabelled::get_label(exposure)),
                   expl, sjlabelled::get_label(exposure))
    names(final)[1] <- resp
    names(final)[4] <- expl
    final
  }
  else {
    condit <- vars[3]
    strat <- data[[condit]]
    n1 <- length(levels(exposure))
    n2 <- length(levels(strat))
    n <- n1 * n2
    first <- tapply(outcome, survival::strata(exposure, strat),
                    bst, stat, ...)
    second <- numeric(n * 3)
    dim(second) <- c(n, 3)
    pos <- c(2, 4, 5)
    for (i in 1:n) second[i, ] <- as.numeric(first[[i]][pos])
    df <- as.data.frame(second)
    names(df) <- c(resp, "CI_low", "CI_high")
    gp1 <- factor(rep(levels(exposure), each = n2), levels = levels(exposure))
    gp2 <- factor(rep(levels(strat), n1), levels = levels(strat))
    final <- data.frame(df, gp1, gp2)
    resp <- ifelse(is.null(sjlabelled::get_label(outcome)),
                   resp, sjlabelled::get_label(outcome))
    expl <- ifelse(is.null(sjlabelled::get_label(exposure)),
                   expl, sjlabelled::get_label(exposure))
    strat <- ifelse(is.null(sjlabelled::get_label(strat)),
                    expl, sjlabelled::get_label(strat))
    names(final)[1] <- resp
    names(final)[4] <- expl
    names(final)[5] <- condit
    final
  }
}

#' Constructs effect plots on for cases when the response (outcome) is a numerical variable and the explanatory (exposure) is a categorical variable.
#'
#' \code{effect_plot} constructs plots with error bars on which the upper and lower limits
#' of the bars are given by the variables \code{CI_high} and \code{CI_low}, respectively.
#' @param data A data frame.
#' @param x A variable in the data frame.
#' @param y A variable in the data frame.
#' @param orientation A string indicating the orientation of the bars.
#' @return An error bars plot, for example, those to show the effects of a categorical variable
#' on a numerical variable.
#' @examples
#' data(Fentress)
#' require(dplyr, quietly = TRUE)
#' require(sjlabelled, quietly = TRUE)
#'
#' pain_cis <- Fentress |>
#'   gen_bst_df(pain ~ group)
#' names(pain_cis)[1] = "Mean"
#' pain_cis = pain_cis |>
#'   var_labels(Mean = "Pain reduction")
#'
#' pain_cis |>
#'   effect_plot(x = Cohort, y = Mean)
#' @export
effect_plot <- function(
  data, x, y, orientation = "x"
){
  if(orientation == "x"){
    ggplot2::ggplot(
      data=data, ggplot2::aes(x={{x}}, y={{y}})
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = data[["CI_low"]], ymax = data[["CI_high"]]),
      width = 0.2
    ) +
    ggplot2::geom_point()
  } else {
    ggplot2::ggplot(
      data=data, ggplot2::aes(x={{x}}, y={{y}})
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = data[["CI_low"]], xmax = data[["CI_high"]]),
      width = 0.2
    ) +
    ggplot2::geom_point()
  }

}
