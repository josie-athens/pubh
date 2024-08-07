% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{box_plot}
\alias{box_plot}
\title{Construct box plots.}
\usage{
box_plot(
  object = NULL,
  formula = NULL,
  data = NULL,
  fill = "indianred3",
  alpha = 0.7,
  outlier.shape = 20,
  outlier.size = 1,
  ...
)
}
\arguments{
\item{object}{When chaining, this holds an object produced in the earlier portions of the chain. Most users can safely ignore this argument. See details and examples.}

\item{formula}{A formula with shape: \code{y ~ x} where \code{y} is a numerical variable and \code{x} is a factor.}

\item{data}{A data frame where the variables in the \code{formula} can be found.}

\item{fill}{Colour used for the box passed to \code{\link[ggformula]{gf_boxplot}}.}

\item{alpha}{Opacity (0 = invisible, 1 = opaque).}

\item{outlier.shape}{Shape (\code{pch}) used as symbol for the outliers.}

\item{outlier.size}{Size of the outlier symbol.}

\item{...}{Further arguments passed to \code{\link[ggformula]{gf_boxplot}}.}
}
\description{
\code{box_plot} is a wrap function that calls \code{gf_boxplot} to construct more aesthetic box plots.
}
\examples{
data(kfm, package = "ISwR")
require(sjlabelled, quietly = TRUE)
kfm <- kfm |>
  var_labels(
    dl.milk = "Breast-milk intake (dl/day)",
    sex = "Sex",
    weight = "Child weight (kg)",
    ml.suppl = "Milk substitute (ml/day)",
    mat.weight = "Maternal weight (kg)",
    mat.height = "Maternal height (cm)"
  )

kfm |>
  box_plot(dl.milk ~ sex, fill = "thistle", alpha = 0.8)

t.test(dl.milk ~ sex, data = kfm)

kfm |>
  box_plot(dl.milk ~ sex, fill = "thistle", alpha = 0.8) |>
  gf_star(1, 10.9, 2, 11, 11.4, legend = "p = 0.035", size = 2.5)
}
