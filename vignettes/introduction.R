## ----eval=FALSE---------------------------------------------------------------
#  y ~ x, data = my_data

## ----eval=FALSE---------------------------------------------------------------
#  y ~ x|z, data = my_data

## ----eval=FALSE---------------------------------------------------------------
#  data |>
#    f1(...) |>
#    f2(...) |>
#    f3(...)

## ----message=FALSE, results='hide'--------------------------------------------
rm(list = ls())
library(dplyr)
library(rstatix)
library(crosstable)
library(pubh)
library(sjlabelled)

## -----------------------------------------------------------------------------
data(Oncho)
Oncho |> head()

## -----------------------------------------------------------------------------
crosstable_options(
  total = "row",
  percent_pattern="{n} ({p_col})",
  percent_digits = 1,
  funs = c("Mean (std)" = meansd, "Median [IQR]" = mediqr)
)

## -----------------------------------------------------------------------------
Oncho |>
  select(mf, area) |>
  mutate(
    mf = relevel(mf, ref = "Infected")
  ) |>
  copy_labels(Oncho) |>
  crosstable(by = area) |>
  ctf()

## -----------------------------------------------------------------------------
Oncho |>
  select(- c(id, mfload)) |>
  mutate(
    mf = relevel(mf, ref = "Infected")
  ) |>
  copy_labels(Oncho) |>
  crosstable(by = area) |>
  ctf()

## -----------------------------------------------------------------------------
data(Hodgkin)
Hodgkin <- Hodgkin |>
  mutate(Ratio = CD4/CD8) |>
  var_labels(
    Ratio = "CD4+ / CD8+ T-cells ratio"
    )

Hodgkin |> head()

## -----------------------------------------------------------------------------
Hodgkin |>
  estat(~ CD4)

## -----------------------------------------------------------------------------
Hodgkin |>
  estat(~ Ratio|Group)

## -----------------------------------------------------------------------------
Hodgkin |>
  mutate(
    Group = relevel(Group, ref = "Hodgkin")
  ) |>
  copy_labels(Hodgkin) |>
  crosstable(by = Group) |>
  ctf()

## -----------------------------------------------------------------------------
var.test(Ratio ~ Group, data = Hodgkin)

## -----------------------------------------------------------------------------
Hodgkin |>
  qq_plot(~ Ratio|Group) 

## -----------------------------------------------------------------------------
wilcox.test(Ratio ~ Group, data = Hodgkin)

## -----------------------------------------------------------------------------
Hodgkin |>
  strip_error(Ratio ~ Group)

## -----------------------------------------------------------------------------
Hodgkin |>
  strip_error(Ratio ~ Group) |>
  gf_star(x1 = 1, y1 = 4, x2 = 2, y2 = 4.05, y3 = 4.1, "**")

## -----------------------------------------------------------------------------
data(birthwt, package = "MASS")
birthwt <- birthwt |>
  mutate(
    smoke = factor(smoke, labels = c("Non-smoker", "Smoker")),
    Race = factor(race > 1, labels = c("White", "Non-white")),
    race = factor(race, labels = c("White", "Afican American", "Other"))
    ) |>
  var_labels(
    bwt = 'Birth weight (g)',
    smoke = 'Smoking status',
    race = 'Race',
  )

## -----------------------------------------------------------------------------
birthwt |>
  bar_error(bwt ~ smoke)

## -----------------------------------------------------------------------------
birthwt |>
  qq_plot(~ bwt|smoke)

## -----------------------------------------------------------------------------
birthwt |> 
  t_test(bwt ~ smoke, detailed = TRUE) |> 
  as.data.frame()

## -----------------------------------------------------------------------------
birthwt |>
  bar_error(bwt ~ smoke) |>
  gf_star(x1 = 1, x2 = 2, y1 = 3400, y2 = 3500, y3 = 3550, "**")

## -----------------------------------------------------------------------------
birthwt |>
  bar_error(bwt ~ smoke, fill = ~ Race) 

## -----------------------------------------------------------------------------
birthwt |>
  bar_error(bwt ~ smoke|Race)

## -----------------------------------------------------------------------------
birthwt |>
  strip_error(bwt ~ smoke, pch = ~ Race, col = ~ Race)

## -----------------------------------------------------------------------------
model_bwt <- lm(bwt ~ smoke + race, data = birthwt)

## -----------------------------------------------------------------------------
model_bwt |>
  glm_coef(labels = model_labels(model_bwt))

## -----------------------------------------------------------------------------
multiple(model_bwt, ~ race)$df

## -----------------------------------------------------------------------------
multiple(model_bwt, ~ race)$fig_ci |>
  gf_labs(x = "Difference in birth weights (g)")

## -----------------------------------------------------------------------------
multiple(model_bwt, ~ race)$fig_pval |>
  gf_labs(y = " ")

