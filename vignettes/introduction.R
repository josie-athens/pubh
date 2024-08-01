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
library(pubh)
library(sjlabelled)

theme_set(see::theme_lucid(base_size = 10))
theme_update(legend.position = "top")
options('huxtable.knit_print_df' = FALSE)
options('huxtable.autoformat_number_format' = list(numeric = "%5.2f"))
knitr::opts_chunk$set(comment = NA)

## -----------------------------------------------------------------------------
data(Oncho)
Oncho |> head()

## -----------------------------------------------------------------------------
Oncho |>
  mutate(
    mf = relevel(mf, ref = "Infected")
  ) |>
  copy_labels(Oncho) |>
  select(mf, area) |> 
  cross_tbl(by = "area") |>
  theme_pubh(2)

## -----------------------------------------------------------------------------
Oncho |>
  select(- c(id, mfload)) |>
  mutate(
    mf = relevel(mf, ref = "Infected")
  ) |>
  copy_labels(Oncho) |>
  cross_tbl(by = "area") |>
  theme_pubh(2)

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
  estat(~ CD4) |>
  as_hux() |> theme_pubh()

## -----------------------------------------------------------------------------
Hodgkin |>
  estat(~ Ratio|Group) |>
  as_hux() |> theme_pubh()

## -----------------------------------------------------------------------------
Hodgkin |>
  mutate(
    Group = relevel(Group, ref = "Hodgkin")
  ) |>
  copy_labels(Hodgkin) |>
  cross_tbl(by = "Group", bold = FALSE) |>
  theme_pubh(2) |>
  add_footnote("Median (IQR)", font_size = 9)

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
  bar_error(bwt ~ smoke, fill = ~ Race) |>
  gf_refine(ggsci::scale_fill_jama()) 

## -----------------------------------------------------------------------------
birthwt |>
  bar_error(bwt ~ smoke|Race)

## -----------------------------------------------------------------------------
birthwt |>
  strip_error(bwt ~ smoke, pch = ~ Race, col = ~ Race) |>
  gf_refine(ggsci::scale_color_jama())

## -----------------------------------------------------------------------------
model_bwt <- lm(bwt ~ smoke + race, data = birthwt)

## -----------------------------------------------------------------------------
model_bwt |> 
  tbl_regression() |> 
  cosm_reg() |> theme_pubh() |> 
  add_footnote(get_r2(model_bwt), font_size = 9)

## -----------------------------------------------------------------------------
model_bwt |>
  glm_coef(labels = model_labels(model_bwt)) |>
  as_hux() |> theme_pubh() |> 
  set_align(everywhere, 2:3, "right") |>
  add_footnote(get_r2(model_bwt), font_size = 9)

## -----------------------------------------------------------------------------
multiple(model_bwt, ~ race)$df

## -----------------------------------------------------------------------------
multiple(model_bwt, ~ race)$fig_ci |>
  gf_labs(x = "Difference in birth weights (g)")

## -----------------------------------------------------------------------------
multiple(model_bwt, ~ race)$fig_pval |>
  gf_labs(y = " ")

## -----------------------------------------------------------------------------
data(Bernard)
Bernard |> head()

## -----------------------------------------------------------------------------
Bernard |>
  mutate(
    fate = relevel(fate, ref = "Dead"),
    treat = relevel(treat, ref = "Ibuprofen")
  ) |>
  copy_labels(Bernard) |>
  select(fate, treat) |> 
  cross_tbl(by = "treat") |>
  theme_pubh(2)

## -----------------------------------------------------------------------------
Bernard |>
  contingency(fate ~ treat) 

## ----eval=FALSE---------------------------------------------------------------
#  outcome ~ stratum/exposure, data = my_data

## -----------------------------------------------------------------------------
data(oswego, package = "epitools")

oswego <- oswego |>
  mutate(
    ill = factor(ill, labels = c("No", "Yes")),
    sex = factor(sex, labels = c("Female", "Male")),
    chocolate.ice.cream = factor(chocolate.ice.cream, labels = c("No", "Yes"))
  ) |>
  var_labels(
    ill = "Developed illness",
    sex = "Sex",
    chocolate.ice.cream = "Consumed chocolate ice cream"
  )

## -----------------------------------------------------------------------------
oswego |>
  mhor(ill ~ sex/chocolate.ice.cream)

## -----------------------------------------------------------------------------
Freq <- c(1739, 8, 51, 22)
BCG <- gl(2, 1, 4, labels=c("Negative", "Positive"))
Xray <- gl(2, 2, labels=c("Negative", "Positive"))
tb <- data.frame(Freq, BCG, Xray)
tb <- expand_df(tb)

diag_test(BCG ~ Xray, data = tb)

## -----------------------------------------------------------------------------
diag_test2(22, 51, 8, 1739)

