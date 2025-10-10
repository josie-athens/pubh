## ----eval=FALSE---------------------------------------------------------------
# y ~ x, data = my_data

## ----eval=FALSE---------------------------------------------------------------
# y ~ x|z, data = my_data

## ----eval=FALSE---------------------------------------------------------------
# data |>
#   f1(...) |>
#   f2(...) |>
#   f3(...)

## ----message=FALSE, results='hide'--------------------------------------------
rm(list = ls())
library(dplyr)
library(rstatix)
library(pubh)
library(sjlabelled)

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
  cross_tbl(by = "Group", bold = FALSE) |>
  theme_pubh(2) |>
  add_footnote("Median (IQR)", font_size = 9)

## -----------------------------------------------------------------------------
var.test(Ratio ~ Group, data = Hodgkin)

## -----------------------------------------------------------------------------
Hodgkin |>
  qq_plot(y = Ratio) +
  facet_grid(cols = vars(Group))

## -----------------------------------------------------------------------------
wilcox.test(Ratio ~ Group, data = Hodgkin)

## -----------------------------------------------------------------------------
Hodgkin |>
  strip_error(x = Group, y = Ratio)

## -----------------------------------------------------------------------------
p1 <- Hodgkin |>
  strip_error(x = Group, y = Ratio)

gg_star(p1, x1 = 1, y1 = 4, x2 = 2, y2 = 4.05, y3 = 4.1, "**") +
  ylim(0, 4.2)

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
  bar_error(x = smoke, y = bwt)

## -----------------------------------------------------------------------------
birthwt |>
  bar_error(x = smoke, y = bwt) +
  facet_grid(cols = vars(Race))

## -----------------------------------------------------------------------------
birthwt |>
  strip_error(x = smoke, y = bwt) +
  facet_grid(cols = vars(Race))

## -----------------------------------------------------------------------------
model_bwt <- lm(bwt ~ smoke + race, data = birthwt)

## -----------------------------------------------------------------------------
model_bwt |>
  glm_coef(labels = model_labels(model_bwt))

## -----------------------------------------------------------------------------
multiple(model_bwt, ~ race)

## -----------------------------------------------------------------------------
multiple(model_bwt, ~ race) |>
  effect_plot(x = estimate, y = contrast, orientation = "y") +
  geom_vline(xintercept = 0, lty = 2, col = bmj[2]) +
  xlab("Difference in birth weight (g)") + 
  ylab("")

