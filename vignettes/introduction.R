## ----eval=FALSE---------------------------------------------------------------
#  y ~ x, data = my_data

## ----eval=FALSE---------------------------------------------------------------
#  y ~ x|z, data = my_data

## ----message=FALSE, results='hide'--------------------------------------------
rm(list = ls())
library(car)
library(broom)
library(kableExtra)
library(tidyverse)
library(ggfortify)
library(mosaic)
library(huxtable)
library(jtools)
library(latex2exp)
library(pubh)
library(sjlabelled)
library(sjPlot)
library(sjmisc)

theme_set(sjPlot::theme_sjplot2(base_size = 10))
theme_update(legend.position = "top")
# options('huxtable.knit_print_df' = FALSE)
options('huxtable.knit_print_df_theme' = theme_article)
options('huxtable.autoformat_number_format' = list(numeric = "%5.2f"))
knitr::opts_chunk$set(collapse = TRUE, comment = NA)

## -----------------------------------------------------------------------------
data(Oncho)
Oncho %>% head

## -----------------------------------------------------------------------------
Oncho %>%
  cross_tab(mf ~ area)

## -----------------------------------------------------------------------------
Oncho %>%
  select(- c(id, mfload)) %>%
  cross_tab(mf ~ .)

## -----------------------------------------------------------------------------
data(Hodgkin)
Hodgkin <- Hodgkin %>%
  mutate(Ratio = CD4/CD8) %>%
  var_labels(
    Ratio = "CD4+ / CD8+ T-cells ratio"
    )

Hodgkin %>% head

## -----------------------------------------------------------------------------
Hodgkin %>%
  estat(~ CD4)

## -----------------------------------------------------------------------------
Hodgkin %>%
  estat(~ Ratio|Group)

## -----------------------------------------------------------------------------
Hodgkin %>%
  cross_tab(Group ~ .)

## -----------------------------------------------------------------------------
var.test(Ratio ~ Group, data = Hodgkin)

## -----------------------------------------------------------------------------
Hodgkin %>%
  qq_plot(~ Ratio|Group) %>%
  axis_labs()

## -----------------------------------------------------------------------------
wilcox.test(Ratio ~ Group, data = Hodgkin)

## -----------------------------------------------------------------------------
Hodgkin %>%
  strip_error(Ratio ~ Group) %>%
  axis_labs()

## -----------------------------------------------------------------------------
Hodgkin %>%
  strip_error(Ratio ~ Group) %>%
  axis_labs() %>%
  gf_star(x1 = 1, y1 = 4, x2 = 2, y2 = 4.05, y3 = 4.1, "**")

## -----------------------------------------------------------------------------
data(birthwt, package = "MASS")
birthwt <- birthwt %>%
  mutate(
    smoke = factor(smoke, labels = c("Non-smoker", "Smoker")),
    Race = factor(race > 1, labels = c("White", "Non-white")),
    race = factor(race, labels = c("White", "Afican American", "Other"))
    ) %>%
  var_labels(
    bwt = 'Birth weight (g)',
    smoke = 'Smoking status',
    race = 'Race',
  )

## -----------------------------------------------------------------------------
birthwt %>%
  bar_error(bwt ~ smoke) %>%
  axis_labs()

## -----------------------------------------------------------------------------
birthwt %>%
  qq_plot(~ bwt|smoke) %>%
  axis_labs()

## -----------------------------------------------------------------------------
t.test(bwt ~ smoke, data = birthwt)

## -----------------------------------------------------------------------------
birthwt %>%
  bar_error(bwt ~ smoke) %>%
  axis_labs() %>%
  gf_star(x1 = 1, x2 = 2, y1 = 3400, y2 = 3500, y3 = 3550, "**")

## -----------------------------------------------------------------------------
birthwt %>%
  bar_error(bwt ~ smoke, fill = ~ Race) %>%
  gf_refine(ggsci::scale_fill_jama()) %>%
  axis_labs()

## -----------------------------------------------------------------------------
birthwt %>%
  bar_error(bwt ~ smoke|Race) %>%
  axis_labs()

## -----------------------------------------------------------------------------
birthwt %>%
  strip_error(bwt ~ smoke, pch = ~ Race, col = ~ Race) %>%
  gf_refine(ggsci::scale_color_jama()) %>%
  axis_labs()

## -----------------------------------------------------------------------------
model_bwt <- lm(bwt ~ smoke + race, data = birthwt)

model_bwt %>%
  glm_coef(labels = model_labels(model_bwt))

## -----------------------------------------------------------------------------
tab_model(model_bwt,  collapse.ci = TRUE)

## -----------------------------------------------------------------------------
multiple(model_bwt, ~ race)$df

## -----------------------------------------------------------------------------
multiple(model_bwt, ~ race)$fig_ci %>%
  gf_labs(x = "Difference in birth weights (g)")

## -----------------------------------------------------------------------------
multiple(model_bwt, ~ race)$fig_pval %>%
  gf_labs(y = " ")

## -----------------------------------------------------------------------------
data(Bernard)
head(Bernard)

## -----------------------------------------------------------------------------
Bernard %>%
  cross_tab(treat ~ fate)

## -----------------------------------------------------------------------------
dat <- matrix(c(84, 140 , 92, 139), nrow = 2, byrow = TRUE)
epiR::epi.2by2(as.table(dat))

## -----------------------------------------------------------------------------
Bernard %>%
  contingency(fate ~ treat) 

## ----eval=FALSE---------------------------------------------------------------
#  outcome ~ stratum/exposure, data = my_data

## -----------------------------------------------------------------------------
data(oswego, package = "epitools")

oswego <- oswego %>%
  mutate(
    ill = factor(ill, labels = c("No", "Yes")),
    sex = factor(sex, labels = c("Female", "Male")),
    chocolate.ice.cream = factor(chocolate.ice.cream, labels = c("No", "Yes"))
  ) %>%
  var_labels(
    ill = "Developed illness",
    sex = "Sex",
    chocolate.ice.cream = "Consumed chocolate ice cream"
  )

## -----------------------------------------------------------------------------
oswego %>%
  cross_tab(ill ~ sex + chocolate.ice.cream, na_include = TRUE)

## -----------------------------------------------------------------------------
oswego %>%
  mhor(ill ~ sex/chocolate.ice.cream)

## -----------------------------------------------------------------------------
data(Oncho)

odds_trend(mf ~ agegrp, data = Oncho, angle = 0, hjust = 0.5)$fig

## -----------------------------------------------------------------------------
Freq <- c(1739, 8, 51, 22)
BCG <- gl(2, 1, 4, labels=c("Negative", "Positive"))
Xray <- gl(2, 2, labels=c("Negative", "Positive"))
tb <- data.frame(Freq, BCG, Xray)
tb <- expand_df(tb)

diag_test(BCG ~ Xray, data = tb)

## -----------------------------------------------------------------------------
diag_test2(22, 51, 8, 1739)

