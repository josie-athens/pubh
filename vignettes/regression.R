## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = NA, size = "small", 
                      message = FALSE, warning = FALSE)

## ---- message=FALSE-----------------------------------------------------------
library(broom)
library(car)
library(Hmisc)
library(MASS)
library(kableExtra)
library(tidyverse)
library(mosaic)
library(latex2exp)
library(moonBook)
library(pubh)
library(sjlabelled)
library(sjPlot)

theme_set(sjPlot::theme_sjplot2(base_size = 10))
theme_update(legend.position = "top")
options(knitr.table.format = "pandoc")

## -----------------------------------------------------------------------------
data(birthwt)
birthwt <- birthwt %>%
  mutate(
    smoke = factor(smoke, labels = c("Non-smoker", "Smoker")),
    race = factor(race, labels = c("White", "African American", "Other"))
    ) %>%
  var_labels(
    bwt = 'Birth weight (g)',
    smoke = 'Smoking status',
    race = 'Race'
    )

## -----------------------------------------------------------------------------
birthwt %>%
  group_by(race, smoke) %>%
  summarise(
    n = n(),
    Mean = round(mean(bwt, na.rm = TRUE), 2),
    SD = round(sd(bwt, na.rm = TRUE), 2),
    Median = round(median(bwt, na.rm = TRUE), 2),
    CV = round(rel_dis(bwt), 2)
  ) %>%
  kable(caption = "Descriptive statistics of birth weight (g) by ethnicity
        and smoking status.")

## -----------------------------------------------------------------------------
birthwt %>%
  gen_bst_df(bwt ~ race|smoke) %>%
  kable(caption = "Mean birth weight (g) by ethnicity
        and smoking status with 95% CIs.")

## -----------------------------------------------------------------------------
birthwt %>%
  bar_error(bwt ~ race, fill = ~ smoke) %>%
  axis_labs() %>%
  gf_labs(fill = "Smoking status:")

## -----------------------------------------------------------------------------
model_norm <- lm(bwt ~ smoke + race, data = birthwt)

## -----------------------------------------------------------------------------
model_norm %>%
  Anova() 

## -----------------------------------------------------------------------------
model_norm %>%
  summary()

## -----------------------------------------------------------------------------
model_norm %>%
  glm_coef() 

## -----------------------------------------------------------------------------
model_norm %>%
  glm_coef(se_rob = FALSE, labels = c("Constant",
                                      "Smoker - No smoker",
                                      "African American - White",
                                      "Other - White")) %>%
  kable(caption = "Table of coefficients using naive standard errors.",
        align = 'r')

## -----------------------------------------------------------------------------
model_norm %>%
  glm_coef(labels = c("Constant",
                      "Smoker - No smoker",
                      "African American - White",
                      "Other - White")) %>%
  kable(caption = "Table of coefficients using robust standard errors.")

## -----------------------------------------------------------------------------
model_norm %>% glance %>% round(digits = 2)

## -----------------------------------------------------------------------------
plot_model(model_norm, "pred", terms = ~race|smoke, dot.size = 2, title = "")

## -----------------------------------------------------------------------------
emmip(model_norm, smoke ~ race) %>%
  gf_labs(y = get_label(birthwt$bwt), x = "", col = "Smoking status")

## -----------------------------------------------------------------------------
data(diet, package = "Epi")
diet <- diet %>%
  mutate(
    chd = factor(chd, labels = c("No CHD", "CHD"))
  ) %>%
  var_labels(
    chd = "Coronary Heart Disease",
    fibre = "Fibre intake (10 g/day)"
    )

## -----------------------------------------------------------------------------
diet %>% estat(fibre ~ chd) %>% kable

## -----------------------------------------------------------------------------
diet %>%
  gf_boxploth(chd ~ fibre, fill = "indianred3", alpha = 0.7) %>%
  axis_labs()

## -----------------------------------------------------------------------------
model_binom <- glm(chd ~ fibre, data = diet, family = binomial)
model_binom %>%
  glm_coef(labels = c("Constant", "Fibre intake (g/day)")) %>%
  kable(caption = "Parameter estimates from logistic regression.", align = "r")

## -----------------------------------------------------------------------------
model_binom %>% glance %>% round(digits = 2)

## -----------------------------------------------------------------------------
plot_model(model_binom, "pred", terms = "fibre [all]", title = "")

## -----------------------------------------------------------------------------
data(bdendo, package = "Epi") 
bdendo <- bdendo %>%
  mutate(
    cancer = factor(d, labels = c('Control', 'Case')),
    gall = factor(gall, labels = c("No GBD", "GBD")),
    est = factor(est, labels = c("No oestrogen", "Oestrogen"))
  ) %>%
  var_labels(
    cancer = 'Endometrial cancer',
    gall = 'Gall bladder disease',
    est = 'Oestrogen'
  )

## -----------------------------------------------------------------------------
mytable(cancer ~ est + gall, data = bdendo, show.total = TRUE) %>%
  mytable2df() %>%
  kable()

## -----------------------------------------------------------------------------
model_clogit <- clogit(cancer == 'Case'  ~ est * gall + strata(set), data = bdendo)
model_clogit %>%
  glm_coef(labels = c("Oestrogen/No oestrogen", "GBD/No GBD",
                      "Oestrogen:GBD Interaction")) %>%
  kable()

## -----------------------------------------------------------------------------
model_clogit %>% glance %>% round(digits = 2)

## -----------------------------------------------------------------------------
require(ggeffects)
bdendo_pred <- ggemmeans(model_clogit, terms = c('gall', 'est'))

## -----------------------------------------------------------------------------
bdendo_pred %>%
  gf_pointrange(predicted + conf.low + conf.high ~ x|group, col = ~ x) %>%
  gf_labs(y = "P(cancer)", x = "", col = get_label(bdendo$gall))

## -----------------------------------------------------------------------------
library(ordinal, warn.conflicts = FALSE)
data(housing)
housing <- housing %>%
  var_labels(
    Sat = "Satisfaction",
    Infl = "Perceived influence",
    Type = "Type of rental",
    Cont = "Contact"
    )

## -----------------------------------------------------------------------------
model_clm <- clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)

## -----------------------------------------------------------------------------
labs_ord <- c("Constant: Low/Medium satisfaction",
              "Constant: Medium/High satisfaction",
              "Perceived influence: Medium/Low",
              "Perceived influence: High/Low",
              "Accommodation: Apartment/Tower",
              "Accommodation: Atrium/Tower",
              "Accommodation: Terrace/Tower",
              "Afforded: High/Low")
kable(glm_coef(model_clm, labels = labs_ord), 
      caption = "Parameter estimates on satisfaction of householders.")

## -----------------------------------------------------------------------------
model_clm %>% glance %>% round(digits = 2)

## -----------------------------------------------------------------------------
plot_model(model_clm, dot.size = 1, title = "")

## -----------------------------------------------------------------------------
plot_model(model_clm, type = "pred", terms = c("Infl", "Cont"), 
           dot.size = 1, title = "") %>%
  gf_theme(axis.text.x = element_text(angle = 45, hjust = 1))

## -----------------------------------------------------------------------------
plot_model(model_clm, type = "pred", terms = c("Infl", "Type"), 
           dot.size = 1, title = "") %>%
  gf_theme(axis.text.x = element_text(angle = 45, hjust = 1))

## -----------------------------------------------------------------------------
emmip(model_clm, Infl ~ Type |Cont) %>%
  gf_labs(x = "Type of rental", col = "Perceived influence")

## -----------------------------------------------------------------------------
data(quine)
levels(quine$Eth) <- c("Aboriginal", "White")
levels(quine$Sex) <- c("Female", "Male")

## -----------------------------------------------------------------------------
quine <- quine %>%
  var_labels(
    Days = "Number of absent days",
    Eth = "Ethnicity",
    Age = "Age group"
    )

## -----------------------------------------------------------------------------
quine %>%
  group_by(Eth, Sex, Age) %>%
  summarise(
    n = n(),
    Mean = round(mean(Days, na.rm = TRUE), 2),
    SD = round(sd(Days, na.rm = TRUE), 2),
    Median = round(median(Days, na.rm = TRUE), 2)
  ) %>%
  kable()

## -----------------------------------------------------------------------------
model_pois <- glm(Days ~ Eth + Sex + Age, family = poisson, data = quine)
glm_coef(model_pois) %>% kable

## -----------------------------------------------------------------------------
model_pois %>% glance %>% round(digits = 2)

## -----------------------------------------------------------------------------
deviance(model_pois) / df.residual(model_pois)

## -----------------------------------------------------------------------------
model_negbin <- glm.nb(Days ~ Eth + Sex + Age, data = quine)
glm_coef(model_negbin, 
         labels=c(
           "Constant",
           "Race: Aboriginal/White",
           "Sex: Female/Male",
           "F1/Primary",
           "F2/Primary",
           "F3/Primary")
         ) %>%
  kable()

## -----------------------------------------------------------------------------
model_negbin %>% glance %>% round(digits = 2)

## -----------------------------------------------------------------------------
model_negbin %>%
  Anova() 

## -----------------------------------------------------------------------------
plot_model(model_negbin, type = "pred", terms = c("Age", "Eth"), 
           dot.size = 1, title = "") 

## -----------------------------------------------------------------------------
emmip(model_negbin, Eth ~ Age|Sex) %>%
  gf_labs(y = "Number of absent days", x = "Age group", col = "Ethnicity")

## -----------------------------------------------------------------------------
multiple(model_negbin, ~ Age|Eth)$df 

## -----------------------------------------------------------------------------
multiple(model_negbin, ~ Age + Sex|Eth)$fig_ci %>%
  gf_labs(y = "Age group", x = "Number of absent days")

## -----------------------------------------------------------------------------
data(bladder)
bladder <- bladder %>%
  mutate(times = stop,
         rx = factor(rx, labels=c("Placebo", "Thiotepa"))
         ) %>%
  var_labels(times = "Survival time",
             rx = "Treatment")

## -----------------------------------------------------------------------------
model_surv <- survreg(Surv(times, event) ~ rx, data = bladder)

## -----------------------------------------------------------------------------
model_surv %>%
  glm_coef(labels = c("Treatment: Thiotepa/Placebo", "Scale")) %>%
  kable()

## -----------------------------------------------------------------------------
model_exp <- survreg(Surv(times, event) ~ rx, data = bladder, dist = "exponential")
model_exp %>%
  glm_coef(labels = c("Treatment: Thiotepa/Placebo")) %>%
  kable()

## -----------------------------------------------------------------------------
model_exp %>% glance %>% round(digits = 2)

## -----------------------------------------------------------------------------
model_exp %>%
  glm_coef(labels = c("Treatment: Thiotepa/Placebo"), se_rob = FALSE) %>%
  kable()

## -----------------------------------------------------------------------------
plot_model(model_exp, type = "pred", terms = ~ rx, dot.size = 1, title = "") %>%
  gf_labs(y = "Survival time", x = "Treatment", title = "")

## -----------------------------------------------------------------------------
model_cox <-  coxph(Surv(times, event) ~ rx, data = bladder)

## -----------------------------------------------------------------------------
model_cox %>%
  glm_coef(labels = c("Treatment: Thiotepa/Placebo")) %>%
  kable()

## -----------------------------------------------------------------------------
model_cox %>% glance %>% round(digits = 2)

## -----------------------------------------------------------------------------
plot_model(model_cox, type = "pred", terms = ~ rx, dot.size = 1, 
           title = "") %>%
  gf_labs(x = "Treatment", title = "")

## -----------------------------------------------------------------------------
library(nlme)
data(Orthodont)

## -----------------------------------------------------------------------------
Orthodont <- Orthodont %>%
  var_labels(
    distance = "Pituitary distance (mm)",
    age = "Age (years)"
    )

## -----------------------------------------------------------------------------
model_lme <- lme(distance ~ Sex * I(age - mean(age, na.rm = TRUE)), random = ~ 1|Subject, 
                 method = "ML", data = Orthodont)

## -----------------------------------------------------------------------------
model_lme %>%
  glm_coef(labels = c("Constant", "Sex: female-male", "Age (years)",
                      "Sex:Age interaction")) %>%
  kable()

## -----------------------------------------------------------------------------
model_lme %>% glance

## -----------------------------------------------------------------------------
plot_model(model_lme, type = "pred", terms = age ~ Sex) %>%
  gf_labs(y = get_label(Orthodont$distance), x = "Age (years)", title = "")

