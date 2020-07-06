## ----message=FALSE, results = 'hide'------------------------------------------
rm(list = ls())
library(car)
library(broom)
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
data(birthwt, package = "MASS")
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
    Mean = mean(bwt, na.rm = TRUE),
    SD = sd(bwt, na.rm = TRUE),
    Median = median(bwt, na.rm = TRUE),
    CV = rel_dis(bwt)
  ) 

## -----------------------------------------------------------------------------
birthwt %>%
  gen_bst_df(bwt ~ race|smoke)

## -----------------------------------------------------------------------------
birthwt %>%
  bar_error(bwt ~ race, fill = ~ smoke) %>%
  axis_labs() %>%
  gf_labs(fill = "Smoking status:")

## -----------------------------------------------------------------------------
model_norm <- lm(bwt ~ smoke + race, data = birthwt)

## -----------------------------------------------------------------------------
model_norm %>% Anova() %>% tidy()

## -----------------------------------------------------------------------------
model_norm %>% tidy()

## -----------------------------------------------------------------------------
model_norm %>% 
  glm_coef(labels = model_labels(model_norm))

## -----------------------------------------------------------------------------
model_norm %>%
  glm_coef(se_rob = TRUE, labels = model_labels(model_norm))

## -----------------------------------------------------------------------------
model_norm %>% glance()

## -----------------------------------------------------------------------------
model_norm %>%
  plot_model("pred", terms = ~race|smoke, dot.size = 1.5, title = "")

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
diet %>% estat(~ fibre|chd)

## ----warning=FALSE------------------------------------------------------------
diet %>%
  gf_boxploth(chd ~ fibre, fill = "indianred3", alpha = 0.7) %>%
  axis_labs()

## -----------------------------------------------------------------------------
model_binom <- glm(chd ~ fibre, data = diet, family = binomial)

model_binom %>%
  glm_coef(labels = model_labels(model_binom))

## -----------------------------------------------------------------------------
model_binom %>%
  plot_model("pred", terms = "fibre [all]", title = "")

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
bdendo %>%
  mutate(
    cancer = relevel(cancer, ref = "Case"),
    est = relevel(est, ref = "Oestrogen"),
    gall = relevel(gall, ref = "GBD")
  ) %>%
  copy_labels(bdendo) %>%
  cross_tab(cancer ~ est + gall) %>%
  theme_article()

## -----------------------------------------------------------------------------
library(survival)
model_clogit <- clogit(cancer == 'Case'  ~ est * gall + strata(set), data = bdendo)

model_clogit %>%
  glm_coef(labels = c("Oestrogen/No oestrogen", "GBD/No GBD",
                      "Oestrogen:GBD Interaction")) 

## -----------------------------------------------------------------------------
require(ggeffects)
bdendo_pred <- ggemmeans(model_clogit, terms = c('gall', 'est'))

## -----------------------------------------------------------------------------
bdendo_pred %>%
  gf_pointrange(predicted + conf.low + conf.high ~ x|group, col = ~ x) %>%
  gf_labs(y = "P(cancer)", x = "", col = get_label(bdendo$gall))

## -----------------------------------------------------------------------------
library(ordinal)
data(housing, package = "MASS")
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
model_clm %>%
  glm_coef(labels = model_labels(model_clm, intercept = FALSE))

## -----------------------------------------------------------------------------
model_clm %>%
  plot_model(type = "pred", terms = c("Infl", "Cont"), 
           dot.size = 1, title = "") %>%
  gf_theme(axis.text.x = element_text(angle = 45, hjust = 1))

## -----------------------------------------------------------------------------
model_clm %>%
  plot_model(type = "pred", terms = c("Infl", "Type"), 
           dot.size = 1, title = "") %>%
  gf_theme(axis.text.x = element_text(angle = 45, hjust = 1))

## -----------------------------------------------------------------------------
emmip(model_clm, Infl ~ Type |Cont) %>%
  gf_labs(x = "Type of rental", col = "Perceived influence")

## -----------------------------------------------------------------------------
data(quine, package = "MASS")
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
    Mean = mean(Days, na.rm = TRUE),
    SD = sd(Days, na.rm = TRUE),
    Median = median(Days, na.rm = TRUE),
    CV = rel_dis(Days)
  ) 

## -----------------------------------------------------------------------------
model_pois <- glm(Days ~ Eth + Sex + Age, family = poisson, data = quine)

model_pois %>%
  glm_coef(labels = model_labels(model_pois), se_rob = TRUE)

## -----------------------------------------------------------------------------
model_pois %>% glance()

## -----------------------------------------------------------------------------
deviance(model_pois) / df.residual(model_pois)

## -----------------------------------------------------------------------------
library(MASS)
model_negbin <- glm.nb(Days ~ Eth + Sex + Age, data = quine)

model_negbin %>%
  glm_coef(labels = model_labels(model_negbin), se_rob = TRUE) 

## -----------------------------------------------------------------------------
model_negbin %>% glance()

## -----------------------------------------------------------------------------
model_negbin %>% Anova()

## -----------------------------------------------------------------------------
model_negbin %>%
  plot_model(type = "pred", terms = c("Age", "Eth"), 
           dot.size = 1.5, title = "") 

## -----------------------------------------------------------------------------
emmip(model_negbin, Eth ~ Age|Sex) %>%
  gf_labs(y = "Number of absent days", x = "Age group", col = "Ethnicity")

## -----------------------------------------------------------------------------
multiple(model_negbin, ~ Age|Eth)$df 

## -----------------------------------------------------------------------------
multiple(model_negbin, ~ Age|Eth)$fig_ci %>%
  gf_labs(x = "IRR")

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
  glm_coef(labels = c("Treatment: Thiotepa/Placebo", "Scale"), se_rob = TRUE)

## -----------------------------------------------------------------------------
model_exp <- survreg(Surv(times, event) ~ rx, data = bladder, dist = "exponential")

## -----------------------------------------------------------------------------
model_exp %>%
  glm_coef(labels = c("Treatment: Thiotepa/Placebo"), se_rob = TRUE)

## -----------------------------------------------------------------------------
model_exp %>%
  glm_coef(labels = c("Treatment: Thiotepa/Placebo"))

## -----------------------------------------------------------------------------
model_exp %>%
  plot_model(type = "pred", terms = ~ rx, dot.size = 1.5, title = "") %>%
  gf_labs(y = "Survival time", x = "Treatment", title = "")

## -----------------------------------------------------------------------------
model_cox <-  coxph(Surv(times, event) ~ rx, data = bladder)

## -----------------------------------------------------------------------------
model_cox %>%
  glm_coef(labels = c("Treatment: Thiotepa/Placebo"))

## -----------------------------------------------------------------------------
model_cox %>%
  plot_model(type = "pred", terms = ~ rx, dot.size = 1.5, 
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
  glm_coef(labels = c(
    "Constant", 
    "Sex: female-male", 
    "Age (years)",
    "Sex:Age interaction"
    )) 

## -----------------------------------------------------------------------------
model_lme %>%
  plot_model("pred", terms = age ~ Sex, 
           show.data = TRUE, jitter = 0.1, dot.size = 1.5) %>%
  gf_labs(y = get_label(Orthodont$distance), x = "Age (years)", title = "")

