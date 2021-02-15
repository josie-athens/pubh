## ----message=FALSE, results = 'hide'------------------------------------------
rm(list = ls())
library(car)
library(broom)
library(mosaic)
library(tidyverse)
library(ggfortify)
library(huxtable)
library(jtools)
library(latex2exp)
library(pubh)
library(sjlabelled)
library(sjPlot)
library(sjmisc)

theme_set(sjPlot::theme_sjplot2(base_size = 10))
theme_update(legend.position = "top")
options('huxtable.knit_print_df' = FALSE)
options('huxtable.autoformat_number_format' = list(numeric = "%5.2f"))
knitr::opts_chunk$set(comment = NA)

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
  ) %>%
  as_hux() %>% theme_pubh(1)

## -----------------------------------------------------------------------------
birthwt %>%
  box_plot(bwt ~ smoke, fill = ~ race) %>%
  axis_labs()

## -----------------------------------------------------------------------------
birthwt %>%
  gen_bst_df(bwt ~ race|smoke) %>%
  as_hux() %>% theme_pubh(1)

## -----------------------------------------------------------------------------
model_norm <- lm(bwt ~ smoke + race, data = birthwt)

## -----------------------------------------------------------------------------
model_norm %>% Anova()

## -----------------------------------------------------------------------------
model_norm %>% 
  summ(confint = TRUE, model.info = FALSE)

## -----------------------------------------------------------------------------
model_norm %>% 
  glm_coef(labels = model_labels(model_norm)) %>%
  as_hux() %>% set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(get_r2(model_norm), font_size = 9)

## -----------------------------------------------------------------------------
model_norm %>%
  glm_coef(se_rob = TRUE, labels = model_labels(model_norm)) %>%
  as_hux() %>% set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(paste(
    get_r2(model_norm), "\n",
    "CIs and p-values estimated with robust standard errors."),
    font_size = 9)

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
diet %>% estat(~ fibre|chd) %>%
  as_hux() %>% theme_pubh(1)

## -----------------------------------------------------------------------------
diet %>% na.omit() %>%
  copy_labels(diet) %>%
  box_plot(fibre ~ chd) %>%
  axis_labs()

## -----------------------------------------------------------------------------
model_binom <- glm(chd ~ fibre, data = diet, family = binomial)

## -----------------------------------------------------------------------------
model_binom %>% 
  summ(confint = TRUE, model.info = FALSE, exp = TRUE)

## -----------------------------------------------------------------------------
model_binom %>%
  glm_coef(labels = model_labels(model_binom)) %>%
  as_hux() %>% set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(get_r2(model_binom), font_size = 9)

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
t1 = bdendo %>%
  filter(est == "Oestrogen") %>%
  mutate(
    cancer = relevel(cancer, ref = "Case"),
    gall = relevel(gall, ref = "GBD")
  ) %>%
  copy_labels(bdendo) %>%
  cross_tab(cancer ~ gall,
            label = "Oestrogen users")

t2 = bdendo %>%
  filter(est == "No oestrogen") %>%
  mutate(
    cancer = relevel(cancer, ref = "Case"),
    gall = relevel(gall, ref = "GBD")
  ) %>%
  copy_labels(bdendo) %>%
  cross_tab(cancer ~ gall,
            label = "Non-oestrogen users")

rbind(t1, t2) %>%
  theme_pubh(c(3, 6, 9))

## -----------------------------------------------------------------------------
bdendo %>%
  gf_percents(~ cancer|gall, fill = ~ est, position = "dodge", alpha = 0.6) %>%
  gf_labs(
    y = "Percent",
    x = "",
    fill = ""
  )

## -----------------------------------------------------------------------------
require(survival, quietly = TRUE)
model_clogit <- clogit(cancer == 'Case' ~ est * gall + strata(set), data = bdendo)

model_clogit %>%
  glm_coef(labels = c("Oestrogen/No oestrogen", "GBD/No GBD",
                      "Oestrogen:GBD Interaction")) %>%
  as_hux() %>% set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(get_r2(model_clogit), font_size = 9)

## -----------------------------------------------------------------------------
require(ggeffects, quietly = TRUE)
bdendo_pred <- ggemmeans(model_clogit, terms = c('gall', 'est'))

## -----------------------------------------------------------------------------
bdendo_pred %>%
  gf_pointrange(predicted + conf.low + conf.high ~ x|group, col = ~ x) %>%
  gf_labs(y = "P(cancer)", x = "", col = get_label(bdendo$gall))

## ----message=FALSE------------------------------------------------------------
require(MASS, quietly = TRUE)
data(housing)
housing <- housing %>%
  var_labels(
    Sat = "Satisfaction",
    Infl = "Perceived influence",
    Type = "Type of rental",
    Cont = "Contact"
    )

## -----------------------------------------------------------------------------
model_ord <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing, Hess = TRUE)

## -----------------------------------------------------------------------------
model_ord %>%
  glm_coef(labels = model_labels(model_ord, intercept = FALSE)) %>%
  as_hux() %>% set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(get_r2(model_ord), font_size = 9)

## -----------------------------------------------------------------------------
model_ord %>%
  plot_model(type = "pred", terms = c("Infl", "Cont"), 
           dot.size = 1, title = "") %>%
  gf_theme(axis.text.x = element_text(angle = 45, hjust = 1))

## -----------------------------------------------------------------------------
model_ord %>%
  plot_model(type = "pred", terms = c("Infl", "Type"), 
           dot.size = 1, title = "") %>%
  gf_theme(axis.text.x = element_text(angle = 45, hjust = 1))

## -----------------------------------------------------------------------------
emmip(model_ord, Infl ~ Type |Cont) %>%
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

## ----message=FALSE------------------------------------------------------------
quine %>%
  group_by(Eth, Sex, Age) %>%
  summarise(
    n = n(),
    Mean = mean(Days, na.rm = TRUE),
    SD = sd(Days, na.rm = TRUE),
    Median = median(Days, na.rm = TRUE),
    CV = rel_dis(Days)
  ) %>%
  as_hux() %>% theme_pubh(1)

## -----------------------------------------------------------------------------
quine %>%
  box_plot(Days ~ Age|Sex, fill = ~ Eth) %>%
  axis_labs() %>% gf_labs(fill = "")

## -----------------------------------------------------------------------------
model_pois <- glm(Days ~ Eth + Sex + Age, family = poisson, data = quine)

model_pois %>%
  glm_coef(labels = model_labels(model_pois), se_rob = TRUE) %>%
  as_hux() %>% set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(get_r2(model_pois), font_size = 9)

## -----------------------------------------------------------------------------
model_pois %>% glance()

## -----------------------------------------------------------------------------
deviance(model_pois) / df.residual(model_pois)

## ----message=FALSE------------------------------------------------------------
model_negbin <- glm.nb(Days ~ Eth + Sex + Age, data = quine)

model_negbin %>%
  glm_coef(labels = model_labels(model_negbin), se_rob = TRUE)  %>%
  as_hux() %>% set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(get_r2(model_negbin), font_size = 9)

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
multiple(model_negbin, ~ Age|Eth)$df %>%
  as_hux() %>% theme_pubh(1)

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
  glm_coef(labels = c("Treatment: Thiotepa/Placebo", "Scale"), se_rob = TRUE) %>%
  as_hux() %>% set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(get_r2(model_surv), font_size = 9)

## -----------------------------------------------------------------------------
model_exp <- survreg(Surv(times, event) ~ rx, data = bladder, dist = "exponential")

## -----------------------------------------------------------------------------
model_exp %>%
  glm_coef(labels = c("Treatment: Thiotepa/Placebo"), se_rob = TRUE) %>%
  as_hux() %>% set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(get_r2(model_exp), font_size = 9)

## -----------------------------------------------------------------------------
model_exp %>%
  glm_coef(labels = c("Treatment: Thiotepa/Placebo")) %>%
  as_hux() %>% set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(get_r2(model_exp), font_size = 9)

## -----------------------------------------------------------------------------
model_exp %>%
  plot_model(type = "pred", terms = ~ rx, dot.size = 1.5, title = "") %>%
  gf_labs(y = "Survival time", x = "Treatment", title = "")

## -----------------------------------------------------------------------------
model_cox <-  coxph(Surv(times, event) ~ rx, data = bladder)

## -----------------------------------------------------------------------------
model_cox %>%
  glm_coef(labels = c("Treatment: Thiotepa/Placebo")) %>%
  as_hux() %>% set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(get_r2(model_cox), font_size = 9)

## -----------------------------------------------------------------------------
model_cox %>%
  plot_model(type = "pred", terms = ~ rx, dot.size = 1.5, 
           title = "") %>%
  gf_labs(x = "Treatment", title = "")

## ---- message=FALSE-----------------------------------------------------------
require(lme4, quietly = TRUE)
data(Orthodont, package = "nlme")

## -----------------------------------------------------------------------------
Orthodont <- Orthodont %>%
  var_labels(
    distance = "Pituitary distance (mm)",
    age = "Age (years)"
    )

## -----------------------------------------------------------------------------
model_mix <- lmer(distance ~ Sex * age + (1|Subject), data = Orthodont)

## -----------------------------------------------------------------------------
model_mix %>%
  summ(center = TRUE, confint = TRUE, model.info = FALSE)

## -----------------------------------------------------------------------------
model_mix %>%
  center_mod() %>%
  glm_coef(labels = c(
    model_labels(model_mix),
    "Sex:Age interaction"
    )) %>%
  as_hux() %>% set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(get_r2(model_mix), font_size = 9)

## -----------------------------------------------------------------------------
model_mix %>%
  plot_model("eff", terms = age ~ Sex, 
           show.data = TRUE, jitter = 0.1, dot.size = 1.5) %>%
  gf_labs(y = get_label(Orthodont$distance), x = "Age (years)", title = "")

