## ----message=FALSE, results = 'hide'------------------------------------------
rm(list = ls())
library(dplyr)
library(rstatix)
library(easystats)
library(pubh)
library(sjlabelled)

theme_set(see::theme_lucid(base_size = 10))
theme_update(legend.position = "top")
options('huxtable.knit_print_df' = FALSE)
options('huxtable.autoformat_number_format' = list(numeric = "%5.2f"))
knitr::opts_chunk$set(comment = NA)

## -----------------------------------------------------------------------------
data(birthwt, package = "MASS")
birthwt <- birthwt |>
  mutate(
    smoke = factor(smoke, labels = c("Non-smoker", "Smoker")),
    race = factor(race, labels = c("White", "African American", "Other"))
    ) |>
  var_labels(
    bwt = 'Birth weight (g)',
    smoke = 'Smoking status',
    race = 'Race'
    )

## -----------------------------------------------------------------------------
birthwt |>
  box_plot(bwt ~ smoke, fill = ~ race)

## -----------------------------------------------------------------------------
birthwt |>
  gen_bst_df(bwt ~ race|smoke) |>
  as_hux() |> theme_pubh()

## -----------------------------------------------------------------------------
model_norm <- lm(bwt ~ smoke + race, data = birthwt)

## -----------------------------------------------------------------------------
model_norm |> Anova()

## -----------------------------------------------------------------------------
model_norm |> parameters()

## -----------------------------------------------------------------------------
model_norm |> performance()

## -----------------------------------------------------------------------------
model_norm |>
  tbl_regression() |> 
  cosm_reg() |> theme_pubh() |> 
  add_footnote(get_r2(model_norm), font_size = 9)

## -----------------------------------------------------------------------------
model_norm |> 
  glm_coef(labels = model_labels(model_norm)) |>
  as_hux() |> set_align(everywhere, 2:3, "right") |>
  theme_pubh() |>
  add_footnote(get_r2(model_norm), font_size = 9)

## -----------------------------------------------------------------------------
model_norm |>
  glm_coef(se_rob = TRUE, labels = model_labels(model_norm)) |>
  as_hux() |> set_align(everywhere, 2:3, "right") |>
  theme_pubh() |>
  add_footnote(paste(
    get_r2(model_norm), "\n",
    "CIs and p-values estimated with robust standard errors."),
    font_size = 9)

## -----------------------------------------------------------------------------
model_norm |>
  estimate_means(c("race", "smoke"))

## -----------------------------------------------------------------------------
model_norm |> 
  estimate_means(by = c("race", "smoke")) |> 
  plot() |> 
  gf_labs(
    x = "", y = "Birth weight (g)", title = "",
    col = "Smoking status"
  )

## -----------------------------------------------------------------------------
model_norm |>
  estimate_means(c("race", "smoke")) |>
  gf_point(Mean ~ race) |>
  gf_errorbar(CI_low + CI_high ~ race, width = 0.3) |>
  gf_facet_wrap(~smoke) |>
  gf_labs(
    x = "Race",
    y = "Birth weight (g)"
  )

## -----------------------------------------------------------------------------
emmip(model_norm, smoke ~ race) |>
  gf_labs(y = get_label(birthwt$bwt), x = "", col = "Smoking status")

## -----------------------------------------------------------------------------
data(diet, package = "Epi")
diet <- diet |>
  mutate(
    chd = factor(chd, labels = c("No CHD", "CHD"))
  ) |>
  var_labels(
    chd = "Coronary Heart Disease",
    fibre = "Fibre intake (10 g/day)"
    )

## -----------------------------------------------------------------------------
diet |> estat(~ fibre|chd) |>
  as_hux() |> theme_pubh()

## -----------------------------------------------------------------------------
diet |> na.omit() |>
  copy_labels(diet) |>
  box_plot(fibre ~ chd)

## -----------------------------------------------------------------------------
model_binom <- glm(chd ~ fibre, data = diet, family = binomial)

## -----------------------------------------------------------------------------
model_binom |> parameters(exponentiate = TRUE)

## -----------------------------------------------------------------------------
model_binom |> performance()

## -----------------------------------------------------------------------------
model_binom |> 
  tbl_regression(exponentiate = TRUE) |> 
  cosm_reg() |> theme_pubh() |> 
  add_footnote(get_r2(model_binom), font_size = 9)

## -----------------------------------------------------------------------------
model_binom |>
  estimate_means(by = "fibre") |>
  gf_line(Probability ~ fibre) |>
  gf_ribbon(CI_low + CI_high ~ fibre) |>
  gf_labs(
    x = "Fibre intake (10 g/day)",
    y = "Probability of CHD"
  )

## -----------------------------------------------------------------------------
data(bdendo, package = "Epi") 
bdendo <- bdendo |>
  mutate(
    cancer = factor(d, labels = c('Control', 'Case')),
    gall = factor(gall, labels = c("No GBD", "GBD")),
    est = factor(est, labels = c("No oestrogen", "Oestrogen"))
  ) |>
  var_labels(
    cancer = 'Endometrial cancer',
    gall = 'Gall bladder disease',
    est = 'Oestrogen'
  )

## -----------------------------------------------------------------------------
bdendo |> 
  mutate(
    cancer = relevel(cancer, ref = "Case"),
    gall = relevel(gall, ref = "GBD"),
    est = relevel(est, ref = "Oestrogen")
  ) |>
  copy_labels(bdendo) |>
  select(cancer, gall, est) |> 
  tbl_strata(
    strata = est,
    .tbl_fun = ~ .x |>
      tbl_summary(by = gall)
  ) |> 
  cosm_sum(bold = TRUE, head_label = " ") |> 
  theme_pubh(2) |> 
  set_align(1, everywhere, "center")

## -----------------------------------------------------------------------------
bdendo |>
  gf_percents(~ cancer|gall, fill = ~ est, position = "dodge", alpha = 0.6) |>
  gf_labs(
    y = "Percent",
    x = "",
    fill = ""
  )

## -----------------------------------------------------------------------------
require(survival, quietly = TRUE)
model_clogit <- clogit(cancer == 'Case' ~ est * gall + strata(set), data = bdendo)

## -----------------------------------------------------------------------------
model_clogit |>
  glm_coef(labels = c("Oestrogen/No oestrogen", "GBD/No GBD",
                      "Oestrogen:GBD Interaction")) |>
  as_hux() |> set_align(everywhere, 2:3, "right") |>
  theme_pubh() |>
  add_footnote(get_r2(model_clogit), font_size = 9)

## ----message=FALSE------------------------------------------------------------
#require(MASS, quietly = TRUE)
data(housing, package = "MASS")
housing <- housing |>
  var_labels(
    Sat = "Satisfaction",
    Infl = "Perceived influence",
    Type = "Type of rental",
    Cont = "Contact"
    )

## -----------------------------------------------------------------------------
model_ord <- MASS::polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing, Hess = TRUE)

## -----------------------------------------------------------------------------
model_ord |> 
  tbl_regression(exponentiate = TRUE) |> 
  cosm_reg() |> theme_pubh() |> 
  add_footnote(get_r2(model_ord), font_size = 9)

## -----------------------------------------------------------------------------
model_ord |>
  estimate_means(by = c("Infl", "Cont", "Type")) |>
  mutate(
    Mean = inv_logit(Mean),
    CI_low = inv_logit(CI_low),
    CI_high = inv_logit(CI_high)
  ) |>
  gf_errorbar(CI_low + CI_high ~ Infl, width = 0.2) |>
  gf_point(Mean ~ Infl, size = 0.7) |>
  gf_facet_grid(~ Type + Cont) |>
  gf_labs(
    x = "Perceived influence",
    y = "Satisfaction"
  )

## -----------------------------------------------------------------------------
emmip(model_ord, Infl ~ Type |Cont) |>
  gf_labs(x = "Type of rental", col = "Perceived influence")

## -----------------------------------------------------------------------------
data(quine, package = "MASS")
levels(quine$Eth) <- c("Aboriginal", "White")
levels(quine$Sex) <- c("Female", "Male")

## -----------------------------------------------------------------------------
quine <- quine |>
  var_labels(
    Days = "Number of absent days",
    Eth = "Ethnicity",
    Age = "Age group"
    )

## -----------------------------------------------------------------------------
quine |> 
  cross_tbl(by = "Eth") |> 
  theme_pubh(2) |> 
  add_footnote("n (%); Median (IQR)", font_size = 9)

## -----------------------------------------------------------------------------
quine |>
  box_plot(Days ~ Age|Sex, fill = ~ Eth)

## -----------------------------------------------------------------------------
model_pois <- glm(Days ~ Eth + Sex + Age, family = poisson, data = quine)

## -----------------------------------------------------------------------------
model_pois |> 
  tbl_regression(exponentiate = TRUE) |> 
  cosm_reg() |> theme_pubh() |> 
  add_footnote(get_r2(model_pois), font_size = 9)

## -----------------------------------------------------------------------------
model_pois |> performance()

## -----------------------------------------------------------------------------
model_pois |> check_overdispersion()

## -----------------------------------------------------------------------------
model_negbin <- MASS::glm.nb(Days ~ Eth + Sex + Age, data = quine)

## -----------------------------------------------------------------------------
model_negbin |> 
  tbl_regression(exponentiate = TRUE) |> 
  cosm_reg() |> theme_pubh() |> 
  add_footnote(get_r2(model_negbin), font_size = 9)

## -----------------------------------------------------------------------------
model_negbin |> performance()

## -----------------------------------------------------------------------------
model_negbin |>
  estimate_means(by = c("Age", "Eth", "Sex")) |>
  gf_errorbar(CI_low + CI_high ~ Age, width = 0.2) |>
  gf_point(Mean ~ Age, size = 1) |>
  gf_facet_wrap(~ Eth + Sex) |>
  gf_labs(
    x = "Age group",
    y = "Number of absent days"
  )

## -----------------------------------------------------------------------------
emmip(model_negbin, Eth ~ Age|Sex) |>
  gf_labs(y = "Number of absent days", x = "Age group", col = "Ethnicity")

## -----------------------------------------------------------------------------
multiple(model_negbin, ~ Age|Eth)$df

## -----------------------------------------------------------------------------
multiple(model_negbin, ~ Age|Eth)$fig_ci |>
  gf_labs(x = "IRR")

## -----------------------------------------------------------------------------
data(bladder)
bladder <- bladder |>
  mutate(times = stop,
         rx = factor(rx, labels=c("Placebo", "Thiotepa"))
         ) |>
  var_labels(times = "Survival time",
             rx = "Treatment")

## -----------------------------------------------------------------------------
model_surv <- survreg(Surv(times, event) ~ rx, data = bladder)

## -----------------------------------------------------------------------------
model_surv |>
  glm_coef(labels = c("Treatment: Thiotepa/Placebo", "Scale"), se_rob = TRUE) |>
  as_hux() |> set_align(everywhere, 2:3, "right") |>
  theme_pubh(1) |>
  add_footnote(get_r2(model_surv), font_size = 9)

## -----------------------------------------------------------------------------
model_exp <- survreg(Surv(times, event) ~ rx, data = bladder, dist = "exponential")

## -----------------------------------------------------------------------------
model_exp |>
  glm_coef(labels = c("Treatment: Thiotepa/Placebo"), se_rob = TRUE) |>
  as_hux() |> set_align(everywhere, 2:3, "right") |>
  theme_pubh() |>
  add_footnote(get_r2(model_exp), font_size = 9)

## -----------------------------------------------------------------------------
model_exp |>
  glm_coef(labels = c("Treatment: Thiotepa/Placebo")) |>
  as_hux() |> set_align(everywhere, 2:3, "right") |>
  theme_pubh() |>
  add_footnote(get_r2(model_exp), font_size = 9)

## -----------------------------------------------------------------------------
model_exp |>
  estimate_means(by = "rx") |>
  gf_errorbar(CI_low + CI_high ~ rx, width = 0.3) |>
  gf_point(Mean ~ rx, size = 1) |>
  gf_labs(
    x = "Treatment",
    y = "Survival time"
  )

## -----------------------------------------------------------------------------
model_cox <-  coxph(Surv(times, event) ~ rx, data = bladder)

## -----------------------------------------------------------------------------
model_cox |> 
  tbl_regression(exponentiate = TRUE) |> 
  cosm_reg() |> theme_pubh() |> 
  add_footnote(get_r2(model_cox), font_size = 9)

## -----------------------------------------------------------------------------
model_cox |>
  estimate_means(by = "rx")

