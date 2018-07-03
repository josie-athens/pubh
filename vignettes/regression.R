## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE)

## ---- message=FALSE------------------------------------------------------
library(pubh, warn.conflicts = FALSE)
library(car, warn.conflicts = FALSE)
library(descr, warn.conflicts = FALSE)
library(effects, warn.conflicts = FALSE)
library(multcomp, warn.conflicts = FALSE)
library(pander, warn.conflicts = FALSE)

set.alignment("right", row.names = "left", permanent = TRUE)

## ---- message=FALSE------------------------------------------------------
data(birthwt)
birthwt$smoke <- factor(birthwt$smoke, labels=c("Non-smoker", "Smoker"))
birthwt$race <- factor(birthwt$race > 1, labels=c("White", "Non-white"))
model_norm <- glm(bwt ~ smoke + race, data = birthwt)

## ------------------------------------------------------------------------
pander(Anova(model_norm))

## ------------------------------------------------------------------------
pander(summary(model_norm))

## ------------------------------------------------------------------------
glm_coef(model_norm)

## ------------------------------------------------------------------------
pander(glm_coef(model_norm, labels=c("Constant", "Smoker - Non-smoker", "Non-white - White"),
                se.rob = FALSE), split.table=Inf, caption="Table of coeffients using naive 
       standard errors.")

## ------------------------------------------------------------------------
pander(glm_coef(model_norm, labels=c("Constant", "Smoker - Non-smoker", "Non-white - White")),
       split.table = Inf, caption="Table of coeffients using robust standard errors.")

## ------------------------------------------------------------------------
plot(Effect(c("smoke", "race"), model_norm), multiline = TRUE, main = NULL, 
     ylab = "Birth weight (g)", xlab = "Smoking status", symbols = list(pch=16),
     confint = list(style="auto"), aspect = 3/4, lines = list(col=c(2,4), lwd=1.5))

## ------------------------------------------------------------------------
data(diet, package = "Epi")
model_binom <- glm(chd ~ fibre, data = diet, family = binomial)

## ------------------------------------------------------------------------
pander(glm_coef(model_binom, labels = c("Constant", "Fibre intake (g/day)")), split.table = Inf,
       caption = "Parameter estimates from logistic regression.")

## ------------------------------------------------------------------------
plot(Effect("fibre", model_binom), type = "response", rug = FALSE, aspect = 3/4,
       ylab = "P (CHD)", xlab = "Fibre (g/day)", lwd = 2, confint = list(style="none"),
     main = NULL)

## ------------------------------------------------------------------------
data(bdendo, package = "Epi")
levels(bdendo$gall) <- c("No GBD", "GBD")
levels(bdendo$est) <- c("No oestrogen", "Oestrogen")

## ------------------------------------------------------------------------
model_clogit <- clogit(d ~ est * gall + strata(set), data = bdendo)
glm_coef(model_clogit)

## ------------------------------------------------------------------------
pander(glm_coef(model_clogit, labels = c("Oestrogen/No oestrogen", "GBD/No GBD", 
                                         "Oestrogen:GBD Interaction")), 
       split.table = Inf, caption = "Parameter estimates from conditional logistic regression.")

## ------------------------------------------------------------------------
bdendo_grid <- with(bdendo, expand.grid(
  gall = levels(gall),
  est = levels(est),
  set = sample(1:63, 1)
))

## ------------------------------------------------------------------------
bdendo_grid$pred <- inv_logit(predict(model_clogit, bdendo_grid, type = "lp"))

## ------------------------------------------------------------------------
xyplot(pred  ~ gall, data = bdendo_grid, groups = est, type = "l", 
    lwd = 2, xlab = "Gall blader disease", ylab = "P (cancer)", 
    auto.key = list(title = "Oestrogen", space = "right", cex = 0.8))

## ------------------------------------------------------------------------
library(ordinal, warn.conflicts = FALSE)
data(housing)
model_clm <- clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
glm_coef(model_clm)

## ------------------------------------------------------------------------
labs_ord <- c("Constant: Low/Medium satisfaction",
              "Constant: Medium/High satisfaction",
              "Perceived influence: Medium/Low",
              "Perceived influence: High/Low",
              "Accommodation: Apartment/Tower",
              "Accommodation: Atrium/Tower",
              "Accommodation: Terrace/Tower",
              "Afforded: High/Low")
pander(glm_coef(model_clm, labels = labs_ord), split.table = Inf,
       caption = "Parameter estimates on satisfaction of householders.")

## ------------------------------------------------------------------------
plot(Effect(c("Infl", "Type", "Cont"), model_clm), main = NULL, aspect = 3/4, rotx = 45, 
     ylab = "Satisfaction (probability)", lines = list(lwd = 1.5, multiline = TRUE),
     confint = list(style="none"), symbols = list(pch = rep(20, 3)),
     ylim = c(0, 1))

## ---- message=FALSE------------------------------------------------------
library(nnet)
model_multi <- multinom(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)

## ------------------------------------------------------------------------
glm_coef(model_multi)

## ------------------------------------------------------------------------
plot(Effect(c("Infl", "Type", "Cont"), model_multi), main = NULL, aspect = 3/4, rotx = 45, 
     ylab = "Satisfaction (probability)", lines = list(lwd = 1.5, multiline = TRUE),
     confint = list(style="none"), symbols = list(pch = rep(20, 3)),
     ylim = c(0, 1))

## ------------------------------------------------------------------------
data(quine)
levels(quine$Eth) <- list(White = "N", Aboriginal = "A")
levels(quine$Sex) <- list(Male = "M", Female = "F")
model_pois <- glm(Days ~ Eth + Sex + Age, family = poisson, data = quine)
glm_coef(model_pois)

## ------------------------------------------------------------------------
pander(estat(~ Days|Eth, data = quine, label = "Days of school absences"), split.table=Inf)

## ------------------------------------------------------------------------
deviance(model_pois) / df.residual(model_pois)

## ------------------------------------------------------------------------
model_negbin <- glm.nb(Days ~ Eth + Sex + Age, data = quine)
unadj <- glm_coef(model_negbin, labels=c("Constant",
                                   "Race: Aboriginal/White",
                                   "Sex: Female/Male",
                                   "F1/Primary",
                                   "F2/Primary",
                                   "F3/Primary"))

## ------------------------------------------------------------------------
pander(Anova(model_negbin))

## ------------------------------------------------------------------------
pander(unadj, split.table = Inf, caption = "Parameter estimates with unadjusted CIs and p-values.")

## ------------------------------------------------------------------------
plot(Effect(c("Age", "Eth"), model_negbin), lines = list(lwd = 1.5, multiline = TRUE),
     confint = list(style="none"), symbols = list(pch = rep(20, 2)), main = NULL, 
     aspect = 3/4)

## ---- message=FALSE------------------------------------------------------
model_glht <- glht(model_negbin, linfct  = mcp(Age = "Tukey"))
age_glht <- xymultiple(model_glht, Exp = TRUE, plot = FALSE)

## ---- fig.cap="Parameter estimates on the effect of age group on the number of days absent from school. Bars represent 95% CIs adjusted by the method of Westfall for multiple comparisons."----
xymultiple(model_glht, Exp = TRUE)

## ------------------------------------------------------------------------
final <- unadj
final[, 5] <- as.character(final[, 5])
age_glht[, 5] <- as.character(age_glht[, 5])
final[4:6, 3:5] <- age_glht[1:3, 3:5]

## ------------------------------------------------------------------------
pander(final, split.table = Inf, caption = "Parameter estimates. CIs and
       p-values for age group were adjusted for multiple comparisons by the 
       method of Westfall.")

## ------------------------------------------------------------------------
data(bladder)
bladder$times <- bladder$stop
bladder$rx <- factor(bladder$rx, labels=c("Placebo", "Thiotepa"))

## ------------------------------------------------------------------------
model_surv <- survreg(Surv(times, event) ~ rx, data = bladder)

## ------------------------------------------------------------------------
glm_coef(model_surv)

## ------------------------------------------------------------------------
pander(glm_coef(model_surv, labels = c("Treatment: Thiotepa/Placebo", "Scale")),
       split.table = Inf)

## ------------------------------------------------------------------------
model_exp <- survreg(Surv(times, event) ~ rx, data = bladder, dist = "exponential")
pander(glm_coef(model_exp, labels = "Treatment: Thiotepa/Placebo"),
       split.table = Inf)

## ------------------------------------------------------------------------
pander(glm_coef(model_exp, se.rob = FALSE, labels = "Treatment: Thiotepa/Placebo"),
       split.table = Inf)

## ------------------------------------------------------------------------
bladder_grid <- with(bladder, expand.grid(
  rx = levels(rx)
))

## ------------------------------------------------------------------------
bladder_pred <- predict(model_exp, bladder_grid, se.fit = TRUE, type = "response")
bladder_grid$fit <- bladder_pred$fit
bladder_grid$se <- bladder_pred$se
bladder_grid$lo <- bladder_pred$fit - 1.96 * bladder_pred$se
bladder_grid$up <- bladder_pred$fit + 1.96 * bladder_pred$se

## ------------------------------------------------------------------------
xyplot(cbind(fit, lo, up) ~ rx, data = bladder_grid, pch = 20, panel = panel.errbars,
       ylab = "Survival time", xlab = "Treatment", aspect = 3/4)

## ------------------------------------------------------------------------
model_cox <-  coxph(Surv(times, event) ~ rx, data = bladder)

## ------------------------------------------------------------------------
pander(glm_coef(model_cox, labels = c("Treatment: Thiotepa/Placebo")), split.table = Inf)

## ------------------------------------------------------------------------
cox_grid <- with(bladder, expand.grid(
  rx = levels(rx)
))

## ------------------------------------------------------------------------
cox_pred <- predict(model_cox, cox_grid, se.fit = TRUE, type = "risk")
cox_grid$fit <- cox_pred$fit
cox_grid$se <- cox_pred$se
cox_grid$lo <- cox_pred$fit - 1.96 * cox_pred$se
cox_grid$up <- cox_pred$fit + 1.96 * cox_pred$se

## ------------------------------------------------------------------------
xyplot(cbind(fit, lo, up) ~ rx, data = cox_grid, pch = 20, panel = panel.errbars,
       ylab = "Hazard", xlab = "Treatment", aspect = 3/4)

## ------------------------------------------------------------------------
library(nlme, warn.conflicts = FALSE)
data(Orthodont)

## ------------------------------------------------------------------------
model_lme <- lme(distance ~ Sex * I(age - mean(age, na.rm = TRUE)), random = ~ 1|Subject, 
                 method = "ML", data = Orthodont)
glm_coef(model_lme)

## ------------------------------------------------------------------------
pander(glm_coef(model_lme, labels = c("Constant", "Sex: female-male", "Age (years)", 
                                      "Sex:Age interaction")), split.table = Inf)

## ---- warning=FALSE------------------------------------------------------
plot(Effect(c("age", "Sex"), model_lme, residuals = TRUE), rug = FALSE, xlab = "Age (years)", 
     ylab = "Distance (mm)", main = NULL, aspect = 3/4, partial.residuals = list(pch = 20),
     lines = list(col = c(2, 4), lwd = 1.5))

## ------------------------------------------------------------------------
library(gee, warn.conflicts = FALSE)
model_gee_norm <- gee(distance ~ Sex * I(age - mean(age, na.rm = TRUE)), id = Subject, 
                      data = Orthodont, corstr = "AR-M")

## ------------------------------------------------------------------------
pander(glm_coef(model_gee_norm, labels = c("Constant", "Sex: female-male", "Age (years)", 
                                      "Sex:Age interaction")), split.table = Inf)

## ------------------------------------------------------------------------
Orthodont$fit <- model_gee_norm$fitted.values

## ------------------------------------------------------------------------
xyplot(distance ~ age|Sex, data = Orthodont, groups = Subject, pch = 20, xlab = "Age (years)", 
     ylab = "Distance (mm)", aspect = 3/4) +
  xyplot(fit ~ age|Sex, data = Orthodont, type = "l", lwd = 2)

## ------------------------------------------------------------------------
data(Thall)

c1 <- cbind(Thall[, c(1:5)], count = Thall$y1)[, c(1:4, 6)]
c2 <- cbind(Thall[, c(1:4, 6)], count = Thall$y2)[, c(1:4, 6)]
c3 <- cbind(Thall[, c(1:4, 7)], count = Thall$y3)[, c(1:4, 6)]
c4 <- cbind(Thall[, c(1:4, 8)], count = Thall$y3)[, c(1:4, 6)]
epilepsy <- rbind(c1, c2, c3, c4)

## ----results = "hide"----------------------------------------------------
model_gee <- gee(count ~ treat + base + I(age - mean(age, na.rm = TRUE)), id = factor(id), 
                 data = epilepsy, family = poisson, corstr = "exchangeable", scale.fix = TRUE)

## ------------------------------------------------------------------------
pander(glm_coef(model_gee, labels = c("Constant", "Treatment (Prograbide/Control)", 
                               "Baseline count", "Age (years)")), split.table = Inf)

## ----results = "hide"----------------------------------------------------
library(lme4, warn.conflicts = FALSE)
model_glmer <- glmer(count ~ treat + base + I(age - mean(age, na.rm = TRUE)) + 
                       (1|id), data=epilepsy, family=poisson)

## ------------------------------------------------------------------------
pander(glm_coef(model_glmer, labels = c("Constant", "Treatment (Prograbide/Control)", 
                               "Baseline count", "Age (years)")), split.table = Inf)

## ------------------------------------------------------------------------
plot(Effect(c("age", "treat"), model_glmer), rug = FALSE, lwd = 2, main = NULL,
     xlab = "Age (years)", ylab = "Events", aspect = 3/4, multiline = TRUE)

## ------------------------------------------------------------------------
pander(estat(~ count|treat, data = epilepsy, label = "Number of seizures"))

## ----results = "hide"----------------------------------------------------
model_quasi <- gee(count ~ treat + base + I(age - mean(age, na.rm = TRUE)), id = factor(id), 
                 data = epilepsy, family = quasi(variance = "mu^2", link = "log"), 
                 corstr = "exchangeable")

## ------------------------------------------------------------------------
pander(glm_coef(model_quasi, labels = c("Constant", "Treatment (Prograbide/Control)", 
                               "Baseline count", "Age (years)")), split.table = Inf)

