## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  y ~ x, data = my_data

## ----eval=FALSE----------------------------------------------------------
#  y ~ x|z, data = my_data

## ---- message=FALSE------------------------------------------------------
library(pubh, warn.conflicts = FALSE)
library(car, warn.conflicts = FALSE)
library(descr, warn.conflicts = FALSE)
library(multcomp, warn.conflicts = FALSE)
library(pander, warn.conflicts = FALSE)
library(visreg, warn.conflicts = FALSE)

set.alignment("right", row.names = "left", permanent = TRUE)

## ------------------------------------------------------------------------
data(Hodgkin)
Hodgkin$Ratio <- Hodgkin$CD4/Hodgkin$CD8
pander(head(Hodgkin))

## ------------------------------------------------------------------------
estat(~ CD4, data = Hodgkin)

## ------------------------------------------------------------------------
pander(estat(~ CD4, data = Hodgkin, label = "CD4^+^ T cells (cells / mm^3^)"))

## ------------------------------------------------------------------------
pander(estat(~ Ratio|Group, data = Hodgkin, label = "CD4^+^/CD8^+^ T cells"))

## ------------------------------------------------------------------------
pander(estat(Ratio ~ Group, data = Hodgkin, label = "CD4^+^/CD8^+^ T cells"))

## ------------------------------------------------------------------------
pander(var.test(Ratio ~ Group, data = Hodgkin), split.table = Inf)

## ------------------------------------------------------------------------
qq_plot(~ Ratio|Group, data = Hodgkin, ylab = "CD4+ / CD8+ T cells", aspect = 1)

## ------------------------------------------------------------------------
pander(wilcox.test(Ratio ~ Group, data = Hodgkin), split.table = Inf)

## ------------------------------------------------------------------------
strip_error(Ratio ~ Group, data = Hodgkin, ylab = "CD4+ / CD8+ T cells")

## ------------------------------------------------------------------------
strip_error(Ratio ~ Group, data = Hodgkin, ylab = "CD4+ / CD8+ T cells", ylim = c(0, 4.5)) +
  layer(panel.segments(1, 4, 2, 4, lwd=1.5)) +
  layer(panel.text(1.5, 4.1, "**"))

## ------------------------------------------------------------------------
data(birthwt)
birthwt$smoke <- factor(birthwt$smoke, labels = c("Non-smoker", "Smoker"))
birthwt$race <- factor(birthwt$race, labels = c("White", "African American", "Other"))

## ------------------------------------------------------------------------
pander(gen_bst_df(bwt ~ smoke, data = birthwt))

## ------------------------------------------------------------------------
bar_error(bwt ~ smoke, data = birthwt, ylab = "Birth weight (g)", ylim = c(0, 3500))

## ------------------------------------------------------------------------
qq_plot(~ bwt|smoke, data = birthwt, ylab = "Birth weight (g)", aspect = 1)

## ------------------------------------------------------------------------
pander(t.test(bwt ~ smoke, data = birthwt))

## ------------------------------------------------------------------------
bar_error(bwt ~ smoke, data = birthwt, ylab = "Birth weight (g)", ylim = c(0, 3500)) +
  layer(panel.segments(1, 3300, 2, 3300, lwd=1.5)) +
  layer(panel.text(1.5, 3350, "**"))

## ------------------------------------------------------------------------
bar_error(bwt ~ smoke|race, data = birthwt, ylab = "Birth weight (g)", ylim = c(0, 3800),
          col = c("gray80", "gray30"))

## ------------------------------------------------------------------------
strip_error(bwt ~ smoke|race, data = birthwt, ylab = "Birth weight (g)", cex = 0.3)

## ------------------------------------------------------------------------
model_bwt <- lm(bwt ~ smoke + race, data = birthwt)
glm_coef(model_bwt)

## ------------------------------------------------------------------------
unadj <- glm_coef(model_bwt, labels = c("Constant",
                                      "Smoking: smoker - non-smoker",
                                      "Race: African American - White",
                                      "Race: Other - White"))
pander(unadj, split.table = Inf)

## ------------------------------------------------------------------------
model_glht <- glht(model_bwt, linfct  = mcp(race = "Tukey"))
xymultiple(model_glht)

## ------------------------------------------------------------------------
race_glht <- xymultiple(model_glht, plot = FALSE)
pander(race_glht)

## ------------------------------------------------------------------------
final <- unadj
final[, 5] <- as.character(final[, 5])
race_glht[, 5] <- as.character(race_glht[, 5])
final[3:4, 3:5] <- race_glht[1:2, 3:5]

## ------------------------------------------------------------------------
pander(final, split.table = Inf)

## ------------------------------------------------------------------------
visreg(model_bwt, "race", by = "smoke", partial = FALSE, overlay = TRUE, band = FALSE, 
       rug = FALSE, ylab = "Birth weight (g)", xlab = "Race/Ethnicity")

## ------------------------------------------------------------------------
data(Bernard)
pander(head(Bernard), split.table = Inf)

## ------------------------------------------------------------------------
pander(with(Bernard, crosstab(treat, fate, prop.r = TRUE, plot = FALSE)), digits = 2)

## ------------------------------------------------------------------------
dat <- matrix(c(84, 140 , 92, 139), nrow = 2, byrow = TRUE)
epiR::epi.2by2(as.table(dat))

## ------------------------------------------------------------------------
contingency(fate ~ treat, data = Bernard)

## ----eval=FALSE----------------------------------------------------------
#  outcome ~ stratum/exposure, data = my_data

## ------------------------------------------------------------------------
data(oswego, package = "epitools")
oswego$ill <- factor(oswego$ill)
oswego$sex <- factor(oswego$sex)
oswego$chocolate.ice.cream <- factor(oswego$chocolate.ice.cream)
mhor(ill ~ sex/chocolate.ice.cream, data = oswego)

## ------------------------------------------------------------------------
odds_trend(Cancer ~ Age, data = Macmahon)

## ------------------------------------------------------------------------
Freq <- c(1739, 8, 51, 22)
BCG <- gl(2, 1, 4, labels=c("Negative", "Positive"))
Xray <- gl(2, 2, labels=c("Negative", "Positive"))
tb <- data.frame(Freq, BCG, Xray)
tb <- expand_df(tb)

diag_test(BCG ~ Xray, data=tb)

## ------------------------------------------------------------------------
diag_test2(22, 51, 8, 1739)

