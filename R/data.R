# Data sets

#' Prevalence of Helicobacter pylori infection in preschool children.
#'
#' A data set containing the prevalence of Helicobacter pylori infection in preschool children
#' according to parental history of duodenal or gastric ulcer.
#' @format A labelled tibble with 863 rows and 2 variables:
#' \describe{
#' \item{ulcer}{History of duodenal or gastric ulcer, factor with levels "No" and "Yes".}
#' \item{infected}{Infected with Helicobacter pylori, factor with levels "No" and "Yes".}
#' }
#' @source Brenner H, Rothenbacher D, Bode G, Adler G (1998) Parental history of gastric or duodenal ulcer
#' and prevalence of Helicobacter pylori infection in preschool children: population based study.
#' BMJ 316:665.
#' @examples
#' data(Brenner)
#'
#' Brenner %>%
#'   cross_tab(infected ~ ulcer)
#'
#' contingency(infected ~ ulcer, data = Brenner, method = "cross.sectional")
"Brenner"

#' Extracorporeal membrane oxygenation in neonates.
#'
#' A clinical trial on the value of extracorporeal membrane oxygenation for term neonates with severe
#' respiratory failure. RCT compares active treatment against conventional management.
#' @format A labelled tibble with 185 rows and 2 variables:
#' \describe{
#' \item{emo}{Extracorporeal membrane oxygenation treatment, factor with levels "No" and "Yes".}
#' \item{survived}{One year survival, factor with levels "No" and "Yes".}
#' }
#' @source Roberts, TE (1998) Extracorporeal Membrane Oxygenation Economics Working Group. Economic
#' evaluation and randomised controlled trial of extracorporeal membrane oxygenation: UK collaborative
#' trial. Brit Med J 317:911-16.
#' @examples
#' data(Roberts)
#'
#' Roberts %>%
#'   cross_tab(survived ~ emo)
"Roberts"

#' Breast cancer and age of childbirth.
#'
#' An international case-control study to test the hypothesis that breast
#' cancer is related to the age that a woman gives childbirth.
#' @format A labelled tibble with 185 rows and 2 variables:
#' \describe{
#' \item{cancer}{Diagnosed with breast cancer, a factor with levels "No" and "Yes".}
#' \item{age}{Age mother gives childbirth, factor with levels "<20", "20-24",
#' "25-29", "30-34" and ">34".}
#' }
#' @source Macmahon, B. et al. (1970). Age at first birth and breast cancer risk. Bull WHO 43, 209-221.
#' @examples
#' data(Macmahon)
#'
#' Macmahon %>%
#'   cross_tab(cancer ~ age)
#'
#' odds_trend(cancer ~ age, data = Macmahon)$df
#' odds_trend(cancer ~ age, data = Macmahon)$fig
"Macmahon"

#' Passive smoking in adulthood and cancer risk.
#'
#' A case-control study to investigate the effects of passive smoking
#' on cancer. Passive smoking was defined as exposure to the cigarette smoke
#' of a spouse who smoked at least one cigarette per day for at least 6 months.
#' @format A labelled tibble with 998 rows and 3 variables:
#' \describe{
#' \item{passive}{Passive smoker, factor with levels "No" and "Yes".}
#' \item{cancer}{Diagnosed with cancer, factor with levels "No" and "Yes".}
#' \item{smoke}{Active smoker, factor with levels "No" and "Yes".}
#' }
#' @source Sandler, DP, Everson, RB, Wilcox, AJ (1985). Passive smoking in adulthood and cancer risk. Amer J Epidem, 121: 37-48.
#' @examples
#' data(Sandler)
#'
#' Sandler %>%
#'   cross_tab(cancer ~ passive)
#'
#' cross_tab(cancer ~ passive + smoke, data = Sandler)
#'
#' mhor(cancer ~ smoke/passive, data = Sandler)
"Sandler"

#' Smoking and mortality in Whickham, England.
#'
#' Data represents women participating in a health survey in Whickham, England in 1972-1974.
#' @format A labelled tibble with 1314 rows and 3 variables:
#' \describe{
#' \item{vstatus}{Vitality status, factor with levels "Alive" and "Death".}
#' \item{smoker}{Smoking status, factor with levels "Non-smoker" and "Smoker".}
#' \item{agegrp}{Age group, factor with levels "18-44", "45-64" and "64+".}
#' }
#' @source Vanderpump, MP, et al (1996) Thyroid, 6:155-160.
#' @source Appleton, DR, French, JM and Vanderpump, PJ (1996) Ignoring a covariate:
#' An example of Simpson's paradox. The American Statistician 50:340-341.
#' @source Vittinghoff, E, Glidden, DV, Shiboski, SC and McCulloh, CE (2005) Regression methods in
#' Biostatistics. Springer.
#' @examples
#' data(Vanderpump)
#'
#' Vanderpump %>%
#'   cross_tab(vstatus ~ .)
#'
#' mhor(vstatus ~ agegrp/smoker, data = Vanderpump)
"Vanderpump"

#' Oral contraceptives and stroke.
#'
#' A case-control study of oral contraceptives and stroke in young women with presence or absence of hypertension.
#' Cases represent thrombotic stroke and controls are hospital controls. The group of no hypertension includes
#' normal blood pressure (<140/90 mm Hg) and borderline hypertension (140-159/90-94 mm Hg). Hypertension group
#' includes moderate hypertension (160-179/95-109 mm Hg) and severe hypertension (180+/110+ mm Hg). This data has
#' been used as an example of join exposure by Rothman for measuring interactions (see examples).
#' @format A labelled tibble with 477 rows and 3 variables:
#' \describe{
#' \item{stroke}{Thrombotic stroke, factor with levels "No" and "Yes".}
#' \item{oc}{Current user of oral contraceptives, factor with levels "Non-user" and "User".}
#' \item{ht}{Hypertension, factor with levels "No" (<160/95 mm Hg) and "Yes".}
#' }
#' @source Collaborative Group for the Study of Stroke in Young Women (1975) Oral contraceptives and stroke in
#' young women. JAMA 231:718-722.
#' @source Rothman, KJ (2002) Epidemiology. An Introduction. Oxford University Press.
#' @examples
#' data(Rothman)
#'
#' cross_tab(stroke ~ oc + ht, data = Rothman)
#'
#' mhor(stroke ~ ht/oc, data = Rothman)
#'
#' ## Model with standard interaction term:
#' model1 <- glm(stroke ~ ht*oc, data = Rothman, family = binomial)
#' glm_coef(model1)
#'
#' ## Model considering join exposure:
#' Rothman$join <- 0
#' Rothman$join[Rothman$oc == "Non-user" & Rothman$ht == "Yes"] <- 1
#' Rothman$join[Rothman$oc == "User" & Rothman$ht == "No"] <- 2
#' Rothman$join[Rothman$oc == "User" & Rothman$ht == "Yes"] <- 3
#' Rothman$join <- factor(Rothman$join, labels=c("Unexposed", "Hypertension", "OC user",
#'                        "OC and hypertension"))
#'
#' require(sjlabelled)
#' Rothman$join <- set_label(Rothman$join, label = "Exposure")
#'
#' cross_tab(stroke ~ join, data = Rothman)
#'
#' model2 <- glm(stroke ~ join, data = Rothman, family = binomial)
#' glm_coef(model2)
#' odds_trend(stroke ~ join, data = Rothman)$df
#' odds_trend(stroke ~ join, data = Rothman)$fig
"Rothman"

#' T-cell counts from Hodgkin's disease patients.
#'
#' Number of CD4+ T-cells and CD8+ T-cells in blood samples from patients in remission from Hodgkin's disease or
#' in remission from disseminated malignancies.
#' @format A labelled tibble with 40 rows and 3 variables:
#' \describe{
#' \item{CD4}{Concentration of CD4+ T-cells (cells / mm^3).}
#' \item{CD8}{Concentration of CD8+ T-cells (cells / mm^3).}
#' \item{Group}{Group, factor with levels "Non-Hodgkin" and "Hodgkin".}
#' }
#' @source Shapiro, CM, et al (1986) Immunologic status of patients in remission from Hodgkin's disease and
#' disseminated malignancies. Am J Med Sci 293:366-370.
#' @source Altman, DA (1991) Practical statistics for medical research. Chapman & Hall/CRC.
#' @examples
#' data(Hodgkin)
#' require(dplyr)
#' require(sjlabelled)
#'
#' Hodgkin <- Hodgkin %>%
#'   mutate(
#'     Ratio = CD4/CD8
#'   ) %>%
#'   var_labels(
#'     Ratio = "CD4+ / CD8+ T-cells"
#'   )
#'
#' estat(~ Ratio|Group, data = Hodgkin)
#'
#' Hodgkin %>%
#'   qq_plot(~ Ratio|Group) %>%
#'   axis_labs()
#'
#' Hodgkin$Ratio <- Hodgkin$CD4/Hodgkin$CD8
#' estat(~ Ratio|Group, data = Hodgkin)
#'
#' qq_plot(~ Ratio|Group, data = Hodgkin) %>%
#' axis_labs()
"Hodgkin"

#' Survival of patients with sepsis.
#'
#' A randomised, double-blind, placebo-controlled trial of intravenous ibuprofen in 455 patients who
#' had sepsis, defined as fever, tachycardia, tachypnea, and acute failure of at least one organ system.
#' @format A labelled tibble with 455 rows and 9 variables:
#' \describe{
#' \item{id}{Patient ID}
#' \item{treat}{Treatment, factor with levels "Placebo" and "Ibuprofen".}
#' \item{race}{Race/ethnicity, factor with levels "White", "African American" and "Other".}
#' \item{fate}{Mortality status at 30 days, factor with levels "Alive" and "Dead".}
#' \item{apache}{Baseline APACHE score.}
#' \item{o2del}{Oxygen delivery at baseline.}
#' \item{followup}{Follow-up time in hours.}
#' \item{temp0}{Baseline temperature in centigrades.}
#' \item{temp10}{Temperature after 36 hr in centigrades.}
#' }
#' @source  Bernard, GR, et al. (1997) The effects of ibuprofen on the physiology and survival of patients
#' with sepsis, N Engl J Med 336: 912–918.
#' @examples
#' data(Bernard)
#'
#' cross_tab(fate ~ treat, data = Bernard)
#'
#' contingency(fate ~ treat, data = Bernard)
"Bernard"

#' Peak knee velocity in walking at flexion and extension.
#'
#' Data of peak knee velocity in walking at flexion and extension in studies about functional performance
#' in cerebral palsy.
#' @format A labelled tibble with 18 rows and 2 variables:
#' \describe{
#' \item{flexion}{Peak knee velocity in gait: flexion (degree/s).}
#' \item{extension}{Peak knee velocity in gait: extension (degree/s).}
#' }
#' @source Tuzson, AE, Granata, KP, and Abel, MF (2003) Spastic velocity threshold constrains functional
#' performance in cerebral palsy. Arch Phys Med Rehabil 84: 1363-1368.
#' @examples
#' data(Tuzson)
#'
#' Tuzson %>%
#'   gf_point(flexion ~ extension) %>%
#'   axis_labs()
#'
#' cor.test(~ flexion + extension, data = Tuzson)
"Tuzson"

#' Measured and self-reported weight in New Zealand.
#'
#' Data on measured and self-reported weight from 40–50 year old participants in the 1989/1990
#' Life In New Zealand Survey.
#' @format A tibble with 343 rows and 4 variables:
#' \describe{
#' \item{srweight}{Self-reported weight in kg.}
#' \item{weight}{Measured weight in kg.}
#' \item{srbmi}{Body mass index calculated from self-reported weight and self-reported height in kg/m^2.}
#' \item{mbmi}{Body mass index calculated from measured weight and measured height in kg/m^2.}
#' }
#' @source Sharples, H, et al. (2012) Agreement between measured and self-reported height, weight and BMI
#' in predominantly European middle-aged New Zealanders: findings from a nationwide 1989 survey. New Zealand
#' Med J 125: 60-69.
#' @examples
#' Sharples %>%
#'   bland_altman(srweight ~ weight, transform = TRUE) %>%
#'   gf_labs(x = "Mean of weights (kg)", y = "Measured weight / Self-reported weight") %>%
#'   gf_theme(theme = sjPlot::theme_sjplot2(base_size = 9))
"Sharples"

#' Migraine pain reduction.
#'
#' Randomised control trial on children suffering from frequent and severe migraine. Control group represents
#' untreated children. The active treatments were either relaxation alone or relaxation with biofeedback.
#' @format A labelled tibble with 18 rows and 2 variables:
#' \describe{
#' \item{pain}{Reduction in weekly headache activity expressed as percentage of baseline data.}
#' \item{group}{Group, a factor with levels "Untreated", "Relaxation" (alone) and "Biofeedback" (relaxation
#' and biofeedback).}
#' }
#' @source Fentress, DW, et al. (1986) Biofeedback and relaxation-response in the treatment of pediatric
#' migraine. Dev Med Child Neurol 28:1 39-46.
#' @source Altman, DA (1991) Practical statistics for medical research. Chapman & Hall/CRC.
#' @examples
#' data(Fentress)
#'
#' Fentress %>%
#'   strip_error(pain ~ group) %>%
#'   axis_labs()
"Fentress"

#' Body weight and plasma volume.
#'
#' Body weight and plasma volume in eight healthy men.
#' @format A labelled data frame with 8 rows and 3 variables:
#' \describe{
#' \item{subject}{Subject ID.}
#' \item{weight}{Body weight in kg.}
#' \item{volume}{Plasma volume in litres.}
#' }
#' @source Kirkwood, BR and Sterne, JAC (2003) Essential Medical Statistics. Second Edition. Blackwell.
#' @examples
#' data(Kirkwood)
#'
#' Kirkwood %>%
#'   gf_point(volume ~ weight) %>%
#'   gf_lm(col = "indianred3", interval = "confidence", fill = "indianred3") %>%
#'   axis_labs()
"Kirkwood"

#' Onchocerciasis in Sierra Leone.
#'
#' Study of onchocerciasis ("river blindness") in Sierra Leone, in which subjects were classified according
#' to whether they lived in villages in savannah or rainforest area.
#' @format A labelled tibble with 1302 rows and 7 variables:
#' \describe{
#' \item{id}{Subject ID.}
#' \item{mf}{Infected with Onchocerciasis volvulus, factor with levels "Not-infected" and "Infected".}
#' \item{area}{Area of residence, factor with levels "Savannah" and "Rainforest".}
#' \item{agegrp}{Age group in years, factor with levels "5-9", "10-19", "20-39" and "40+".}
#' \item{sex}{Subject sex, factor with levels "Male" and "Female".}
#' \item{mfload}{Microfiliariae load.}
#' \item{lesions}{Severe eye lesions, factor with levels "No" and "Yes".}
#' }
#' @source McMahon, JE, Sowa, SIC, Maude, GH and Kirkwood BR (1988) Onchocerciasis in Sierra Leone 2:
#' a comparison of forest and savannah villages. Trans Roy Soc Trop Med Hyg 82: 595-600.
#' @source Kirkwood, BR and Sterne, JAC (2003) Essential Medical Statistics. Second Edition. Blackwell.
#' @examples
#' data(Oncho)
#'
#' odds_trend(mf ~ agegrp, data = Oncho)$df
#' odds_trend(mf ~ agegrp, data = Oncho)$fig
"Oncho"

#' RCT on the treatment of epilepsy.
#'
#' Randomised control trial of an antiepilectic drug (prograbide), in which the number of seizures of 59 patients
#' at baseline and other four follow-up visits were recorded.
#' @format A tibble with 59 rows and 8 variables:
#' \describe{
#' \item{id}{Subject ID.}
#' \item{treat}{Treatment, factor with levels "Control" and "Prograbide".}
#' \item{base}{Number of seizures at baseline.}
#' \item{age}{Age in years at baseline.}
#' \item{y1}{Number of seizures at year one follow-up.}
#' \item{y2}{Number of seizures at year two follow-up.}
#' \item{y3}{Number of seizures at year three follow-up.}
#' \item{y4}{Number of seizures at year four follow-up.}
#' }
#' @source Thall, PF and Vail, SC (1990) Some covariance models for longitudinal count data with over-dispersion.
#' Biometrics, 46: 657-671.
#' @source Stukel, TA (1993) Comparison of methods for the analysis of longitudinal data. Statistics Med 12:
#' 1339-1351.
#' @source Shoukri, MM and Chaudhary, MA (2007) Analysis of correlated data with SAS and R. Third Edition.
#' Chapman & Hall/CRC.
#' @examples
#' data(Thall)
#'
#' c1 <- cbind(Thall[, c(1:5)], count = Thall$y1)[, c(1:4, 6)]
#' c2 <- cbind(Thall[, c(1:4, 6)], count = Thall$y2)[, c(1:4, 6)]
#' c3 <- cbind(Thall[, c(1:4, 7)], count = Thall$y3)[, c(1:4, 6)]
#' c4 <- cbind(Thall[, c(1:4, 8)], count = Thall$y3)[, c(1:4, 6)]
#' epilepsy <- rbind(c1, c2, c3, c4)
#'
#' require(lme4)
#' model_glmer <- glmer(count ~ treat + base + I(age - mean(age, na.rm = TRUE)) +
#'                  (1|id), data = epilepsy, family = poisson)
#' glm_coef(model_glmer, labels = c("Treatment (Prograbide/Control)",
#'                                "Baseline count", "Age (years)"))
"Thall"
