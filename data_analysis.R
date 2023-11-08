#######################################
### pulling in data sheets 
#######################################

behav_data_complete <- read.csv("/path/to/csv/study3_behav_data_complete.csv")
sMRI_data_baseline <- read.csv("/path/to/csv/study3_sMRI_data_baseline.csv")
demographics_complete <- read.csv("/path/to/csv/study3_demographics_complete.csv")
#family_id <- read.csv("/path/to/csv/abcd_rel_family_id.csv")

##############################
### additional data cleaning
##############################
behav_data_baseline <- behav_data_complete %>% 
  filter(eventname == "baseline_year_1_arm_1")

behav_data_followup <- behav_data_complete %>% 
  filter(eventname %in% c("1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1"))

demographics_baseline <- demographics_complete %>% 
  filter(eventname == "baseline_year_1_arm_1")

demographics_followup <- demographics_complete %>% 
  filter(eventname %in% c("1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1"))

##merge all data
data_list <- list(demographics_complete, behav_data_complete, sMRI_data_baseline)
data_merged <- data_list %>% 
  reduce(full_join, by=c("subjectkey", "eventname"))

##subset baseline and drop NAs (based on key var in each of the 3 frames)
data_complete_baseline <- data_merged %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  drop_na(smri_thick_cdk_banksstslh, abcd_site, interview_age)

data_followup <- left_join(behav_data_followup, demographics_followup,
                           by=c("subjectkey", "eventname"))

##filter out subjects not in baseline 
data_followup_filtered <- semi_join(data_followup, data_complete_baseline,
                                    by = "subjectkey")
##all time points
data_complete <- bind_rows(data_complete_baseline, data_followup_filtered)
  
##add new time variable (factor 0/1/2)
data_complete <- data_complete %>% 
  mutate(eventname_num = recode_factor(eventname, "baseline_year_1_arm_1" = 0, "1_year_follow_up_y_arm_1" = 1, "2_year_follow_up_y_arm_1" = 2))

##add new "all meds" YN variables
data_complete <-
  data_complete %>%
  mutate(allmedssum=rowSums(across(ends_with("YN")), na.rm = TRUE)) %>%
  mutate(allmedsYN=ifelse(allmedssum>0, 1, 0)) 
data_complete$allmedsYN <- as.factor(data_complete$allmedsYN)


data_complete_baseline <-
  data_complete_baseline %>%
  mutate(allmedssum=rowSums(across(ends_with("YN")), na.rm = TRUE)) %>%
  mutate(allmedsYN=ifelse(allmedssum>0, 1, 0))
data_complete_baseline$allmedsYN <- as.factor(data_complete_baseline$allmedsYN)

##other variable transformations
data_complete$sex <- as.factor(data_complete$sex)
data_complete$interview_age_yrs <- data_complete$interview_age /12
data_complete_baseline$sex <- as.factor(data_complete_baseline$sex)
data_complete_baseline$interview_age_yrs <- data_complete_baseline$interview_age /12


data_complete$allmedsYN <- as.factor(data_complete$allmedsYN)
data_complete_baseline$allmedsYN <- as.factor(data_complete_baseline$allmedsYN)

#write.csv(data_complete, "/Users/leepropp/Desktop/Dissertation/ABCD Data/data_complete.csv", row.names = FALSE)
#colSums(is.na(data_followup_filtered))


####################################
## analysis - mixed models 
####################################
#colSums(is.na(baseline_analysis_test))
library(lme4)
library(lmerTest)
library(sjPlot)

#baseline_data_externalizing <- baseline_data_complete %>% 
  #filter(cbcl_scr_syn_external_t >= 60)

##replace the NAs with "no_response" so they don't get dropped when used as control var
##this way it is a factor with 4 levels (lowest income group is the reference group in models)
data_complete_baseline <- data_complete_baseline %>% 
  mutate_at("household.income", ~replace_na(., "[no_response]"))
data_complete_baseline$household.income <- as.factor(data_complete_baseline$household.income)

subcortical_rois_pre <- data_complete_baseline %>% 
  select("smri_vol_scs_caudatelh", "smri_vol_scs_caudaterh", "smri_vol_scs_putamenlh", "smri_vol_scs_putamenrh",
         "smri_vol_scs_amygdalalh", "smri_vol_scs_amygdalarh", "smri_vol_scs_aal", "smri_vol_scs_aar")

###renaming to make figure easier to understand
colnames(subcortical_rois_pre) <- c("L_caudate", "R_caudate", "L_putamen", "R_putamen",
         "L_amygdala", "R_amygdala", "L_accumbens", "R_accumbens")

cortical_rois_pre <- data_complete_baseline %>% 
  select("smri_thick_cdk_insulalh", "smri_thick_cdk_insularh", "smri_thick_cdk_cdacatelh", "smri_thick_cdk_cdacaterh",
         "smri_thick_cdk_rracatelh", "smri_thick_cdk_rracaterh", "smri_thick_cdk_parsopclh", "smri_thick_cdk_parsopcrh",
         "smri_thick_cdk_parsobislh", "smri_thick_cdk_parsobisrh", "smri_thick_cdk_parstgrislh", "smri_thick_cdk_parstgrisrh",
         "smri_thick_cdk_cdmdfrlh", "smri_thick_cdk_cdmdfrrh", "smri_thick_cdk_rrmdfrlh", "smri_thick_cdk_rrmdfrrh")

###renaming to make figure easier to understand
colnames(cortical_rois_pre) <- c("L_insula", "R_insula", "L_CaudAnterCingulate", "R_CaudAnterCingulate",
         "L_RostAnterCingulate", "R_RostAnterCingulate", "L_ParsOpercularis", "R_ParsOpercularis",
         "L_ParsOrbitalis", "R_ParsOrbitalis", "L_ParsTriangularis", "R_ParsTriangularis",
         "L_CaudMidFrontal", "R_CaudMidFrontal", "L_RostMidFrontal", "R_RostMidFrontal")

options(scipen = 999)  


########################################
## FINAL CROSS-SECTIONAL MODELS
########################################
###updates to this set of baseline models are (1) scaling the vars and (2) how it saves the output (saving this way allows for making plots with data and saving in a data frame) (3) added in fdr correction

## baseline CUtotal ~ subcortical volumes
cudata_vol <- data.frame()

for(i in names(subcortical_rois_pre)) {
  cumodels_vol <- c()
  m1 <- lmer(paste0("scale(CUtotal) ~ scale(" ,	i, ")+ sex + interview_age_yrs + allmedsYN + household.income + (1|rel_family_id) + (1|abcd_site)"), data = data_complete_baseline,
             control = lmerControl(optimizer ="bobyqa"))
  cumodels_vol <- c(cumodels_vol, m1)
  
  temp <- bind_cols(summary(m1)$coefficients %>% as.data.frame,
                    confint(m1)[4:11,] %>% as.data.frame)
  
  temp$ROI <- i
  temp$int <- c("i", "roi", "sex", "age", "meds", "ses[>=100K]", "ses[>=50K & <100K]", "ses[no_response]")
  cudata_vol <- bind_rows(cudata_vol, temp[c(2),])
  cudata_vol$adjusted_p <- p.adjust(cudata_vol$`Pr(>|t|)`, method="fdr", n=8)  ##add these to the rest
}

#tab_model(m1)
#qqnorm(residuals(m1))
#sjPlot::plot_residuals(m1)

#res_m1 <- residuals(m1)
#view(res_m1)
#hist(res_m1)

###baseline dysreg ~ subcortical volumes
dysregdata_vol <- data.frame()

for(i in names(subcortical_rois_pre)) {
  dysregmodels_vol <- c()
  m2 <- lmer(paste0("scale(cbcl_dysreg_raw) ~ scale(",	i, ")+ sex + interview_age_yrs + allmedsYN + household.income + (1|rel_family_id) + (1|abcd_site)"), data = data_complete_baseline,
             control = lmerControl(optimizer ="bobyqa"))
  dysregmodels_vol <- c(dysregmodels_vol, m2)
  
  temp <- bind_cols(summary(m2)$coefficients %>% as.data.frame,
                    confint(m2)[4:11,] %>% as.data.frame)
  
  temp$ROI <- i
  temp$int <- c("i", "roi", "sex", "age", "meds", "ses[>=100K]", "ses[>=50K & <100K]", "ses[no_response]")
  dysregdata_vol <- bind_rows(dysregdata_vol, temp[c(2),])
  dysregdata_vol$adjusted_p <- p.adjust(dysregdata_vol$`Pr(>|t|)`, method="fdr", n=8) 
}

##baseline externalizing ~ subcortical volumes
extdata_vol <- data.frame()

for(i in names(subcortical_rois_pre)) {
  extmodels_vol <- c()
  m3 <- lmer(paste0("scale(ksads_ext_symptoms_total) ~ scale(",	i, ")+ sex + interview_age_yrs + allmedsYN + household.income + (1|rel_family_id) + (1|abcd_site)"), data = data_complete_baseline,
             control = lmerControl(optimizer ="bobyqa"))
  extmodels_vol <- c(extmodels_vol, m3)
  
  temp <- bind_cols(summary(m3)$coefficients %>% as.data.frame,
                    confint(m3)[4:11,] %>% as.data.frame)
  
  temp$ROI <- i
  temp$int <- c("i", "roi", "sex", "age", "meds", "ses[>=100K]", "ses[>=50K & <100K]", "ses[no_response]")
  extdata_vol <- bind_rows(extdata_vol, temp[c(2),])
  extdata_vol$adjusted_p <- p.adjust(extdata_vol$`Pr(>|t|)`, method="fdr", n=8)
}


## baseline CUtotal ~ cortical thickness
cudata_thick <- data.frame()

for(i in names(cortical_rois_pre)) {
  cumodels_thick <- c()
  m4 <- lmer(paste0("scale(CUtotal) ~ scale(" ,	i, ")+ sex + interview_age_yrs + allmedsYN + household.income + (1|rel_family_id) + (1|abcd_site)"), data = data_complete_baseline,
             control = lmerControl(optimizer ="bobyqa"))
  cumodels_thick <- c(cumodels_thick, m4)
  
  temp <- bind_cols(summary(m4)$coefficients %>% as.data.frame,
                    confint(m4)[4:11,] %>% as.data.frame)
  
  temp$ROI <- i
  temp$int <- c("i", "roi", "sex", "age", "meds", "ses[>=100K]", "ses[>=50K & <100K]", "ses[no_response]")
  cudata_thick <- bind_rows(cudata_thick, temp[c(2),])
  cudata_thick$adjusted_p <- p.adjust(cudata_thick$`Pr(>|t|)`, method="fdr", n=16) ##is n=18 correct here?
}

###baseline dysreg ~ cortical thickness
dysregdata_thick <- data.frame()

for(i in names(cortical_rois_pre)) {
  dysregmodels_thick <- c()
  m5 <- lmer(paste0("scale(cbcl_dysreg_raw) ~ scale(",	i, ")+ sex + interview_age_yrs + allmedsYN + household.income + (1|rel_family_id) + (1|abcd_site)"), data = data_complete_baseline,
             control = lmerControl(optimizer ="bobyqa"))
  dysregmodels_thick <- c(dysregmodels_thick, m5)
  
  temp <- bind_cols(summary(m5)$coefficients %>% as.data.frame,
                    confint(m5)[4:11,] %>% as.data.frame)
  
  temp$ROI <- i
  temp$int <- c("i", "roi", "sex", "age", "meds", "ses[>=100K]", "ses[>=50K & <100K]", "ses[no_response]")
  dysregdata_thick <- bind_rows(dysregdata_thick, temp[c(2),])
  dysregdata_thick$adjusted_p <- p.adjust(dysregdata_thick$`Pr(>|t|)`, method="fdr", n=16)
}

##baseline externalizing ~ cortical thickness
extdata_thick <- data.frame()

for(i in names(cortical_rois_pre)) {
  extmodels_thick <- c()
  m6 <- lmer(paste0("scale(ksads_ext_symptoms_total) ~ scale(",	i, ")+ sex + interview_age_yrs + allmedsYN + household.income + (1|rel_family_id) + (1|abcd_site)"), data = data_complete_baseline,
             control = lmerControl(optimizer ="bobyqa"))
  extmodels_thick <- c(extmodels_thick, m6)
  
  temp <- bind_cols(summary(m6)$coefficients %>% as.data.frame,
                    confint(m6)[4:11,] %>% as.data.frame)
  
  temp$ROI <- i
  temp$int <- c("i", "roi", "sex", "age", "meds", "ses[>=100K]", "ses[>=50K & <100K]", "ses[no_response]")
  extdata_thick <- bind_rows(extdata_thick, temp[c(2),])
  extdata_thick$adjusted_p <- p.adjust(extdata_thick$`Pr(>|t|)`, method="fdr", n=16)
}


###################################################
## FINAL LONGITUDINAL MIXED EFFECTS MODELS
##################################################

##set up for longitudunal mixed models 

#data_complete$abcd_site <- as.factor(data_complete$abcd_site)
#data_complete$rel_family_id <- as.character(data_complete$rel_family_id)
#data_complete$subjectkey <- as.character(data_complete$subjectkey)

##same as did for baseline  
data_complete <- data_complete %>% 
  mutate_at("household.income", ~replace_na(., "[no_response]")) ##run this before creating dtmerged!
data_complete$household.income <- as.factor(data_complete$household.income)

dt <- data_complete 
dtpre <- data_complete %>% filter(eventname_num == 0) %>% select(subjectkey, starts_with("smri_thick_cdk_"), starts_with("smri_vol_scs"), 
                                                                 household.income, cbcl_dysreg_raw, CUtotal, ksads_ext_symptoms_total) %>% 
  rename_with(~ glue("b_{.}"), starts_with("smri_thick_cdk_") | starts_with("smri_vol_scs") | household.income | cbcl_dysreg_raw | CUtotal | ksads_ext_symptoms_total) 
 
dtmerged <- left_join(dt, dtpre) %>% 
  relocate(323:326, .before = high.educ) 


#test <- lmer(CUtotal ~ b_smri_vol_scs_aar * eventname_num + b_CUtotal + sex + interview_age_yrs + allmedsYN + (1|subjectkey) + (1|rel_family_id) + (1|abcd_site), data = dtmerged)
#summary(test) 

subcortical_rois <- dtmerged %>% 
  select("b_smri_vol_scs_caudatelh", "b_smri_vol_scs_caudaterh", "b_smri_vol_scs_putamenlh", "b_smri_vol_scs_putamenrh",
         "b_smri_vol_scs_amygdalalh", "b_smri_vol_scs_amygdalarh", "b_smri_vol_scs_aal", "b_smri_vol_scs_aar") #

cortical_rois <- dtmerged %>% 
  select("b_smri_thick_cdk_insulalh", "b_smri_thick_cdk_insulalh", "b_smri_thick_cdk_cdacatelh", "b_smri_thick_cdk_cdacaterh",
         "b_smri_thick_cdk_rracatelh", "b_smri_thick_cdk_rracaterh", "b_smri_thick_cdk_parsopclh", "b_smri_thick_cdk_parsopcrh",
         "b_smri_thick_cdk_parsobislh", "b_smri_thick_cdk_parsobisrh", "b_smri_thick_cdk_parstgrislh", "b_smri_thick_cdk_parstgrisrh",
         "b_smri_thick_cdk_cdmdfrlh", "b_smri_thick_cdk_cdmdfrrh", "b_smri_thick_cdk_rrmdfrlh", "b_smri_thick_cdk_rrmdfrrh",
         "b_smri_thick_cdk_rrmdfrlh", "b_smri_thick_cdk_rrmdfrrh")

##CU trajectory ~ subcortical volumes*time
cumod_vol_summaries_long <- list()
for(i in names(subcortical_rois)) {
  cumod_vol_summaries_long[[i]] <- summary(
    lmer(paste0("scale(CUtotal) ~ scale(",	i, ")*eventname_num + scale(b_CUtotal) + sex + interview_age_yrs + allmedsYN + household.income + (1|subjectkey) + (1|rel_family_id) + (1|abcd_site)"), data = dtmerged,
         control = lmerControl(optimizer ="bobyqa")))
}
cumod_vol_summaries_long


cumod_vol_summaries_long_anova <- list()
for(i in names(subcortical_rois)) {
  cumod_vol_summaries_long_anova[[i]] <- anova(
    lmer(paste0("scale(CUtotal) ~ scale(",	i, ")*eventname_num + scale(b_CUtotal) + sex + interview_age_yrs + allmedsYN + household.income + (1|subjectkey) + (1|rel_family_id) + (1|abcd_site)"), data = dtmerged,
         control = lmerControl(optimizer ="bobyqa")))
}
cumod_vol_summaries_long_anova

cumod_vol_summaries_long_anova_coef <- cbind.data.frame(cumod_vol_summaries_long_anova$b_smri_vol_scs_caudatelh$`Pr(>F)`, 
                                      cumod_vol_summaries_long_anova$b_smri_vol_scs_caudaterh$`Pr(>F)`,
                                      cumod_vol_summaries_long_anova$b_smri_vol_scs_putamenlh$`Pr(>F)`,
                                      cumod_vol_summaries_long_anova$b_smri_vol_scs_putamenrh$`Pr(>F)`,
                                      cumod_vol_summaries_long_anova$b_smri_vol_scs_amygdalalh$`Pr(>F)`,
                                      cumod_vol_summaries_long_anova$b_smri_vol_scs_amygdalarh$`Pr(>F)`, 
                                      cumod_vol_summaries_long_anova$b_smri_vol_scs_aal$`Pr(>F)`, 
                                      cumod_vol_summaries_long_anova$b_smri_vol_scs_aar$`Pr(>F)`)

cumod_vol_summaries_long_anova_coef_pval <- as.data.frame(t(cumod_vol_summaries_long_anova_coef[c(8),]))
colnames(cumod_vol_summaries_long_anova_coef_pval) <- "p-values"
cumod_vol_summaries_long_anova_coef_pval$adjusted_p <- p.adjust(cumod_vol_summaries_long_anova_coef_pval$`p-values`, method="fdr", n=8)


##dysreg trajectory ~ subcortical volumes*time
dysregmod_vol_summaries_long <- list()
for(i in names(subcortical_rois)) {
  dysregmod_vol_summaries_long[[i]] <- summary(
    lmer(paste0("scale(cbcl_dysreg_raw) ~ scale(",	i, ")*eventname_num + scale(b_cbcl_dysreg_raw) + sex + interview_age_yrs + allmedsYN + household.income + (1|subjectkey) + (1|rel_family_id) + (1|abcd_site)"), data = dtmerged,
         control = lmerControl(optimizer ="bobyqa"))
  )
}
dysregmod_vol_summaries_long
    
##ext trajectory ~ subcortical volumes*time
extmod_vol_summaries_long_results <- data.frame()
extmod_vol_summaries_long_results <- lapply(names(subcortical_rois), function(i) {
  lmer(paste0("scale(ksads_ext_symptoms_total) ~ scale(",	i, ")*eventname_num + scale(b_ksads_ext_symptoms_total) + sex + interview_age_yrs + allmedsYN + household.income + (1|subjectkey) + (1|rel_family_id) + (1|abcd_site)"), data = dtmerged,
       control = lmerControl(optimizer ="bobyqa"))
})


for(i in names(subcortical_rois)) {
  extmod_vol_summaries_long[[i]] <- summary(
    lmer(paste0("scale(ksads_ext_symptoms_total) ~ scale(",	i, ")*eventname_num + scale(b_ksads_ext_symptoms_total) + sex + interview_age_yrs + allmedsYN + household.income + (1|subjectkey) + (1|rel_family_id) + (1|abcd_site)"), data = dtmerged,
         control = lmerControl(optimizer ="bobyqa"))
  )
}
extmod_vol_summaries_long


extmod_df <- list()
extmod_df <- for(i in names(subcortical_rois)) {
  temp <- bind_cols(summary(extmod_vol_summaries_long)$coefficients %>% as.data.frame)
  extmod_df <- extmod_vol_summaries_long[[i]]$coefficients
}

extmod_df <- bind_rows(extmod_df, temp[c(2),])
cumod_vol_summaries_long_results$adjusted_p <- p.adjust(cumod_vol_summaries_long_results$`Pr(>|t|)`, method="fdr")



##CU trajectory ~ cortical thickness*time
cumod_thick_summaries_long <- list()
for(i in names(cortical_rois)) {
  cumod_thick_summaries_long[[i]] <- summary(
    lmer(paste0("scale(CUtotal) ~ scale(",	i, ")*eventname_num + scale(b_CUtotal) + sex + interview_age_yrs + allmedsYN + household.income + (1|subjectkey) + (1|rel_family_id) + (1|abcd_site)"), data = dtmerged,
         control = lmerControl(optimizer ="bobyqa"))
  )
}
cumod_thick_summaries_long



ksads_thick_test1 <- lmer(paste0("scale(ksads_ext_symptoms_total) ~ scale(b_smri_thick_cdk_rrmdfrlh)*eventname_num + scale(b_ksads_ext_symptoms_total) + sex + interview_age_yrs + allmedsYN + household.income + (1|subjectkey) + (1|rel_family_id) + (1|abcd_site)"), data = dtmerged, control = lmerControl(optimizer ="bobyqa"))
ksads_thick_test2 <- lmer(paste0("scale(ksads_ext_symptoms_total) ~ scale(b_smri_thick_cdk_parstgrislh)*eventname_num + scale(b_ksads_ext_symptoms_total) + sex + interview_age_yrs + allmedsYN + household.income + (1|subjectkey) + (1|rel_family_id) + (1|abcd_site)"), data = dtmerged, control = lmerControl(optimizer ="bobyqa"))

anova(ksads_thick_test1)
anova(ksads_thick_test2)

summary(anova(cu_thick_test))

##dysreg trajectory ~ subcortical volumes*time
dysregmod_thick_summaries_long <- list()
for(i in names(cortical_rois)) {
  dysregmod_thick_summaries_long[[i]] <- summary(
    lmer(paste0("scale(cbcl_dysreg_raw) ~ scale(",	i, ")*eventname_num + scale(b_cbcl_dysreg_raw) + sex + interview_age_yrs + allmedsYN + household.income + (1|subjectkey) + (1|rel_family_id) + (1|abcd_site)"), data = dtmerged,
         control = lmerControl(optimizer ="bobyqa"))
  )
}
dysregmod_thick_summaries_long

##ext trajectory ~ subcortical volumes*time
extmod_thick_summaries_long <- list()
for(i in names(cortical_rois)) {
  extmod_thick_summaries_long[[i]] <- summary(
    lmer(paste0("scale(ksads_ext_symptoms_total) ~ scale(",	i, ")*eventname_num + scale(b_ksads_ext_symptoms_total) + sex + interview_age_yrs + allmedsYN + household.income + (1|subjectkey) + (1|rel_family_id) + (1|abcd_site)"), data = dtmerged,
         control = lmerControl(optimizer ="bobyqa"))
  )
}
anova(extmod_thick_summaries_long$b_smri_thick_cdk_parstgrislh)
