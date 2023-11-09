###########################
## figure 1
###########################


###subcortical volume 

cudata_vol$beh <- "cu"
dysregdata_vol$beh <- "dysreg"
extdata_vol$beh <- "ext"


##error bar fig with the model coefficients (in progress)
vol_fig <- bind_rows(cudata_vol,dysregdata_vol,extdata_vol) %>% 
  filter(int == "roi") %>% 
  mutate(signif = ifelse(adjusted_p <= 0.05, "sig", "nonsig"))

vol_fig$signif <- as.factor(vol_fig$signif)

##try to figure out how to change shape fill if signif, right now not working properly##
###going to rename the ROI labels so that the figure is easier to digest
vol_fig$ROI <- sub("smri_vol_scs_", "", vol_fig[,8]) 

##changing the names of the ROIs to make it easier in the figure 
vol_fig$ROI <- c("L_caudate", "R_caudate", "L_putamen", "R_putamen",
                           "L_amygdala", "R_amygdala", "L_accumbens", "R_accumbens", "L_caudate", "R_caudate", "L_putamen", "R_putamen",
                           "L_amygdala", "R_amygdala", "L_accumbens", "R_accumbens", "L_caudate", "R_caudate", "L_putamen", "R_putamen",
                           "L_amygdala", "R_amygdala", "L_accumbens", "R_accumbens")

subcort_vol_LMM <- ggplot(vol_fig, aes(x = ROI, y=Estimate,
        color = beh, shape=signif))+
  geom_point(position = position_dodge(width = 0.5), size = 2.4)+
  scale_fill_manual(values = c("black", "gray"))+
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`),
                position = position_dodge(width = 0.5), 
                linewidth=1)+ ##come back to this something os still off
  geom_hline(yintercept = 0, linetype = "dashed")+
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size=20)) +
  theme(axis.title = element_text(size=20))

###cortical thickness
cudata_thick$beh <- "cu"
dysregdata_thick$beh <- "dysreg"
extdata_thick$beh <- "ext"

##error bar fig with the model coefficients (in progress)
thick_fig <- bind_rows(cudata_thick,dysregdata_thick,extdata_thick) %>% 
  filter(int == "roi") %>% 
  mutate(signif = ifelse(adjusted_p <= 0.05, "sig", "nonsig"))

thick_fig$signif <- as.factor(thick_fig$signif)

##try to figure out how to change shape fill if signif, right now not working properly##
## going to change the names of the ROI labels
thick_fig$ROI <- sub("smri_thick_cdk_", "", thick_fig[,8]) 

##renaming the ROI labels
thick_fig$ROI <- c("L_insula", "R_insula", "L_CaudAnterCingulate", "R_CaudAnterCingulate",
                                 "L_RostAnterCingulate", "R_RostAnterCingulate", "L_ParsOpercularis", "R_ParsOpercularis",
                                 "L_ParsOrbitalis", "R_ParsOrbitalis", "L_ParsTriangularis", "R_ParsTriangularis",
                                 "L_CaudMidFrontal", "R_CaudMidFrontal", "L_RostMidFrontal", "R_RostMidFrontal",
                   "L_insula", "R_insula", "L_CaudAnterCingulate", "R_CaudAnterCingulate",
                   "L_RostAnterCingulate", "R_RostAnterCingulate", "L_ParsOpercularis", "R_ParsOpercularis",
                   "L_ParsOrbitalis", "R_ParsOrbitalis", "L_ParsTriangularis", "R_ParsTriangularis",
                   "L_CaudMidFrontal", "R_CaudMidFrontal", "L_RostMidFrontal", "R_RostMidFrontal",
                   "L_insula", "R_insula", "L_CaudAnterCingulate", "R_CaudAnterCingulate",
                   "L_RostAnterCingulate", "R_RostAnterCingulate", "L_ParsOpercularis", "R_ParsOpercularis",
                   "L_ParsOrbitalis", "R_ParsOrbitalis", "L_ParsTriangularis", "R_ParsTriangularis",
                   "L_CaudMidFrontal", "R_CaudMidFrontal", "L_RostMidFrontal", "R_RostMidFrontal")


cort_thick_LMM <- ggplot(thick_fig, aes(x = ROI, y=Estimate,
                    color=beh, shape=signif))+
  geom_point(position = position_dodge(width = 0.5), size = 2.4)+
  scale_fill_manual(values = c("black", "gray"))+
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`),
                position = position_dodge(width = 0.5),
                linewidth=1)+ ##come back to this something os still off
  geom_hline(yintercept = 0, linetype = "dashed")+
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size=20)) +
  theme(axis.title = element_text(size=20))


####merging the subcortical volume and cortical thickness figures 
fig1_LMM_EX_CT_SV <- gridExtra::grid.arrange(
  grobs = list(subcort_vol_LMM, cort_thick_LMM),
  widths = c(1, 1),
  heights = c(1),
  layout_matrix = rbind(c(1,2)))


###########################
## figure 2
###########################
library(emmeans)
library(sjPlot)
library(modelbased)

cbcl_cdmdfrrh <- lmer(paste0("scale(cbcl_dysreg_raw) ~ scale(b_smri_thick_cdk_cdmdfrrh)*eventname_num + scale(b_cbcl_dysreg_raw) + sex + interview_age_yrs + allmedsYN + household.income + (1|subjectkey) + (1|rel_family_id) + (1|abcd_site)"), data = dtmerged)
sjPlot::tab_model(cbcl_cdmdfrrh)

dtmerged %>% select(CUtotal, b_smri_thick_cdk_rrmdfrrh) %>% summary()

emm_options(pbkrtest.limit = 29749)

emm_options(disable.pbkrtest = TRUE) 

emm_options(emmeans = list(type = "response"),
            contrast = list(infer = c(TRUE, TRUE)))

em <- emmeans(cu_thick_test, ~b_smri_thick_cdk_rrmdfrrh*eventname_num,
              non.nuisance = c("CUTotal", "b_smri_thick_cdk_rrmdfrrh", "eventname_num"),
              at = list(b_smri_thick_cdk_rrmdfrrh = c(2.571,2.640,2.709),
               disable.pbkrtest = TRUE        
              ))







