###loading any necessary packages 
library(ggplot2)
library(plotly)
library(sparklyr)
library(magrittr)
library(lubridate)
library(dplyr)
library(tidyverse)
library(psych)
library(stringr)
library(reshape2)

#######################################
### pulling in data sheets 
#######################################

###behavioural variables
CBCL_data_scales <- read.csv("/path/to/csv/abcd_cbcls01.csv")
CBCL_data_raw <- read.csv("/path/to/csv/abcd_cbcl_raw01.csv")
SDQ_prosocial <- read.csv("/path/to/csv/psb_sdq01.csv")
ksads_full <- read.csv("/path/to/csv/abcd_ksad01.csv")

###imaging variables
sMRI_data <- read.csv("/path/to/csv/abcd_smrip10201.csv")
MRI_info <- read.csv("/path/to/csv/abcd_mri01.csv")
imaging_info <- read.csv("/path/to/csv/abcd_imgincl01.csv")
exclusion_criteria <- read.csv("/path/to/csv/exclusion_criteria.csv")

###demographic and admin variables
family_relationships <- read.csv("/path/to/csv/family_relationships.csv")
abcd_site <- read.csv("/path/to/csv/abcd_site.csv")
meds_data <- read.csv("/path/to/csv/abcd_medsy01.csv")
demographics_survey <- read.csv("/path/to/csv/abcd_pdem02.csv")
#family_id <- read.csv("/path/to/csv/abcd_rel_family_id.csv")

######################################
### cleaning sMRI data
######################################

##filter for baseline only and select variables of interest (exclusion)
sMRI_exclusion_criteria_baseline <- exclusion_criteria %>% 
  filter(event_name == "baseline_year_1_arm_1") %>%
  select("src_subject_id", "event_name", "mrif_score", "fsqc_qc", "iqc_t1_ok_ser")

##rename to match id and events variables 
names(sMRI_exclusion_criteria_baseline)[names(sMRI_exclusion_criteria_baseline) == "src_subject_id"] <- "subjectkey"
names(sMRI_exclusion_criteria_baseline)[names(sMRI_exclusion_criteria_baseline) == "event_name"] <- "eventname"

##inclusion variable (baseline only)
sMRI_dairc_inclusion_baseline <- imaging_info %>% 
  filter(eventname == "baseline_year_1_arm_1") %>%
  select("subjectkey", "eventname", "imgincl_t1w_include")

##manufacturer info (baseline only)
MRI_info_manufacturer_baseline <- MRI_info %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select("subjectkey", "eventname", "mri_info_manufacturersmn")
  
##cortical thickness and subcortical volumes (baseline only)
sMRI_data_CorticalThickness_subVols_baseline <- sMRI_data %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  select("subjectkey", "eventname", starts_with("smri_thick_cdk_"), starts_with("smri_vol_scs"), "smri_vol_cdk_total")

##merge all imaging exclusion criteria (this drops n=94 subjects who are missing dairc inclusion criteria var) 
sMRI_exclusion_inclusion_baseline <- inner_join(sMRI_exclusion_criteria_baseline, sMRI_dairc_inclusion_baseline, 
                                              by=c("subjectkey", "eventname"))

##merge exclusion criteria with manufacturer info (for all subjects who have full exclusion criteria)
sMRI_exclusion_manufacturer_baseline <- left_join(sMRI_exclusion_inclusion_baseline, MRI_info_manufacturer_baseline, 
                                                  by=c("subjectkey", "eventname"))

##merge all sMRI variables and drop any subjects with NAs in any var except manufacturer info (in the end that NA overlapped anyways)
sMRI_all_vars_baseline <- left_join(sMRI_data_CorticalThickness_subVols_baseline, sMRI_exclusion_manufacturer_baseline,
                                    by=c("subjectkey", "eventname")) %>% 
  drop_na(-mri_info_manufacturersmn)

## (1) removing participants with incidental findings (2) removing participants with poor T1 
## (3) removing participants who failed freesurfer qc (4) removing participants who failed dairc qc
sMRI_complete_qc_baseline <- sMRI_all_vars_baseline %>% 
  filter(!mrif_score %in% c("Consider clinical referral", "Consider immediate clinical referral")) %>%
  filter(iqc_t1_ok_ser != 0) %>%
  filter(fsqc_qc != "reject") %>%
  filter(imgincl_t1w_include == 1) 
  
##colSums(is.na(sMRI_complete_qc_baseline))
###lower n at end than Hajer - check about this


#####################################
### demographics and admin
#####################################

demographics_all <- exclusion_criteria %>% 
  filter(event_name %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1")) %>%
  select("src_subject_id", "event_name", "female", "household.income", "high.educ", 
         "race_ethnicity", "race.6level","sex_at_birth")

##rename to match id and events variables 
names(demographics_all)[names(demographics_all) == "src_subject_id"] <- "subjectkey"
names(demographics_all)[names(demographics_all) == "event_name"] <- "eventname"

family_data_baseline <- family_relationships %>% ##there is only baseline data for these vars (take family id from the site datasheet bc it is complete)
  filter(eventname == "baseline_year_1_arm_1") %>%
  select("subjectkey", "eventname", "rel_relationship") 

abcd_site_all <- abcd_site %>%
  filter(event_name %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1")) %>%
  select("src_subject_id", "event_name", "abcd_site", "rel_family_id")

##rename to match id and events variables 
names(abcd_site_all)[names(abcd_site_all) == "src_subject_id"] <- "subjectkey"
names(abcd_site_all)[names(abcd_site_all) == "event_name"] <- "eventname"

##meds data (pulled separately so merge on it's own)
meds_all <- meds_data %>% 
  filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1")) %>%
  select(subjectkey,
         eventname,
         med1_rxnorm_p,
         med2_rxnorm_p,
         med3_rxnorm_p,
         med4_rxnorm_p,
         med5_rxnorm_p,
         med6_rxnorm_p,
         med7_rxnorm_p,
         med8_rxnorm_p,
         med9_rxnorm_p,
         med10_rxnorm_p,
         med11_rxnorm_p,
         med12_rxnorm_p,
         med13_rxnorm_p,
         med14_rxnorm_p,
         med15_rxnorm_p) %>% 
  mutate_all(na_if, "")

# Define search function for stimulant medications (as per list of meds from Shoval et al, 2021)        
searchfunction_stimulants <- function(column) {
  (str_detect(column, regex("Methylphenidate", ignore_case=TRUE))
   | str_detect(column, regex("Ritalin", ignore_case=TRUE))
   | str_detect(column, regex("Quillivant", ignore_case=TRUE))
   | str_detect(column, regex("Quillichew", ignore_case=TRUE))
   | str_detect(column, regex("Methylin", ignore_case=TRUE))
   | str_detect(column, regex("Concerta", ignore_case=TRUE))
   | str_detect(column, regex("Metadate", ignore_case=TRUE))
   | str_detect(column, regex("Aptensio", ignore_case=TRUE))
   | str_detect(column, regex("Daytrana", ignore_case=TRUE))
   | str_detect(column, regex("Dexmethylphenidate", ignore_case=TRUE))
   | str_detect(column, regex("Focalin", ignore_case=TRUE))
   | str_detect(column, regex("Dextroamphetamine", ignore_case=TRUE))
   | str_detect(column, regex("Amphetamine", ignore_case=TRUE))
   | str_detect(column, regex("Adderall", ignore_case=TRUE))
   | str_detect(column, regex("Adzenys", ignore_case=TRUE))
   | str_detect(column, regex("Dynavel", ignore_case=TRUE))
   | str_detect(column, regex("Evekeo", ignore_case=TRUE))
   | str_detect(column, regex("Zenzedi", ignore_case=TRUE))
   | str_detect(column, regex("Dexedrine ", ignore_case=TRUE))
   | str_detect(column, regex("Procentra", ignore_case=TRUE))
   | str_detect(column, regex("Lisdexamfetamine", ignore_case=TRUE))
   | str_detect(column, regex("Vyvanse", ignore_case=TRUE)))
}

# Create med1 to med15 columns and include a 1 where any stimulant medication from above is listed
meds_all <-
  mutate(meds_all, stim1=as.numeric(searchfunction_stimulants(med1_rxnorm_p)))%>%
  mutate(stim2=as.numeric(searchfunction_stimulants(med2_rxnorm_p))) %>%
  mutate(stim3=as.numeric(searchfunction_stimulants(med3_rxnorm_p))) %>%
  mutate(stim4=as.numeric(searchfunction_stimulants(med4_rxnorm_p))) %>%
  mutate(stim5=as.numeric(searchfunction_stimulants(med5_rxnorm_p))) %>%
  mutate(stim6=as.numeric(searchfunction_stimulants(med6_rxnorm_p))) %>%
  mutate(stim7=as.numeric(searchfunction_stimulants(med7_rxnorm_p))) %>%
  mutate(stim8=as.numeric(searchfunction_stimulants(med8_rxnorm_p))) %>%
  mutate(stim9=as.numeric(searchfunction_stimulants(med9_rxnorm_p))) %>%
  mutate(stim10=as.numeric(searchfunction_stimulants(med10_rxnorm_p))) %>%
  mutate(stim11=as.numeric(searchfunction_stimulants(med11_rxnorm_p))) %>%
  mutate(stim12=as.numeric(searchfunction_stimulants(med12_rxnorm_p))) %>%
  mutate(stim13=as.numeric(searchfunction_stimulants(med13_rxnorm_p))) %>%
  mutate(stim14=as.numeric(searchfunction_stimulants(med14_rxnorm_p))) %>%
  mutate(stim15=as.numeric(searchfunction_stimulants(med15_rxnorm_p)))

# Create number of stimulants column and stimulants yes / no column, taking account of missing variables
meds_all <-
  meds_all %>%
  mutate(stimsum=rowSums(across(stim1:stim15), na.rm = TRUE)) %>%
  mutate(stimYN=ifelse(stimsum>0, 1, 0)) 

# Check on numbers
table(meds_all$stimsum)
table(meds_all$stimYN)

##search function for alpha agonists (all rest repeat of above)
searchfunction_alphaagonist <- function(column) {
(str_detect(column, regex("Guanfacine", ignore_case=TRUE))
| str_detect(column, regex("Intuniv", ignore_case=TRUE))
| str_detect(column, regex("Tenex", ignore_case=TRUE))
| str_detect(column, regex("Clonidine", ignore_case=TRUE))
| str_detect(column, regex("Catapres", ignore_case=TRUE))
| str_detect(column, regex("Nexiclon", ignore_case=TRUE)))
}

meds_all <-
  mutate(meds_all, aa1=as.numeric(searchfunction_alphaagonist(med1_rxnorm_p)))%>%
  mutate(aa2=as.numeric(searchfunction_alphaagonist(med2_rxnorm_p))) %>%
  mutate(aa3=as.numeric(searchfunction_alphaagonist(med3_rxnorm_p))) %>%
  mutate(aa4=as.numeric(searchfunction_alphaagonist(med4_rxnorm_p))) %>%
  mutate(aa5=as.numeric(searchfunction_alphaagonist(med5_rxnorm_p))) %>%
  mutate(aa6=as.numeric(searchfunction_alphaagonist(med6_rxnorm_p))) %>%
  mutate(aa7=as.numeric(searchfunction_alphaagonist(med7_rxnorm_p))) %>%
  mutate(aa8=as.numeric(searchfunction_alphaagonist(med8_rxnorm_p))) %>%
  mutate(aa9=as.numeric(searchfunction_alphaagonist(med9_rxnorm_p))) %>%
  mutate(aa10=as.numeric(searchfunction_alphaagonist(med10_rxnorm_p))) %>%
  mutate(aa11=as.numeric(searchfunction_alphaagonist(med11_rxnorm_p))) %>%
  mutate(aa12=as.numeric(searchfunction_alphaagonist(med12_rxnorm_p))) %>%
  mutate(aa13=as.numeric(searchfunction_alphaagonist(med13_rxnorm_p))) %>%
  mutate(aa14=as.numeric(searchfunction_alphaagonist(med14_rxnorm_p))) %>%
  mutate(aa15=as.numeric(searchfunction_alphaagonist(med15_rxnorm_p)))

meds_all <-
  meds_all %>%
  mutate(aasum=rowSums(across(aa1:aa15), na.rm = TRUE)) %>%
  mutate(aaYN=ifelse(aasum>0, 1, 0)) 

# Check on numbers
table(meds_all$aasum)
table(meds_all$aaYN)

##same but for atomoxitine
searchfunction_atomox <- function(column) {
  (str_detect(column, regex("Atomoxetine", ignore_case=TRUE))
   | str_detect(column, regex("Strattera", ignore_case=TRUE)))
}

meds_all <-
  mutate(meds_all, atomox1=as.numeric(searchfunction_atomox(med1_rxnorm_p)))%>%
  mutate(atomox2=as.numeric(searchfunction_atomox(med2_rxnorm_p))) %>%
  mutate(atomox3=as.numeric(searchfunction_atomox(med3_rxnorm_p))) %>%
  mutate(atomox4=as.numeric(searchfunction_atomox(med4_rxnorm_p))) %>%
  mutate(atomox5=as.numeric(searchfunction_atomox(med5_rxnorm_p))) %>%
  mutate(atomox6=as.numeric(searchfunction_atomox(med6_rxnorm_p))) %>%
  mutate(atomox7=as.numeric(searchfunction_atomox(med7_rxnorm_p))) %>%
  mutate(atomox8=as.numeric(searchfunction_atomox(med8_rxnorm_p))) %>%
  mutate(atomox9=as.numeric(searchfunction_atomox(med9_rxnorm_p))) %>%
  mutate(atomox10=as.numeric(searchfunction_atomox(med10_rxnorm_p))) %>%
  mutate(atomox11=as.numeric(searchfunction_atomox(med11_rxnorm_p))) %>%
  mutate(atomox12=as.numeric(searchfunction_atomox(med12_rxnorm_p))) %>%
  mutate(atomox13=as.numeric(searchfunction_atomox(med13_rxnorm_p))) %>%
  mutate(atomox14=as.numeric(searchfunction_atomox(med14_rxnorm_p))) %>%
  mutate(atomox15=as.numeric(searchfunction_atomox(med15_rxnorm_p)))

meds_all <-
  meds_all %>%
  mutate(atomoxsum=rowSums(across(atomox1:atomox15), na.rm = TRUE)) %>%
  mutate(atomoxYN=ifelse(atomoxsum>0, 1, 0)) 

# Check on numbers
table(meds_all$atomoxsum)
table(meds_all$atomoxYN)

##same but for antipsychotics
searchfunction_antipsych <- function(column) {
  (str_detect(column, regex("Aripiprazole", ignore_case=TRUE))
   | str_detect(column, regex("Abilify", ignore_case=TRUE))
   | str_detect(column, regex("Risperidone", ignore_case=TRUE))
   | str_detect(column, regex("Risperdal", ignore_case=TRUE))
   | str_detect(column, regex("Lurasidone", ignore_case=TRUE))
   | str_detect(column, regex("Latuda", ignore_case=TRUE))
   | str_detect(column, regex("Ziprasidone", ignore_case=TRUE))
   | str_detect(column, regex("Clozapine", ignore_case=TRUE))
   | str_detect(column, regex("Quetiapine", ignore_case=TRUE))
   | str_detect(column, regex("Seroquel", ignore_case=TRUE)))
}

meds_all <-
  mutate(meds_all, antipsych1=as.numeric(searchfunction_antipsych(med1_rxnorm_p)))%>%
  mutate(antipsych2=as.numeric(searchfunction_antipsych(med2_rxnorm_p))) %>%
  mutate(antipsych3=as.numeric(searchfunction_antipsych(med3_rxnorm_p))) %>%
  mutate(antipsych4=as.numeric(searchfunction_antipsych(med4_rxnorm_p))) %>%
  mutate(antipsych5=as.numeric(searchfunction_antipsych(med5_rxnorm_p))) %>%
  mutate(antipsych6=as.numeric(searchfunction_antipsych(med6_rxnorm_p))) %>%
  mutate(antipsych7=as.numeric(searchfunction_antipsych(med7_rxnorm_p))) %>%
  mutate(antipsych8=as.numeric(searchfunction_antipsych(med8_rxnorm_p))) %>%
  mutate(antipsych9=as.numeric(searchfunction_antipsych(med9_rxnorm_p))) %>%
  mutate(antipsych10=as.numeric(searchfunction_antipsych(med10_rxnorm_p))) %>%
  mutate(antipsych11=as.numeric(searchfunction_antipsych(med11_rxnorm_p))) %>%
  mutate(antipsych12=as.numeric(searchfunction_antipsych(med12_rxnorm_p))) %>%
  mutate(antipsych13=as.numeric(searchfunction_antipsych(med13_rxnorm_p))) %>%
  mutate(antipsych14=as.numeric(searchfunction_antipsych(med14_rxnorm_p))) %>%
  mutate(antipsych15=as.numeric(searchfunction_antipsych(med15_rxnorm_p)))

meds_all <-
  meds_all %>%
  mutate(antipsychsum=rowSums(across(antipsych1:antipsych15), na.rm = TRUE)) %>%
  mutate(antipsychYN=ifelse(antipsychsum>0, 1, 0)) 

# Check on numbers
table(meds_all$antipsychsum)
table(meds_all$antipsychYN)

##same but for anti depressants
searchfunction_antidep <- function(column) {
  (str_detect(column, regex("Fluoxetine", ignore_case=TRUE))
   | str_detect(column, regex("Prozac", ignore_case=TRUE))
   | str_detect(column, regex("Sarafem", ignore_case=TRUE))
   | str_detect(column, regex("Sertraline", ignore_case=TRUE))
   | str_detect(column, regex("Zoloft", ignore_case=TRUE))
   | str_detect(column, regex("Fluvoxamine", ignore_case=TRUE))
   | str_detect(column, regex("Citalopram", ignore_case=TRUE))
   | str_detect(column, regex("Celexa", ignore_case=TRUE))
   | str_detect(column, regex("Escitalopram", ignore_case=TRUE))
   | str_detect(column, regex("Lexapro", ignore_case=TRUE))
   | str_detect(column, regex("Amitriptyline", ignore_case=TRUE))
   | str_detect(column, regex("Amitriptyline", ignore_case=TRUE))
   | str_detect(column, regex("Imipramine", ignore_case=TRUE))
   | str_detect(column, regex("Imipramine", ignore_case=TRUE))
   | str_detect(column, regex("Venlafaxine", ignore_case=TRUE))
   | str_detect(column, regex("Desvenlafaxine", ignore_case=TRUE))
   | str_detect(column, regex("Mirtazapine", ignore_case=TRUE))
   | str_detect(column, regex("Remeron", ignore_case=TRUE))
   | str_detect(column, regex("Trazodone", ignore_case=TRUE))
   | str_detect(column, regex("Buspirone", ignore_case=TRUE))
   | str_detect(column, regex("Buspar", ignore_case=TRUE))
   | str_detect(column, regex("Bupropion", ignore_case=TRUE))
   | str_detect(column, regex("Wellbutrin", ignore_case=TRUE))
   | str_detect(column, regex("Buproban", ignore_case=TRUE)))
}

meds_all <-
  mutate(meds_all, antidep1=as.numeric(searchfunction_antidep(med1_rxnorm_p)))%>%
  mutate(antidep2=as.numeric(searchfunction_antidep(med2_rxnorm_p))) %>%
  mutate(antidep3=as.numeric(searchfunction_antidep(med3_rxnorm_p))) %>%
  mutate(antidep4=as.numeric(searchfunction_antidep(med4_rxnorm_p))) %>%
  mutate(antidep5=as.numeric(searchfunction_antidep(med5_rxnorm_p))) %>%
  mutate(antidep6=as.numeric(searchfunction_antidep(med6_rxnorm_p))) %>%
  mutate(antidep7=as.numeric(searchfunction_antidep(med7_rxnorm_p))) %>%
  mutate(antidep8=as.numeric(searchfunction_antidep(med8_rxnorm_p))) %>%
  mutate(antidep9=as.numeric(searchfunction_antidep(med9_rxnorm_p))) %>%
  mutate(antidep10=as.numeric(searchfunction_antidep(med10_rxnorm_p))) %>%
  mutate(antidep11=as.numeric(searchfunction_antidep(med11_rxnorm_p))) %>%
  mutate(antidep12=as.numeric(searchfunction_antidep(med12_rxnorm_p))) %>%
  mutate(antidep13=as.numeric(searchfunction_antidep(med13_rxnorm_p))) %>%
  mutate(antidep14=as.numeric(searchfunction_antidep(med14_rxnorm_p))) %>%
  mutate(antidep15=as.numeric(searchfunction_antidep(med15_rxnorm_p)))

meds_all <-
  meds_all %>%
  mutate(antidepsum=rowSums(across(antidep1:antidep15), na.rm = TRUE)) %>%
  mutate(antidepYN=ifelse(antidepsum>0, 1, 0)) 

# Check on numbers
table(meds_all$antidepsum)
table(meds_all$antidepYN)

# Remove unnecessary variables (make sure don't want to search for other med classes before running this)
meds_all_tabulated <- meds_all %>% 
  select("subjectkey", "eventname", 
         stimsum, stimYN, aasum, aaYN, atomoxsum, atomoxYN, antipsychsum, antipsychYN, antidepsum, antidepYN)

##checking on numbers
medsYN_all <- meds_all_tabulated %>% 
  select("subjectkey", "eventname",
         ends_with("YN"))
medsYN_all %>% group_by(eventname) %>% 
  count(stimYN)

##make list of dataframes to merge
demographics_list <- list(demographics_all, abcd_site_all, meds_all_tabulated, family_data_baseline)

##(1) merge all demographics and admin data 
demographics_all_merged <- demographics_list %>%
  reduce(full_join, by=c("subjectkey", "eventname")) 

##subset baseline and drop NAs for all vars (except household income, educ, race)
demographics_all_baseline <- demographics_all_merged %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  drop_na(-household.income, -high.educ, -race_ethnicity, -race.6level)

demographics_all_followup <- demographics_all_merged %>% 
  filter(eventname %in% c("1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1"))

##filter out subjects in the f/u set that are not in baseline (n=28)
demographics_all_followup_filtered <- semi_join(demographics_all_followup, demographics_all_baseline,
                                                by = "subjectkey")
###all time points
demographics_complete <- bind_rows(demographics_all_baseline, demographics_all_followup_filtered)

##colSums(is.na())
 
#####################################
### cleaning behavioural data
#####################################
# colSums(is.na(df))

## (1) extract relevant timepoints and variables of interest (2) replace blanks with NA  
CBCL_data_subscales <- CBCL_data_scales %>%
  filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1")) %>%
  select("subjectkey", "eventname", "interview_age", "sex", 
        starts_with("cbcl_scr_syn_") & (ends_with("_r") | ends_with("_t"))) %>%
  mutate_all(na_if, "")

##fix variable class to be numeric
CBCL_data_subscales[,5:26] <- lapply(CBCL_data_subscales[,5:26], function(x) as.numeric(as.character(x), as.is = TRUE))

##(1) extract CBCL variables for dysregulation profile (2) calculate dysregulation profile (raw and t score) **NAs retained for now (but scale only calculated with complete) 
CBCL_data_subscales_complete <- CBCL_data_subscales %>%
  select("subjectkey", "eventname", "interview_age", "sex",  
         starts_with("cbcl_scr_syn_") & (ends_with("_r") | ends_with("_t"))) %>%
  mutate(cbcl_dysreg_raw = cbcl_scr_syn_anxdep_r + cbcl_scr_syn_attention_r + cbcl_scr_syn_aggressive_r) %>%
  mutate(cbcl_dysreg_t = cbcl_scr_syn_anxdep_t + cbcl_scr_syn_attention_t + cbcl_scr_syn_aggressive_t)
  
##extract relevant timepoints and CBCL question for CU subscale ("Doesn't seem to feel guilty after misbehaving")
CBCL_data_raw_forCU <- CBCL_data_raw %>%
  filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1")) %>%
  select("subjectkey", "eventname", "cbcl_q26_p")

##extract relevant timepoints and SDQ questions for CU subscale
SDQ_prosocial_forCU <- SDQ_prosocial %>%
  filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1")) %>%
  select("subjectkey", "eventname", starts_with("prosocial_"))

##reverse score the 3 SDQ questions
SDQ_prosocial_forCU_rev <- SDQ_prosocial_forCU %>%
  mutate(across (starts_with("prosocial_"), 
         ~ dplyr::recode(.x, '0' = 2 , '1' = 1, '2' = 0),
         ))

##(1) merge raw item scores **subjects with no total score retained for now as NA (bc scale only calculated with complete) 
CU_items_totscore <- inner_join(CBCL_data_raw_forCU, SDQ_prosocial_forCU_rev, by=c("subjectkey", "eventname")) %>%
  mutate(CUtotal = cbcl_q26_p + prosocial_q1_p + prosocial_q2_p + prosocial_q3_p)

CU_items_totscore %>% 
  select("cbcl_q26_p", "prosocial_q1_p", "prosocial_q2_p", "prosocial_q3_p") %>% 
  alpha()

#colSums(is.na(CU_items_totscore))
#describe(CU_items_totscore$CUtotal)
#table(CU_items_totscore$CUtotal)

##extract relevant timepoints and variables of interest
ksads_odd_cd <- ksads_full %>%
  filter(eventname %in% c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1")) %>%
  select("subjectkey", "eventname", 
         starts_with("ksads_15_") | starts_with("ksads_16")) %>%
  mutate_all(na_if, "")

##fix variable class to be numeric
ksads_odd_cd[,3:61] <- lapply(ksads_odd_cd[,3:61], function(x) as.numeric(as.character(x), as.is = TRUE))

#create conduct problems composite var
#(1) extract individual symptom items (2) recode 888 to 0 (3) calculate total symptoms **NAs retained for now (but scale only calculated with complete)
ksads_odd_cd_symptoms_present <- ksads_odd_cd %>%
  select("subjectkey", "eventname",
         ksads_15_95_p, ksads_15_436_p, ksads_15_435_p, ksads_15_433_p, 
         ksads_15_93_p, ksads_15_432_p, ksads_15_91_p, ksads_15_438_p, 
         ksads_15_437_p, ksads_15_434_p,ksads_16_449_p, ksads_16_463_p, 
         ksads_16_453_p, ksads_16_461_p,ksads_16_465_p, ksads_16_98_p, 
         ksads_16_104_p, ksads_16_102_p,ksads_16_457_p, ksads_16_455_p, 
         ksads_16_451_p, ksads_16_106_p, ksads_16_100_p, ksads_16_447_p, 
         ksads_16_459_p) %>%
  mutate(
    across(starts_with("ksads_"),
           ~ dplyr::recode(.x,'888' = 0 , '1' = 1, '0' = 0),
    )) %>% 
  mutate(ksads_ext_symptoms_total = rowSums(across(starts_with("ksads_"))))

ksads_odd_cd_symptoms_present %>% 
  select(ksads_15_95_p, ksads_15_436_p, ksads_15_435_p, ksads_15_433_p, 
         ksads_15_93_p, ksads_15_432_p, ksads_15_91_p, ksads_15_438_p, 
         ksads_15_437_p, ksads_15_434_p,ksads_16_449_p, ksads_16_463_p, 
         ksads_16_453_p, ksads_16_461_p,ksads_16_465_p, ksads_16_98_p, 
         ksads_16_104_p, ksads_16_102_p,ksads_16_457_p, ksads_16_455_p, 
         ksads_16_451_p, ksads_16_106_p, ksads_16_100_p, ksads_16_447_p, 
         ksads_16_459_p) %>% 
  alpha()
  


##merge 3 behavioural measures                                              
behav_list <- list(CBCL_data_subscales_complete, CU_items_totscore, ksads_odd_cd_symptoms_present)
behav_data_all_merged <- behav_list %>%
  reduce(full_join, by=c("subjectkey", "eventname"))

##subset baseline data and drop NAs
behav_data_baseline <- behav_data_all_merged %>% 
  filter(eventname == "baseline_year_1_arm_1") %>%  
  drop_na()

behav_data_followup <- behav_data_all_merged %>% 
  filter(eventname %in% c("1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1"))
##filter out subjects who are in the f/u set that are not in baseline
behav_data_followup_filtered <- semi_join(behav_data_followup, behav_data_baseline,
                                          by = "subjectkey")
##all time points
behav_data_complete <- bind_rows(behav_data_baseline, behav_data_followup_filtered)


##not part of cleaning can ignore
str(ksads_odd_cd)
colSums(is.na(behav_data_baseline))
table(ksads_odd_cd_symptoms_present$ksads_ext_symptoms_total)

#behav_data_all %>%  
 #count(eventname)

##########################################
### csv files to extract 
##########################################

#write.csv(behav_data_complete, "/path/to/csv/study3_behav_data_complete.csv", row.names = FALSE)
#write.csv(sMRI_complete_qc_baseline, "/path/to/csv/study3_sMRI_data_baseline.csv", row.names = FALSE)
#write.csv(demographics_complete, "/path/to/csv/study3_demographics_complete.csv", row.names = FALSE)
