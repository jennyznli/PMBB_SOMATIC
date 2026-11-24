# ========================
# CANCER IDENTIFICATION
# ========================
library(here)
source(here("R", "config.R"))
library(data.table, quietly = T)
source(here("R", "control_selection_f4.R"))
source(here("R", "sample_selection_f4.R"))

pmbb4 <- here("Documentation", "PMBB", "4.0")

# PMBB
cov <- fread(here(pmbb4, "PMBB-Release-2026-4.0_phenotype_covariates.txt"), header = TRUE)
person <- fread(here(pmbb4, "PMBB-Release-2026-4.0_phenotype_person.txt"), header = TRUE)
exome_ids <- read.table(here(pmbb4, "PMBB-Release-2026-4.0_genetic_exome.sample_list.txt"), header = FALSE)$V1

obs <- fread(here(pmbb4, "PMBB-Release-2026-4.0_phenotype_observation.txt"), header = TRUE)
cond <- fread(here(pmbb4, "PMBB-Release-2026-4.0_phenotype_condition_occurrence.txt"), header = TRUE)
substance <- fread(here(pmbb4, "PMBB-Release-2026-4.0_phenotype_substance_use.txt"), header = TRUE)

pmcr_all <- fread(here(pmbb4, "PMBB-Release-2026-4.0_phenotype_cancer_pmcr.txt"), header = TRUE)
brca_all <- fread(here(pmbb4, "PMBB-Release-2026-4.0_phenotype_cancer_brca.txt"), header = TRUE)

# progeny 
progeny_breast <- read_excel(here("Projects", "02-simplexo", "ss", "br_pts_for_exwas_10022025.xlsx")) |> filter(DNA == "D")
progeny_tumor <- read_excel(here("Projects", "03-ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Tumor") 
progeny <- read_excel(here("Projects", "03-ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Data_from_master_table") |> filter(DNA == "D")
progeny_breast$SampNum <- as.numeric(progeny_breast$SampNum)
progeny_tumor$SampNum <- as.numeric(progeny_tumor$SampNum)
progeny$SampNum <- as.numeric(progeny$SampNum)
progeny_tumor$globalid <- as.numeric(progeny_tumor$globalid)
progeny$globalid <- as.numeric(progeny$globalid)

up <- read.csv(here("Projects", "02-simplexo", "data", "simplexo_up_map.csv"))

# id lists
aml <- read_excel(here("Projects", "03-ch", "ss", "20251110_aml_ngs.xlsx"))
cov |> filter(person_id %in% aml$PMBB_ID) |> select(batch) |> unique()
# pilot study has all in freeze 2, so can't really tell...
# PMBB3278234229783 and PMBB9930975965330 are controls
# PMBB9930975965330 - 47.9 

# ========================
# PMBB
# ========================
icd_hem <- read_excel(here("Documentation", "ICD", "ICD.xlsx"), sheet = "Hematologic_Long")
hem <- paste0("^", icd_hem$code)

malig_insitu <- c(
  # malignant neoplasms except non-melanoma skin cancer
  "^C(?!44)",           # ICD10
  "^Z85(?!\\.828)",     # ICD10 personal history
  "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])",  # ICD9 
  "^V10(?!\\.83)",       # ICD9 personal history
  
  # in situ neoplasms
  "^D0[0-9]",     # ICD10
  "^Z86.00",      # ICD10 personal history
  "^(?:23[0-4])"   # ICD-9
)

cancer_codes <- c(malig_insitu, hem)
length(cancer_codes)
# 70 

dim(cov)
# 70925
cancers <- select_samples_f4(
  sample_name = "hem_solid",
  icd_codes = cancer_codes,
  min_instances = 1,
  exclude = FALSE,
  pmbb_dir = here("Documentation", "PMBB", "4.0"),
  data_dir = here("Projects", "00-pmbb", "data"),
  log_dir = here("Projects", "00-pmbb", "log")
)

icd_ids <- sort(unique(cancers$filtered_patients$person_id))
length(icd_ids)
# 37245
# 70925 - 37245 = 33680 controls 

controls <- select_controls_f4(
  control_name = "hem_solid",
  exclude_codes = cancer_codes,
  age_filter = NULL,
  pmbb_dir = here("Documentation", "PMBB", "4.0"),
  data_dir = here("Projects", "00-pmbb", "data"),
  log_dir = here("Projects", "00-pmbb", "log")
)

control_icd_ids <- sort(unique(controls$final_controls$person_id))
length(control_icd_ids)
# 33680

# 70925 - 37245 = 33680 matches up 
dim(obs_only)
# 7839 
# 
# infer age of cancer diagnosis with first date of diagnosis and birth date
icd_age <- cancers$all_patients |> 
  select("person_id", "first_date") |> 
  left_join(person[, c("birth_datetime", "person_id")], by = "person_id") |>
  mutate(icd_age = as.numeric(difftime(first_date, birth_datetime, units = "days")) / 365.25)
colnames(icd_age) <- c("PMBB_ID", "CaDxDate_ICD", "Birth_Date", "CaDxAge_ICD")

write.csv(cancers$filtered_patients, here("Projects", "00-pmbb", "data", "pmbb_icd_hem_solid_cancer_df.csv"), row.names = FALSE)
write.csv(icd_ids, here("Projects", "00-pmbb", "data", "pmbb_icd_hem_solid_cancer_ids.csv"))

# ========================
# PMCR
# ========================
colnames(pmcr_all)
# FirstDxDate

pmcr <- merge_duplicates(pmcr_all, "person_id")
pmcr$CaDxDate_PMCR <- as.Date(sapply(pmcr$FirstDxDate, get_earliest_date, USE.NAMES = FALSE), origin = "1970-01-01")

pmcr_ids <- unique(pmcr$person_id)
length(pmcr_ids)
# 20192

pmcr_age <- pmcr |>
  select(PMBB_ID = person_id, CaDxDate_PMCR) |>
  left_join(person |> select(PMBB_ID = person_id, Birth_Date = birth_datetime), by = "PMBB_ID") |>
  mutate(CaDxAge_PMCR = as.numeric(difftime(CaDxDate_PMCR, as.Date(Birth_Date), units = "days")) / 365.25)
colnames(pmcr_age) <- c("PMBB_ID", "CaDxDate_PMCR", "Birth_Date", "CaDxAge_PMCR")

# ========================
# BRCA
# ========================
colnames(brca_all)
# DateDiagnosis

brca <- merge_duplicates(brca_all, "person_id")
brca$CaDxDate_BRCA <- as.Date(sapply(brca$DateDiagnosis, get_earliest_date, USE.NAMES = FALSE), origin = "1970-01-01")

brca_ids <- unique(brca$person_id)
length(brca_ids)
# 3124

brca_age <- brca |>
  select(PMBB_ID = person_id, CaDxDate_BRCA) |>
  left_join(person |> select(PMBB_ID = person_id, Birth_Date = birth_datetime), by = "PMBB_ID") |>
  mutate(CaDxAge_BRCA = as.numeric(difftime(CaDxDate_BRCA, as.Date(Birth_Date), units = "days")) / 365.25)
colnames(brca_age) <- c("PMBB_ID", "CaDxDate_BRCA", "Birth_Date", "CaDxAge_BRCA")

# ========================
# PROGENY - BRCA1/2
# ========================
### TUMOR ###
progeny_tumor_merged <- merge_duplicates(progeny_tumor, "globalid")

progeny_tumor_merged$CaDxAge_Progeny <- sapply(progeny_tumor_merged$CaDxAge, get_earliest_age)

progeny_tumor_merged$CaDxDate_Progeny <- as.Date(sapply(progeny_tumor_merged$CaDxDate, get_earliest_date, USE.NAMES = FALSE), origin = "1970-01-01")

### MASTER ###
# filter down to PMBB 
progeny_pmbb <- progeny |>
  left_join(up, by = "SampNum") |>
  filter(!is.na(PMBB_ID))
length(unique(progeny_pmbb$PMBB_ID))
# 1278

# merge tumor and master
progeny_pmbb_tumor <- inner_join(progeny_pmbb, progeny_tumor_merged, by = "globalid") |>
  select(PMBB_ID, CaDxDate, CaDxAge, CaDxDate_Progeny, CaDxAge_Progeny)

# ========================
# PROGENY - BREAST (SIMPLEXO)
# ========================
dim(progeny_breast)
# 5046 
progeny_breast_merged <- merge_duplicates(progeny_breast, "SampNum")
dim(progeny_breast_merged)
# 4273

progeny_breast_merged$CaDxAge_Progeny <- sapply(progeny_breast_merged$CaDxAge, get_earliest_age)

progeny_breast_merged$CaDxDate_Progeny <- as.Date(sapply(progeny_breast_merged$CaDxDate, get_earliest_date, USE.NAMES = FALSE), origin = "1970-01-01")

progeny_breast_pmbb <- progeny_breast_merged |>
  left_join(up, by = "SampNum") |>
  filter(!is.na(PMBB_ID)) |> 
  select(PMBB_ID, CaDxDate, CaDxAge, CaDxDate_Progeny, CaDxAge_Progeny)

# ========================
# MERGE PROGENY SHEETS 
# ========================
progeny_pmbb <- rbind(progeny_pmbb_tumor, progeny_breast_pmbb) |>
  merge_duplicates("PMBB_ID") 

progeny_pmbb$CaDxAge_Progeny  <- sapply(progeny_pmbb$CaDxAge_Progeny, get_earliest_age, USE.NAMES = FALSE)
progeny_pmbb$CaDxDate_Progeny <- as.Date(sapply(progeny_pmbb$CaDxDate_Progeny, get_earliest_date, USE.NAMES = FALSE),
                                         origin = "1970-01-01")

# they differ bc exwas only covers breast cancers...

progeny_ids <- progeny_pmbb$PMBB_ID
length(progeny_ids)
# 1739 

# ========================
# WENTING PULL
# ========================
# PMBB
f4 <- fread(here("Documentation", "CREP", "20250928_full_key.csv"), header = TRUE)
# 70925    18

crep <- f4 |> filter(CREP == TRUE) |> select(PMBB_ID, CREP_ID) |> arrange(CREP_ID)
dim(crep)
# 2835 
write.csv(crep, here("Documentation", "CREP", "20261123_f4_crep_sampnum_map.csv"), row.names = FALSE)

# ========================
# FINAL CANCER + CONTROL IDS
# ========================
final_cancer_ids <- sort(unique(c(icd_ids, pmcr_ids, brca_ids, progeny_ids)))
length(final_cancer_ids)
# 38127

control_ids <- setdiff(cov$person_id, final_cancer_ids)
length(control_ids)
# 32800

write.csv(final_cancer_ids, here("Projects", "00-pmbb", "data", "pmbb_all_source_cancer_ids.csv"), row.names = FALSE)
write.csv(control_ids, here("Projects", "00-pmbb", "data", "pmbb_all_source_control_ids.csv"), row.names = FALSE)

# ========================
# COMPILE  SOURCES
# ========================
# create merged dataframe
cancer_df <- data.frame(PMBB_ID = final_cancer_ids)

# diagnosis age based on priority: PMCR > BRCA > Progeny > ICD 
cancer_df$CaDxAge <- coalesce(
  pmcr_age$CaDxAge_PMCR[match(cancer_df$PMBB_ID, pmcr_age$PMBB_ID)],
  brca_age$CaDxAge_BRCA[match(cancer_df$PMBB_ID, brca_age$PMBB_ID)],
  progeny_pmbb$CaDxAge_Progeny[match(cancer_df$PMBB_ID, progeny_pmbb$PMBB_ID)],
  icd_age$CaDxAge_ICD[match(cancer_df$PMBB_ID, icd_age$PMBB_ID)]
)

# ========================
# ADD SAMPLE AGES
# ========================
cancer_df$sample_age <- cov$sample_age[match(cancer_df$PMBB_ID, cov$person_id)]

cancer_df <- cancer_df |> filter(!is.na(sample_age), !is.na(CaDxAge))
dim(cancer_df)
# 38116

write.csv(cancer_df, here("Projects", "00-pmbb", "data", "pmbb_all_sources_cancer_age.csv"))

# ========================
# SMOKING HISTORY
# ========================
### PMBB ###
table(substance$TOBACCO_USER)
# Never Not Asked   Passive      Quit       Yes 
# 6   1664807     43354      8886   1122481    234893

# turn "Not Asked" and "" into NA 
substance <- substance |> mutate(TOBACCO_USER = ifelse(TOBACCO_USER %in% c("Not Asked", ''), NA, TOBACCO_USER))

table(substance$TOBACCO_USER, useNA = "ifany")
# Never Passive    Quit     Yes    <NA> 
#   1664807    8886 1122481  234893   43360 

# collapse down values in this column
merge_duplicates2 <- function(df, group_col) {
  df |>
    group_by(across(all_of(group_col))) |>
    summarise(
      across(-any_of(group_col), ~ {
        non_na_vals <- unique(na.omit(as.character(.x)))
        if(length(non_na_vals) == 0) NA_character_ else paste(non_na_vals, collapse = ";")
      }),
      .groups = "drop"
    )
}
substance_merge <- merge_duplicates2(substance, "person_id")
table(substance_merge$TOBACCO_USER, useNA = "ifany")

# non smokers will only have Never and no other values
pmbb_non_smoker <- substance_merge |> filter(TOBACCO_USER == "Never")
pmbb_non_smoker_ids <- sort(unique(pmbb_non_smoker$person_id))
length(pmbb_non_smoker_ids)
# 25677

# smokers will have anything other than Never and can't be empty
pmbb_smoker <- substance_merge |> filter(!is.na(TOBACCO_USER) & TOBACCO_USER != "Never")
pmbb_smoker_ids <- sort(unique(pmbb_smoker$person_id))
length(pmbb_smoker_ids)
# 31058

### PROGENY ###
colnames(progeny)
dim(progeny)
# 4336

progeny_smk_clean <- progeny |>
  inner_join(up, by = "SampNum")
dim(progeny_smk_clean)
# 1278

# merge to individual level
progeny_smk <- progeny_smk_clean |> merge_duplicates("SampNum")
table(progeny_smk$SmokingEver, useNA = "ifany")
# Never      No Unknown     Yes    <NA>
#     775      62      22     367      52

# get rid of unknowns or NA
progeny_smk_pmbb <- progeny_smk |>
  filter(SmokingEver != "Unknown", !is.na(SmokingEver))

# non smokers have Never, No
progeny_non_smoker <- progeny_smk_pmbb |> filter(SmokingEver %in% c("Never", "No"))
progeny_non_smoker_ids <- sort(unique(progeny_non_smoker$person_id))
length(progeny_non_smoker_ids)
# 837

# smokers have yes
progeny_smoker <- progeny_smk_pmbb |> filter(SmokingEver == "Yes")
progeny_smoker_ids <- sort(unique(progeny_smoker$person_id))
length(progeny_smoker_ids)
# 367

### MERGED SMOKERS ###
cov_smoke <- cov |> mutate(
  PMBB_SmokeHistory = case_when(
    person_id %in% pmbb_smoker_ids ~ 1,
    person_id %in% pmbb_non_smoker_ids ~ 0,
    TRUE ~ NA_real_
  ),
  Progeny_SmokeHistory = case_when(
    person_id %in% progeny_smoker_ids ~ 1,
    person_id %in% progeny_non_smoker_ids ~ 0,
    TRUE ~ NA_real_
  )
)
table(cov_smoke$PMBB_SmokeHistory, useNA = "ifany")
# 0     1  <NA> 
# 36663 31058  3204 

table(cov_smoke$Progeny_SmokeHistory, useNA = "ifany")
# 0     1  <NA>
# 837   367 69722

# investigate concordance
table(PMBB = cov_smoke$PMBB_SmokeHistory,
      Progeny = cov_smoke$Progeny_SmokeHistory,
      useNA = "ifany")
# PMBB       0     1  <NA>
#   0      506    27 36130
# 1       15   210 30833
# <NA>   315   130  2759

### MERGE SOURCES ###
# if any person has smoking history -> smoker
# if they are non smokers in both -> non smoker
all_smoker_ids <- sort(union(pmbb_smoker_ids, progeny_smoker_ids))
all_non_smoker_ids <- setdiff(union(pmbb_non_smoker_ids, progeny_non_smoker_ids), all_smoker_ids)

length(all_smoker_ids)
# 31215
length(all_non_smoker_ids)
# 36952

# ========================
# PREPARE MATCHING DATASET
# ========================
cov1 <- cov %>%
  mutate(
    Smoke_History = case_when(
      person_id %in% all_smoker_ids ~ 1,
      person_id %in% all_non_smoker_ids ~ 0,
      TRUE ~ NA_real_
    ),
    Cancer = ifelse(person_id %in% final_cancer_ids, 1, 0),
  ) %>% left_join(
    cancer_df[1:2], by = c("person_id" = "PMBB_ID")
  )

table(cov1$Smoke_History, useNA = "ifany")
# 0     1  <NA> 
# 36951 31215  2759 

# ========================
# CREATE CANCER STRATA
# ========================
cov2 <- cov1 %>%
  mutate(
    Strata = case_when(
      # Strata 1: No cancer
      Cancer == 0 ~ 1,
      # Strata 2: Cancer developed >1 year AFTER sample collection
      Cancer == 1 & !is.na(CaDxAge) & CaDxAge > (sample_age + 1) ~ 2,
      # Strata 3: Cancer diagnosed BEFORE/AT sample collection
      Cancer == 1 & !is.na(CaDxAge) & CaDxAge <= sample_age ~ 3,
      # Strata 4: Cancer within 1 year of sample OR missing dx age
      Cancer == 1 ~ 4,
      # NA - missing DX age
      TRUE ~ NA_real_
    ),
    Strata_Label = case_when(
      Strata == 1 ~ "No cancer",
      Strata == 2 ~ "Post-sample cancer (>1yr)",
      Strata == 3 ~ "Pre-sample cancer",
      Strata == 4 ~ "Uncertain",
      TRUE ~ NA_character_
    )
  )
# ========================
# CHECK STRATA DISTRIBUTION
# ========================
table(cov2$Strata_Label, useNA = "ifany")
# No cancer Post-sample cancer (>1yr)         Pre-sample cancer
# 27564                      3399                     23596
# Uncertain
# 2611

table(cov2$Class, useNA = "ifany")
# AFR     AMR     EAS     EUR     SAS UNKNOWN
# 12113     761     907   42133     743     513

table(cov2$Smoke_History, useNA = "ifany")
# 0     1  <NA>
# 28657 25798  2715

dim(cov2)
cov3 <- cov2 %>% filter(!is.na(Smoke_History), !is.na(sample_age))
dim(cov3)
# 68157

# ========================
# UNAFFECTED ANALYSIS
# ========================
unaff <- cov3 %>% filter(Strata %in% c(1, 2))

dim(unaff)
# 38595    21

table(unaff$batch)
# 1     2
# 23384  5438

table(unaff$sequenced_gender)
# Female   Male
# 14684  14138

table(unaff$Smoke_History)
# 0     1
# 16000 12822

range(unaff$sample_age)
# 18 - 100.5

# ========================
# IDENTIFY 2 PANEL OF NORMALS
# ========================
table(cov$batch)
# 1     2     3 
# 43834 13246 13845  

# 1 2 3
# 1 = F2, 2 = F3, 3 = F4

### FREEZE 2 ###
# get the youngest 50 people in freeze 2 "1"
youngest_f2 <- cov3 |> filter(Strata == 1, Smoke_History == 0, batch == 1) |>
  arrange(sample_age) |>
  slice_head(n = 50)
youngest_f2_ids <- sort(youngest_f2$person_id)
table(youngest_f2$crep_highrisk_flag)
# join flags and see if any ar ein high risk ccohrots? 
# none are in crep obviously

colnames(youngest_f2)
# check ancestry distribution
prop.table(table(youngest_f2$Class))
# AFR     AMR     EAS     EUR     SAS UNKNOWN
# 0.38    0.02    0.02    0.52    0.02    0.04

prop.table(table((unaff |> filter(batch == 1))$Class))
# AFR        AMR        EAS        EUR        SAS    UNKNOWN
# 0.27373418 0.01492473 0.01821758 0.66412932 0.01761888 0.01137530

# aml_controls <- c("PMBB3278234229783", "PMBB9930975965330")
# f2_pon <- sort(c(youngest_f2_ids[1:78], aml_controls))

range(youngest_f2$sample_age)
# 18.0 18.8

write.csv(youngest_f2_ids, here("Projects", "00-pmbb", "data", "pon_f2_pon_ids.csv"), row.names = FALSE)

### FREEZE 3 / 4 - same platform ###
# exclude CREP just in case... 
youngest_f3 <- cov3 |> filter(Strata == 1, Smoke_History == 0, batch %in% c(2, 3), !crep_highrisk_flag) |>
  arrange(sample_age) |>
  slice_head(n = 50)
youngest_f3_ids <- sort(youngest_f3$person_id)

# check ancestry distribution
prop.table(table(youngest_f3$Class))
# AFR  AMR  EUR  SAS  UNK 
# 0.30 0.06 0.58 0.04 0.02 

prop.table(table((unaff |> filter(batch == 2))$Class))
# AFR        AMR        EAS        EUR        SAS    UNKNOWN
# 0.06914307 0.01930857 0.02427363 0.85582935 0.01838911 0.01305627
# hold on ... 

range(youngest_f3$sample_age)
# 18.1 19.7
write.csv(youngest_f3_ids, here("Projects", "00-pmbb", "data", "pon_f3_pon_ids.csv"), row.names = FALSE)

### save end df
pon_ids <- c(youngest_f2_ids, youngest_f3_ids)

pon_df <- data.frame(
  PMBB_ID = pon_ids,
  batch   = ifelse(pon_ids %in% youngest_f2_ids, "2", "3"),
  stringsAsFactors = TRUE
)
write.csv(pon_df, here("Projects", "00-pmbb", "data", "pon_f23_pon_df.csv"), row.names = FALSE)






