# ========================
# PANEL OF NORMALS: COHORT + LAB SUMMARY
# ========================
library(here)
source(here("R", "config.R"))
library(data.table, quietly = T)

pmbb4 <- here("Documentation", "PMBB", "4.0")
out_dir <- here("Projects", "00-pmbb", "data")

pon_ids <- read.csv(here(out_dir, "pon_f23_pon_df.csv"))$PMBB_ID
length(pon_ids)
# 100

cov <- fread(here(pmbb4, "PMBB-Release-2026-4.0_phenotype_covariates.txt"), header = TRUE)

aml <- read_excel(here("Projects", "03-ch", "ss", "20251110_aml_ngs.xlsx"))
cov |> filter(person_id %in% aml$PMBB_ID) |> select(batch) |> unique()
# freeze 2 only 

# ========================
# COMPARE COVS WITH PMBB
# ========================
pon <- cov |>
  filter(person_id %in% pon_ids) |>
  select(person_id, batch, sequenced_gender, Class, sample_date, sample_age)
dim(pon)
# 100

# sanity: every PON id found in cov
setdiff(pon_ids, pon$person_id)

prop.table(table(pon$sequenced_gender))
# Female   Male
# .69     .31
prop.table(table(cov$sequenced_gender))
# Female      Male
# 0.5206063 0.4793937

prop.table(table(pon$Class))
# AFR  AMR  EAS  EUR  SAS  UNK
# 0.34 0.04 0.01 0.55 0.03 0.03
prop.table(table(cov$Class))
# AFR         AMR         EAS         EUR         SAS         UNK
# 0.210461755 0.014649277 0.018794501 0.731293620 0.015227353 0.009573493

# ========================
# BLOOD COUNTS / LABS
# ========================
labs <- c(
  wbc       = "wbc",
  wbc_urine = "wbc_urine",
  rbc       = "rbc",
  rbc_urine = "rbc_urine",
  platelets = "platelets",
  hemoglobin = "hemoglobin",
  hba1c     = "hba1c"
)

read_lab <- function(lab) {
  # fix column types: fread guesses logical for all-empty columns, which breaks bind_rows across labs
  fread(
    here(pmbb4, paste0("PMBB-Release-2026-4.0_phenotype_labs_", lab, ".txt")),
    header = TRUE,
    colClasses = list(
      character = c("person_id", "measurement_date", "measurement_datetime",
                    "measurement_source_value", "unit_source_value",
                    "value_source_value", "unit_converted"),
      numeric   = c("value_as_number", "value_converted", "range_low", "range_high")
    )
  ) |>
    mutate(across(any_of("outlier"), as.character))
}

pon_lab_rows <- list()     # every measurement for PON people
pon_lab_person <- list()   # one row per PON person per lab
lab_compare <- list()      # PON vs whole PMBB, per lab

for (lab in names(labs)) {
  lab_all <- read_lab(labs[[lab]]) |> filter(!is.na(value_as_number))
  lab_all$unit_converted <- lab_all$unit_converted
  lab_all$outlier <- lab_all$outlier
  
  # urine wbc/rbc files have no outlier column
  if (!"outlier" %in% names(lab_all)) lab_all$outlier <- NA_character_

  # one value per person in the whole cohort (median of their measurements) for comparison
  cohort_person <- lab_all |>
    group_by(person_id) |>
    summarise(median_value = median(value_as_number), .groups = "drop") |>
    mutate(group = ifelse(person_id %in% pon_ids, "PON", "PMBB"))

  # all PON rows, with days from blood draw (sample_date)
  rows <- lab_all |>
    filter(person_id %in% pon_ids) |>
    left_join(pon |> select(person_id, sample_date), by = "person_id") |>
    mutate(
      lab = lab,
      days_from_sample = as.numeric(as.Date(measurement_date) - as.Date(sample_date))
      # out_of_range = case_when(
      #   is.na(value_as_number) ~ NA,
      #   !is.na(range_low)  & value_as_number < range_low  ~ TRUE,
      #   !is.na(range_high) & value_as_number > range_high ~ TRUE,
      #   TRUE ~ FALSE
      # )
    )
  pon_lab_rows[[lab]] <- rows

  # per PON person: count, median, and the measurement closest to the blood draw
  pon_lab_person[[lab]] <- rows |>
    # same test charted twice (same person, same time, same value) -> count once
    distinct(person_id, measurement_datetime, value_as_number, .keep_all = TRUE) |>
    group_by(person_id) |>
    summarise(
      lab = .env$lab,
      n_measurements = n(),
      n_dates = n_distinct(measurement_date),
      n_visits = n_distinct(visit_occurrence_id, na.rm = TRUE),
      median_value = median(value_as_number),
      min_value = min(value_as_number),
      max_value = max(value_as_number),
      # n_out_of_range = sum(out_of_range, na.rm = TRUE),
      n_outlier_flag = sum(outlier == "1", na.rm = TRUE),
      # first() returns NA instead of erroring when days_from_sample is all NA
      closest_value = first(value_as_number[which.min(abs(days_from_sample))]),
      closest_days_from_sample = first(days_from_sample[which.min(abs(days_from_sample))]),
      unit = paste(sort(unique(unit_converted)), collapse = ";"),
      .groups = "drop"
    )

  lab_compare[[lab]] <- cohort_person |>
    group_by(group) |>
    summarise(
      lab = .env$lab,
      n_people = n(),
      median = median(median_value),
      q25 = quantile(median_value, 0.25),
      q75 = quantile(median_value, 0.75),
      .groups = "drop"
    )

  rm(lab_all, cohort_person); gc()
}

pon_lab_rows <- bind_rows(pon_lab_rows)
pon_lab_person <- bind_rows(pon_lab_person)
lab_compare <- bind_rows(lab_compare)

# how many PON people have each lab at all
pon_lab_person |> count(lab)
# lab            n
# <chr>      <int>
# 1 hba1c         40
# 2 hemoglobin    92
# 3 platelets     93
# 4 rbc           91
# 5 wbc           91
# 6 wbc_urine     13

# PON vs PMBB per lab
lab_compare

# ========================
# PON-LEVEL WIDE DF (1 row per PON person)
# ========================
# cols named <lab>_<stat>, e.g. wbc_median, wbc_closest_days, wbc_unit
pon_lab_wide <- pon_lab_person |>
  select(person_id, lab, n_visits, n_dates, median = median_value,
         closest = closest_value, closest_days = closest_days_from_sample, unit) |>
  tidyr::pivot_wider(
    names_from = lab,
    values_from = c(median, closest, closest_days, n_visits, n_dates, unit),
    names_glue = "{lab}_{.value}"
  )

# start from all 100 PON so people with no labs still get a row (NAs)
pon_summary <- pon |>
  left_join(pon_lab_wide, by = "person_id") |>
  select(person_id, batch, sequenced_gender, Class, sample_date, sample_age,
         sort(setdiff(names(pon_lab_wide), "person_id")))

dim(pon_summary)
# want 100 rows
stopifnot(nrow(pon_summary) == nrow(pon), !anyDuplicated(pon_summary$person_id))

write.csv(pon_summary, here(out_dir, "20261124_pon_summary.csv"), row.names = FALSE)

write.csv(pon_lab_rows,   here(out_dir, "20261124_pon_labs_all_rows.csv"),   row.names = FALSE)
write.csv(pon_lab_person, here(out_dir, "20261124_pon_labs_per_person.csv"), row.names = FALSE)
write.csv(lab_compare,    here(out_dir, "20261124_pon_labs_vs_pmbb.csv"),    row.names = FALSE)

# ========================
# CONDITION OCCURRENCES 
# ========================
# autoimmune disease, HIV, cardiovascular disease

library(duckdb)
con <- dbConnect(duckdb())
dbGetQuery(con, sprintf("
  SELECT DISTINCT co.person_id
  FROM read_csv_auto('%s', delim='\t') co
  JOIN read_csv_auto('%s', delim='\t') ca
    ON co.condition_concept_id = ca.descendant_concept_id
  WHERE ca.ancestor_concept_id IN (<parent ids>)",
                        here(pmbb4, "PMBB-Release-2026-4.0_phenotype_condition_occurrence.txt"),
                        here(pmbb4, "PMBB-Release-2026-4.0_phenotype_concept_ancestor.txt")))

# one row per visit (lab rows with a visit id only)
visit_span <- pon_lab_rows |>
  filter(!is.na(visit_occurrence_id)) |>
  group_by(lab, person_id, visit_occurrence_id) |>
  summarise(
    n_dates = n_distinct(measurement_date),
    span_days = as.numeric(max(as.Date(measurement_date)) - min(as.Date(measurement_date))),
    .groups = "drop"
  )

# per lab: how many visits cover 1 day vs several
visit_span |>
  mutate(multi_day = n_dates > 1) |>
  count(lab, multi_day)

# distribution of how many days a multi-day visit covers
visit_span |>
  filter(n_dates > 1) |>
  group_by(lab) |>
  summarise(n_visits = n(), median_span = median(span_days), max_span = max(span_days), .groups = "drop")



