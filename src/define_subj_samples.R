## about ----
## 
## 
## mike freund, 2020-02-19


## setup ----

set.seed(0)

library(here)
library(magrittr)
library(dplyr)
library(data.table)
library(mikeutils)

## variables

dir.ccp.hcp <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS"
dir.subsubj <- "/data/nil-external/ccp/witzermanm/AFNI_ANALYSIS_SUBSUBJECT"
dir.results <- file.path(dir.subsubj, "RESULTS_RUNWISE")

subjs <- list.dirs(dir.results, recursive = FALSE, full.names = FALSE)  ## subjs with fitted glms

subjsum <- fread(
  "C:/Users/mcf/Box/DMCC_Phase2(HCP)/Preprocessed_Data/_wrangled/summary/dmcc2_all-runs-summary.csv"
)
stroop <- fread(
  "C:/Users/mcf/Box/DMCC_Phase2(HCP)/Preprocessed_Data/_wrangled/dmcc2_behavior-and-events_stroop.csv"
)


## get bad subjects ----

pretrial.subjsum <- subjsum %>% filter(session %in% c("bas", "pro"), task == "stp", mb == "four", !grepl("(old)", subj))
subjs.bad <- pretrial.subjsum %>% 
  filter(
    per.censored > 20 | per.error > 20 | per.have.rt.stroop < 80 | missing.rows | num.gii == 0
  ) %>%
  pull(subj) %>% unique


## randomly split related individuals into two sets of unrelated subjs ----


subj.pretrial.twins <- pretrial.subjsum %>% 
  filter(!subj %in% subjs.bad) %>%
  filter(run == 1, session == "bas", !is.na(twin.pair)) %>% ## choose a run & session (doesn't matter which)
  select(subj, twin.type, twin.pair)

subjs.analysis.twins <- subj.pretrial.twins %>%
  filter(!subj %in% subjs.bad) %>%
  group_by(twin.pair) %>%
  sample_n(1) %>%
  pull(subj)

subjs.development <- setdiff(subj.pretrial.twins$subj, subjs.analysis.twins)
subjs.analysis <- pretrial.subjsum %>% filter(!subj %in% c(subjs.development, subjs.bad)) %>% pull(subj) %>% unique


## create columns in data.frames ----

pretrial <- stroop %>% filter(session %in% c("bas", "pro"), subj %in% subjsum$subj)
pretrial$subj.set <- ifelse(
  pretrial$subj %in% subjs.analysis, "analysis", 
  ifelse(
    pretrial$subj %in% subjs.development, "development",
    ifelse(
      pretrial$subj %in% subjs.bad, "bad",
      NA
    )
  )
)


pretrial.subjsum$subj.set <- ifelse(
  pretrial.subjsum$subj %in% subjs.analysis, "analysis", 
  ifelse(
    pretrial.subjsum$subj %in% subjs.development, "development",
    ifelse(
      pretrial.subjsum$subj %in% subjs.bad, "bad",
      NA
    )
  )
)


## write ----


fwrite(pretrial, here("data", "pretrial_behavior.csv"))
fwrite(pretrial.subjsum, here("data", "pretrial_subjsumm.csv"))
