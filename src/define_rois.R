## about ----
## 
## 
## mike freund, 2020-03-19


## setup ----

library(here)
library(magrittr)
library(gifti)
library(cifti)
library(abind)
library(data.table)
library(mikeutils)
library(progress)
library(profvis)

## variables

# dir.ccp.hcp <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS"
# dir.subsubj <- "/data/nil-external/ccp/witzermanm/AFNI_ANALYSIS_SUBSUBJECT"
# dir.results <- file.path(dir.subsubj, "RESULTS_RUNWISE")

pretrial <- fread(here("data", "pretrial_behavior.csv"))
pretrial.subjsum <- fread(here("data", "pretrial_subjsumm.csv"))

pretrial.subjsum[subj.set == "analysis"]$subj

pretrial.subjsum %>% View

r.vn <- readRDS(here("out", "rsa", "rmatrix_vanilla_shaefer400.rds"))  ## vanilla RSA

# str(r.vn)
# sessi <- "baseline"
# sessi.short <- "Bas"
# n.knots <- 6
# glms <- "Congruency_EVENTS_censored"
# regressors <- c("PC50InCon", "biasInCon", "PC50Con", "biasCon")
# measures <- c("corr", "eucl")
# normalizations <- c("raw", "prw")
# rsatypes <- c("vanilla", "crossval")


## ----
