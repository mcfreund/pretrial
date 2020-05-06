## about ----
## contains strings that are commonly used by scripts in this project.
## 
## mike freund, 2020-05-1


## data ----

pretrial <- data.table::fread(here("data", "pretrial_behavior.csv"))
pretrial.subjsum <- data.table::fread(here("data", "pretrial_subjsumm.csv"))

## paths: pointers for atlases and for BOLD timeseries ----

nodename <- Sys.info()["nodename"]

if (nodename == "ccplinux1") {

  dir.atlas <- "/data/nil-external/ccp/freund/atlases"
  dir.schaefer <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"

} else if (nodename == "CCP-FREUND") {
  ## mike freund's (i.e., ccp's) thinkpad
  ## reliant on box drive
  ## assumes box drive location at ./Users/mcf/Box

  dir.atlas <- "C:/local/atlases"
  dir.schaefer <- dir.atlas

}


## factor levels ----

rsmods.vanil.names <- c("PC50", "bias", "InCon", "Con", "runw", "runb", "major")
rsmods.cross.names <- c("PC50", "bias", "InCon", "Con", "major")
modname <- c(paste0("vanil_", rsmods.vanil.names), paste0("cross_", rsmods.cross.names))

conds <- c("PC50InCon", "biasInCon", "PC50Con", "biasCon")
conds.run <- c(
  "PC50InCon_run1", "biasInCon_run1", "PC50Con_run1", "biasCon_run1", 
  "PC50InCon_run2", "biasInCon_run2", "PC50Con_run2", "biasCon_run2"
  )

normal <- c("raw", "prw")
session <- c("baseline", "proactive")
measure <- c("corr", "eucl", "neuc")
method <- c("vanil", "cross")
method.long <- c("vanilla", "crossva")
glmname <- c("Congruency_EVENTS_censored", "Congruency_EVENTS_tentzero01210_censored")
# modname <- c(
#   "cross_bias",
#   "cross_Con",
#   "cross_InCon",
#   "cross_major",
#   "cross_PC50",
#   "vanil_bias",
#   "vanil_Con",
#   "vanil_InCon",
#   "vanil_major",
#   "vanil_minor",
#   "vanil_PC50",
#   "vanil_runw",
#   "vanil_runb"
# )


subjs.analysis <- unique(pretrial.subjsum[subj.set == "analysis"]$subj)
subjs.development <- unique(pretrial.subjsum[subj.set == "development"]$subj)
subjs.bad <- unique(pretrial.subjsum[subj.set == "bad"]$subj)

networks <- c("Vis", "DorsAttn", "SalVentAttn", "Cont", "Default", "SomMot", "Limbic")

## models ----

## read models
# 
# models <- lapply(
#   modname,
#   function(.) data.table::fread(here("out", "rsa", "models", paste0("rsmodel_", ., ".csv")), data.table = FALSE)
# )
# names(models) <- modname


## atlas ----

library(cifti)
parcellation <- read_atlas("schaefer400")


## palettes ----

# colors.profs <- c(
#   incon = "#d95f02ff", targt = "#1b9e77ff", distr = "#7570b3ff", targt.incon = "#26190dff", targt.distr = "#4682b4ff"
# )
