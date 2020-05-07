## about ----
## 
## creates model similarity matrices (for RSA) and writes them to .csv files.
## 
## mike freund, 2020-05-01
##
## 

## setup ----

library(here)
library(mikeutils)
library(magrittr)
library(dplyr)

## varibles

source(here("src", "setup.R"))


## between run models ----

rsmods.vanil <- setNames(vector("list", length(rsmods.vanil.names)), rsmods.vanil.names)

for (ii in seq_along(rsmods.vanil)) {
  
  rsmods.vanil[[ii]] <- matrix(
    0, 
    nrow = length(conds.run), ncol = length(conds.run), 
    dimnames = list(.row = conds.run, .col = conds.run)
    )
  
}

is.incon <- grepl("InCon", conds.run)
is.congr <- grepl("PC50Con|biasCon", conds.run)
is.pc50  <- grepl("PC50", conds.run)
is.bias  <- grepl("bias", conds.run)
is.run1  <- grepl("run1", conds.run)
is.run2  <- grepl("run2", conds.run)

diag(rsmods.vanil$major) <- 1  ## major

# diag(rsmods.vanil$minor[is.run1, is.run2]) <- 1  ## minor
# diag(rsmods.vanil$minor[is.run2, is.run1]) <- 1

rsmods.vanil$InCon[is.incon, is.incon] <- 1  ## incon

rsmods.vanil$Con[is.congr, is.congr] <- 1  ## congr

rsmods.vanil$PC50[is.pc50, is.pc50] <- 1  ## pc50

rsmods.vanil$bias[is.bias, is.bias] <- 1  ## bias

rsmods.vanil$runw[is.run1, is.run1] <- 1  ## run, within
rsmods.vanil$runw[is.run2, is.run2] <- 1

rsmods.vanil$runb[is.run2, is.run1] <- 1  ## run, between
rsmods.vanil$runb[is.run1, is.run2] <- 1


## within run models ----

rsmods.cross <- rsmods.vanil

for (ii in seq_along(rsmods.cross)) {
  
  rsmods.cross[[ii]] <- rsmods.cross[[ii]][1:4, 1:4]
  colnames(rsmods.cross[[ii]]) <- gsub("_run1", "", colnames(rsmods.cross[[ii]]))
  rownames(rsmods.cross[[ii]]) <- gsub("_run1", "", rownames(rsmods.cross[[ii]]))
  
}

rsmods.cross[c("runw", "runb")] <- NULL
# rsmods.cross[c("minor", "run.betw", "run.with")] <- NULL


## check ----

lapply(rsmods.vanil, qcor)
lapply(rsmods.cross, qcor)


## write ----

lapply(
  names(rsmods.vanil), function(name)
    write.csv(rsmods.vanil[[name]], here("out", "rsa", "models", paste0("rsmodel_vanil_", name, ".csv")))
)

lapply(
  names(rsmods.cross), function(name)
    write.csv(rsmods.cross[[name]], here("out", "rsa", "models", paste0("rsmodel_cross_", name, ".csv")))
)
