
library(here)
library(magrittr)
library(dplyr)
library(data.table)
library(mikeutils)
library(cifti)
library(gifti)
library(abind)
library(progress)


## data ----

pretrial <- fread(here("data", "pretrial_behavior.csv"))  ## behavior and events
pretrial.subjsum <- fread(here("data", "pretrial_subjsumm.csv"))  ## summary of behavior and events / subj

## paths ----

dir.ccp.hcp <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS"
dir.subsubj <- file.path("/data/nil-external/ccp/freund/sub-subj-glms/runwise_old")

subjs <- intersect(
  list.dirs(dir.subsubj, recursive = FALSE, full.names = FALSE),
  list.dirs(dir.ccp.hcp, recursive = FALSE, full.names = FALSE)
)

sessi <- c("baseline", "proactive")
sessi.short <- c("Bas", "Pro")
glms <- "Congruency_EVENTS_censored"
n.knots <- 6
regressors <- c("PC50InCon", "biasInCon", "PC50Con", "biasCon")
networks <- c("Vis", "DorsAttn", "SalVentAttn", "Cont", "Default", "SomMot", "Limbic")

## refine list of subjects

## read this guy in just to get the subjs for which we presently have stats:
r.2tr1knot <- abind(
  readRDS(here("out", "rsa", "rmatrix_vanilla_corr_shaefer400_baseline_Congruency_EVENTS_censored.rds")),
  readRDS(here("out", "rsa", "rmatrix_vanilla_corr_shaefer400_proactive_Congruency_EVENTS_censored.rds")),
  rev.along = 0
)
names(dimnames(r.2tr1knot)) <- c("param", "run", "knot", "parcel", "subj", "session")
has.stats <- apply(r.2tr1knot, 6, function(.) !any(is.na(c(.))))
rm(r.2tr1knot)

subjs.with.stats <- names(has.stats)[has.stats]
subjs.bad <- unique(pretrial.subjsum[subj.set == "bad"]$subj)
subjs.good.with.stats <- setdiff(subjs.with.stats, subjs.bad)


## atlases ----

parcellation <- read_atlas("schaefer400")


## loop ----

n.iter <- length(sessi) * length(subjs)
pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = n.iter, clear = FALSE, width = 120
)

a <- array(
  NA,
  dim = c(
    param  = length(regressors),
    measu  = 2,
    knot   = n.knots,
    roi    = length(parcellation$key), 
    subj   = length(subjs),
    sess   = length(sessi)
  ),
  dimnames = list(
    param  = regressors,
    measu  = c("corr", "eucl"),
    knot   = paste0("knot", seq_len(n.knots)),
    roi    = parcellation$key, 
    subj   = subjs,
    sess   = sessi
  )
)


for (sess.i in seq_along(sessi)) {
  # sess.i = 1
  
  name.sess.i <- sessi[sess.i]
  
  for (subj.i in seq_along(subjs)) {
    # subj.i = 1
    
    name.subj.i <- subjs[subj.i]
    
    ## read data ----
    
    ## build paths
    
    dir.subsubj.i <- file.path(
      dir.subsubj, name.subj.i, "RESULTS", "Stroop", name.sess.i, 
      paste0(name.sess.i, "_Congruency_EVENTS_censored")
    )
    dir.ccp.hcp.i <- file.path(
      dir.ccp.hcp, name.subj.i, "SURFACE_RESULTS", "Stroop", 
      paste0(name.sess.i, "_Congruency_EVENTS_censored")
    )
    
    f.betas.subsubj <- c(
      run1_R = file.path(paste0(dir.subsubj.i, "_run1"), paste0("betas_", subjs[subj.i], "_R.func.gii")),
      run1_L = file.path(paste0(dir.subsubj.i, "_run1"), paste0("betas_", subjs[subj.i], "_L.func.gii")),
      run2_R = file.path(paste0(dir.subsubj.i, "_run2"), paste0("betas_", subjs[subj.i], "_R.func.gii")),
      run2_L = file.path(paste0(dir.subsubj.i, "_run2"), paste0("betas_", subjs[subj.i], "_L.func.gii"))
    )
    
    f.betas.ccp.hcp <- c(
      R = file.path(dir.ccp.hcp.i, paste0("STATS_", subjs[subj.i], "_REML_R.func.gii")),
      L = file.path(dir.ccp.hcp.i, paste0("STATS_", subjs[subj.i], "_REML_L.func.gii"))
    )
    
    
    ## check existence
    
    file.is.missing <- any(!file.exists(f.betas.subsubj, f.betas.ccp.hcp))
    if (file.is.missing) next
    
    
    ## read into lists
    
    betas.subsubj.l <- vector("list", 4) %>% setNames(combo_paste(c("run1", "run2"), c("L", "R")))
    betas.ccp.hcp.l <- vector("list", 2) %>% setNames(c("L", "R"))
    
    for (name.hemi.i in c("L", "R")) {
      # name.hemi.i = "L"
      
      for (name.run.i in c("run1", "run2")) {
        # name.run.i = "run1"
        
        name.run.hemi.i <- paste0(name.run.i, "_", name.hemi.i)
        
        betas.subsubj.l[[name.run.hemi.i]] <- collate_surface_params(
          f.betas.subsubj[name.run.hemi.i], 
          pattern = paste0(regressors, c("#[0-9]"), collapse = "|")
        )
        
      }
      
      betas.ccp.hcp.l[[name.hemi.i]] <- collate_surface_params(
        f.betas.ccp.hcp[name.hemi.i], 
        pattern = paste0(regressors, c("#[0-9]_Coef"), collapse = "|")
      )
      
    }
    
    ## wrangle to same format
    
    betas.subsubj <- cbind(
      betas.subsubj.l[["run1_L"]] + betas.subsubj.l[["run2_L"]],
      betas.subsubj.l[["run1_R"]] + betas.subsubj.l[["run2_R"]]
    ) / 2  ## average across run
    
    betas.ccp.hcp <- cbind(betas.ccp.hcp.l[["L"]], betas.ccp.hcp.l[["R"]])
    rownames(betas.ccp.hcp) <- gsub("_Coef", "", rownames(betas.ccp.hcp))
    
    ## check match
    
    are.matched <-
      identical(dim(betas.subsubj), dim(betas.ccp.hcp)) && 
      identical(dimnames(betas.subsubj), dimnames(betas.ccp.hcp))
    
    if (!are.matched) stop("ccphcp and subsubj arrays not matched!")
    
    # params <- rownames(betas.subsubj)
    
    ## calculate similarities ----
    
    for (parcel.i in seq_along(parcellation$key)) {
      # parcel.i = 1
      
      mask.i <- parcellation$atlas == parcel.i
      
      for (knot.i in seq_len(n.knots)) {
        # knot.i = 2
        
        name.knot.i <- paste0(knot.i - 1)
        
        for (name.regressor.i in regressors) {
          # name.regressor.i = regressors[1]
          
          name.row.i <- paste0(name.regressor.i, "#", name.knot.i)
          
          betas.subsubj.i <- betas.subsubj[name.row.i, mask.i]
          betas.ccp.hcp.i <- betas.ccp.hcp[name.row.i, mask.i]
          
          stats <- c(
            cor(betas.subsubj.i, betas.ccp.hcp.i),
            sqrt(sum((betas.subsubj.i - betas.ccp.hcp.i)^2)) / length(betas.subsubj.i)  ## scaled euclidean
          ) 
          
          a[name.regressor.i, c("corr", "eucl"), knot.i, parcel.i, subj.i, sess.i] <- stats
          
        }  ## regressor loop end
        
      }  ## knot loop end
      
    }  ## parcel loop end
    
    pb$tick()  ## progress bar
    
  }  ## subj loop end
  
}  ## session loop end


saveRDS(a, here("out", "qc", paste0("similarity-to-ccphcp-glms_shaefer400_2tr1knot.rds")))

