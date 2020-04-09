## about ----
## reads betas from runwise GLMS, calculates parcel-wise means, writes stats as RDS files.
## 
## mike freund, 09 April 2020


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
library(foreach)
library(doParallel)

## variables

parcellation <- read_atlas("schaefer400")

dir.ccp.hcp <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS"
dir.subsubj <- "/data/nil-external/ccp/freund/sub-subj-glms/"
dir.analysis <- file.path(dir.subsubj, "runwise_old")

subjs <- list.dirs(dir.analysis, recursive = FALSE, full.names = FALSE)
sessi <- c("baseline", "proactive")
sessi.short <- c("Bas", "proactive")
glms <- c("Congruency_EVENTS_censored", "Congruency_EVENTS_tentzero01210_censored")
n.knots <- setNames(c(6, 8), glms)
regressors <- c("PC50InCon", "biasInCon", "PC50Con", "biasCon")


## initialize arrays
.r.un <- list(
  Congruency_EVENTS_censored = array(
    NA,
    dim = c(
      cond   = length(regressors),
      run    = 2,
      knot   = n.knots["Congruency_EVENTS_censored"],
      roi    = length(parcellation$key), 
      subj   = length(subjs)
    ),
    dimnames = list(
      cond   = regressors,
      run    = c("run1", "run2"),
      knot   = paste0("knot", seq_len(n.knots["Congruency_EVENTS_censored"])),
      roi    = parcellation$key, 
      subj   = subjs
    )
  ),
  Congruency_EVENTS_tentzero01210_censored = array(
    NA,
    dim = c(
      cond   = length(regressors),
      run    = 2,
      knot   = n.knots["Congruency_EVENTS_tentzero01210_censored"],
      roi    = length(parcellation$key), 
      subj   = length(subjs)
    ),
    dimnames = list(
      cond   = regressors,
      run    = c("run1", "run2"),
      knot   = paste0("knot", seq_len(n.knots["Congruency_EVENTS_tentzero01210_censored"])),
      roi    = parcellation$key, 
      subj   = subjs
    )
  )
)


## loop ----

n.iter <- length(sessi) * length(subjs) * length(glms) * length(parcellation$key)
pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = n.iter, clear = FALSE, width = 120
)

time.begin <- Sys.time()
for (sess.i in seq_along(sessi)) {
  # sess.i = 1
  
  name.sess.i <- sessi[sess.i]
  
  for (name.glm.i in glms) {
    # name.glm.i <- glms[2]
    
    r.un <- .r.un[[name.glm.i]]  ## copy empty arrays
    
    for (subj.i in seq_along(subjs)) {
      ## subj.i = 1
      
      name.subj.i <- subjs[subj.i]
      
      ## read data ----
      
      dirs <- file.path(
        dir.analysis, name.subj.i, "RESULTS", "Stroop", name.sess.i, paste0(name.sess.i, "_", name.glm.i)
        )
      
      betas <- vector("list", 4) %>% setNames(combo_paste(c("run1", "run2"), c("L", "R")))

      f.betas <- c(
        run1_R = file.path(paste0(dirs, "_run1"), paste0("betas_", subjs[subj.i], "_R.func.gii")),
        run1_L = file.path(paste0(dirs, "_run1"), paste0("betas_", subjs[subj.i], "_L.func.gii")),
        run2_R = file.path(paste0(dirs, "_run2"), paste0("betas_", subjs[subj.i], "_R.func.gii")),
        run2_L = file.path(paste0(dirs, "_run2"), paste0("betas_", subjs[subj.i], "_L.func.gii"))
      )

      file.is.missing <- any(!file.exists(f.betas))
      if (file.is.missing) next
      
      for (name.run.i in c("run1", "run2")) {
        # name.run.i = "run1"
        for (name.hemi.i in c("L", "R")) {
          # name.hemi.i = "L"
          
          name.run.hemi.i <- paste0(name.run.i, "_", name.hemi.i)
          
          betas[[name.run.hemi.i]] <- collate_surface_params(
            f.betas[name.run.hemi.i], 
            pattern = paste0(regressors, c("#[0-9]"), collapse = "|")
          )  ## bottleneck (read betas, not stats, for speed (b/c only need betas))
          
        }
      }
      
      ## wrangle

      betas <- abind(
        run1 = cbind(betas[["run1_L"]], betas[["run1_R"]]),
        run2 = cbind(betas[["run2_L"]], betas[["run2_L"]]),
        rev.along = 0
      )

      ## estimation ----
      
      for (roi.i in seq_along(parcellation$key)) {
        # roi.i = 1
        
        name.roi.i <- parcellation$key[roi.i]
        betas.i <- betas[, parcellation$atlas == roi.i, ]

        for (knot.i in seq_len(n.knots[name.glm.i])) {
          # knot.i = 1
          
          ## extract knot
          
          rows.knot.i <- grep(paste0("#", knot.i - 1), rownames(betas.i))  ## minus one b/c 0-based ind.
          betas.ii <- betas.i[rows.knot.i, , ]
          rownames(betas.ii) <- gsub("#[0-9]", "", rownames(betas.ii))  ## remove knot info
          betas.ii <- betas.ii[rownames(r.un), , ]  ## rearrange row order
          
          r.un[, , knot.i, roi.i, subj.i] <- apply(betas.ii, 3, rowMeans)  ## calculate means
          
          
        }  ## knot loop end
        
        pb$tick()  ## progress bar
        
      }  ## roi loop end
        
    }  ## subj loop end
    
    
    ## save ----    
    
    saveRDS(r.un, here("out", "rsa", paste0("means_shaefer400_", name.sess.i, "_", name.glm.i, ".rds")))
      
  }  ## glm loop end
  
}  ## session loop end

time.run <- Sys.time() - time.begin
print(time.run)



