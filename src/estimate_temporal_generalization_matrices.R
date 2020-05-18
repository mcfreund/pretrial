## about ----
## 
## calculates cross-condition temporal generalization matrices (per parcel, subject) using betas from runwise GLMs.
## 
## - each TR (1 thru 8) and pc condition (pc50, bias) serves as a "training" dataset.
##   training datasets are always run 1.
## - likewise, each TR and pc condition serves as a "test" dataset.
##   test datasets are always run 2.
## - two types of decoders are trained to distinguish I--C: discriminant and template.
##   * discriminant involves projecting test datapoints onto linear discriminant function, and calculating distance 
##     from LDA hyperplane.
##   * template involves correlating each test pattern with two template patterns (I, C), then calculating the
##     difference between correlations.
## - because each TR serves as a training and test dataset, this process yields a TR-byTR confusion matrix.
##   this is called the **temporal generalization** matrix because it indicates the degree to which coding for I--C
##   generalizes across timepoints.
## - because each condition serves as a training and test dataset, two types of temporal generalization matrices are 
##   generated: within-condition matrices (bias--bias, pc50--pc50) and between-condition matrices 
##   (bias--pc50, pc50--bias).
##   the latter indicate the degree to which coding generalizes across both time and conditions.
##   these matrices are referred to as the **cross-condition temporal generalization matrices**, or cctg.
## 
## these matrices are collated into array objects and saved as .RDS for easy I/O.
## 
## additionally, this script reads/writes from ./out/rsa/whitening_matrices.
## in this directory, vertex-by-vertex whitening matrices are saved as .rds files (one for each parcel).
## computing these matrices is necessary for prewhitening, but is a computationally expensive operation.
## prior to estimating them, however, the directory is checked for their existance for a given subject.
## if all Nparcel matrices exist for a given subject, the matrices are read instead of computed.
## else, they are computed and written.

## mike freund, 2020-05-18


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

## functions

whitecor <- function(x, y, W) {
  
  xc <- x - mean(x)
  yc <- y - mean(y)
  
  xcs <- xc / sqrt(sum(xc^2))
  ycs <- yc / sqrt(sum(yc^2))
  
  c(xcs %*% W %*% ycs)
  
}

## variables

n.cores <- detectCores()

parcellation <- read_atlas("schaefer400")

dir.ccp.hcp <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS"
dir.subsubj <- "/data/nil-external/ccp/freund/sub-subj-glms/"
dir.analysis <- file.path(dir.subsubj, "runwise_old")

subjs <- list.dirs(dir.analysis, recursive = FALSE, full.names = FALSE)
sessi <- c("baseline", "proactive")
sessi.short <- c("Bas", "Pro")
# glms <- c("Congruency_EVENTS_tentzero01210_censored", "Congruency_EVENTS_censored")
glms <- c("Congruency_EVENTS_tentzero01210_censored")
name.glm.i <- "Congruency_EVENTS_tentzero01210_censored"
# n.knots <- setNames(c(8, 6), glms)  ## SAME ORDER AS GLMS OBJECT ABOVE!
n.knots <- setNames(c(8), glms)  ## SAME ORDER AS GLMS OBJECT ABOVE!
regressors <- c("PC50InCon", "biasInCon", "PC50Con", "biasCon")
# measures <- c("corr", "eucl", "neuc")
# normalizations <- c("raw", "prw")
normalizations <- c("prw")
# rsatypes <- c("vanilla", "crossva")
projtypes <- c("discriminant", "template")

## for tallying unresponsive voxels:
# counts.silent <- as.data.table(expand.grid(subj = subjs, session = sessi, glm = glms, roi = parcellation$key))
# counts.silent$n <- as.numeric(NA)


## loop ----

n.iter <- length(sessi) * length(subjs) * length(glms)
pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = n.iter, clear = FALSE, width = 120
)

time.begin <- Sys.time()
for (sess.i in seq_along(sessi)) {
  # sess.i = 1
  
  name.sess.i <- sessi[sess.i]
  
  ## initialize arrays
  
  a <- array(
    NA,
    dim = c(
      train.tr   = n.knots[name.glm.i],
      test.tr    = n.knots[name.glm.i],
      train.cond = 2,
      test.cond  = 2,
      projtypes  = length(projtypes),
      roi        = length(parcellation$key), 
      subj       = length(subjs)
    ),
    dimnames = list(
      train.tr = NULL,
      test.tr  = NULL,
      train.cond = c("PC50", "bias"),
      test.cond  = c("PC50", "bias"),
      projtypes  = projtypes,
      roi      = parcellation$key, 
      subj     = subjs
    )
  )
  
  for (subj.i in seq_along(subjs)) {
    ## subj.i = 1
    
    name.subj.i <- subjs[subj.i]
    
    ## read data ----
    
    dirs <- file.path(
      dir.analysis, name.subj.i, "RESULTS", "Stroop", name.sess.i, paste0(name.sess.i, "_", name.glm.i)
    )
    
    betas <- vector("list", 4) %>% setNames(combo_paste(c("run1", "run2"), c("L", "R")))
    resid <- betas
    
    f.betas <- c(
      run1_R = file.path(paste0(dirs, "_run1"), paste0("betas_", subjs[subj.i], "_R.func.gii")),
      run1_L = file.path(paste0(dirs, "_run1"), paste0("betas_", subjs[subj.i], "_L.func.gii")),
      run2_R = file.path(paste0(dirs, "_run2"), paste0("betas_", subjs[subj.i], "_R.func.gii")),
      run2_L = file.path(paste0(dirs, "_run2"), paste0("betas_", subjs[subj.i], "_L.func.gii"))
    )
    f.resid <- c(
      run1_R = file.path(paste0(dirs, "_run1"), paste0("wherr_", subjs[subj.i], "_R.func.gii")),
      run1_L = file.path(paste0(dirs, "_run1"), paste0("wherr_", subjs[subj.i], "_L.func.gii")),
      run2_R = file.path(paste0(dirs, "_run2"), paste0("wherr_", subjs[subj.i], "_R.func.gii")),
      run2_L = file.path(paste0(dirs, "_run2"), paste0("wherr_", subjs[subj.i], "_L.func.gii"))
    )
    
    file.is.missing <- any(!file.exists(f.betas, f.resid))
    if (file.is.missing) next
    
    need.to.calc.invcov <- !all(
      file.exists(
        here(
          "out", "rsa", "whitening_matrices", 
          paste0(
            "glm-", name.glm.i, "_session-", name.sess.i, "_schaefer400-", 1:400, "_subj-",
            name.subj.i, ".RDS"
          )
        )
      )
    )
    
    for (name.run.i in c("run1", "run2")) {
      # name.run.i = "run1"
      for (name.hemi.i in c("L", "R")) {
        # name.hemi.i = "L"
        
        name.run.hemi.i <- paste0(name.run.i, "_", name.hemi.i)
        
        betas[[name.run.hemi.i]] <- collate_surface_params(
          f.betas[name.run.hemi.i], 
          pattern = paste0(regressors, c("#[0-9]"), collapse = "|")
        )  ## bottleneck (read betas, not stats, for speed (b/c only need betas))
        
        ## read residuals only if inverse covariance matrices do not exist (bottleneck)
        if (need.to.calc.invcov)  resid[[name.run.hemi.i]] <- read_gifti2matrix(f.resid[name.run.hemi.i])
        
      }
    }
    
    ## wrangle
    
    betas <- abind(
      run1 = cbind(betas[["run1_L"]], betas[["run1_R"]]),
      run2 = cbind(betas[["run2_L"]], betas[["run2_R"]]),
      rev.along = 0
    )
    if (need.to.calc.invcov)  {
      resid <- abind(
        run1 = cbind(resid[["run1_L"]], resid[["run1_R"]]),
        run2 = cbind(resid[["run2_L"]], resid[["run2_R"]]),
        rev.along = 0
      )
    }
    
    
    ## estimation ----
    
    cl <- makeCluster(n.cores - 1, type = "FORK")
    registerDoParallel(cl)
    
    results.subj.i <- foreach(roi.i = seq_along(parcellation$key)) %dopar% {
      # for (roi.i in seq_along(parcellation$key)) {  ## for debugging.
      # roi.i = 1
      
      name.roi.i <- parcellation$key[roi.i]
      betas.i <- betas[, parcellation$atlas == roi.i, ]
      if (need.to.calc.invcov) resid.i <- resid[, parcellation$atlas == roi.i, ]
      
      ## remove unresponsive vertices (no variance across conditions in any run)
      
      ## NB: use betas instead of resids, bc resids might not be in environment:
      # var.vert <- apply(resid.i, c(2, 3), var)
      # is.silent <- rowSums(is_equal(var.vert, 0)) > 0
      is.silent <- rowSums(colSums(is_equal(betas.i, 0))) > 0
      betas.i <- betas.i[, !is.silent, ]
      if (need.to.calc.invcov) resid.i <- resid.i[, !is.silent, ]
      
      n.vert <- ncol(betas.i)  ## number responsive vertices
      
      ## initialize array slices
      
      a.subj.i.roi.i <- a[, , , , , 1, 1]

      ## get prewhitening matrices (a bottleneck, so save/read results)
      
      fname.mahal <- here(
        "out", "rsa", "whitening_matrices", 
        paste0("glm-", name.glm.i, "_session-", name.sess.i, "_schaefer400-", roi.i, "_subj-", name.subj.i, ".RDS")
      )
      
      if (file.exists(fname.mahal)) {
        
        whitened <- readRDS(fname.mahal)
        
        if (!identical(whitened$inds, which(!is.silent))) 
          stop("non-conformable arrays: betas and whitening matrices (check 1)")
        
      } else {
        
        whitened <- list(
          run1 = whitening(resid.i[, , 1], shrinkage = 0.4),
          run2 = whitening(resid.i[, , 2], shrinkage = 0.4),
          inds = which(!is.silent)
        )
        
        saveRDS(whitened, fname.mahal)
        
      }
      
      dims.are.good <- nrow(whitened$run1$W2) == nrow(whitened$run2$W2) & length(!is.silent)
      if (!dims.are.good) stop("non-conformable arrays: betas and whitening matrices (check 2)")
      
      ## loop over conditions and TRs
      
      for (train.cond.i in dimnames(a.subj.i.roi.i)$train.cond) {
        # train.cond.i = "PC50"
        
        for (test.cond.i in dimnames(a.subj.i.roi.i)$test.cond) {
          # test.cond.i = "PC50"
          
          for (train.tr.i in seq_len(n.knots)) {
            # train.tr.i = 1
            
            rows.train <- c(
              paste0(train.cond.i, "InCon#", train.tr.i - 1),
              paste0(train.cond.i, "Con#", train.tr.i - 1)
            )
            
            betas.ii.train <- betas.i[rows.train, , ]
            
            W2 <- (whitened$run1$W2 + whitened$run1$W2) / 2  ## mean of inverse cov matrices
            
            m_dif <- betas.ii.train[1, , ] - betas.ii.train[2, , ]  ## I - C difference vector per run
            
            m_bar <- (betas.ii.train[1, , ] + betas.ii.train[2, , ]) / 2  ## mean pattern vector per run
            
            for (test.tr.i in seq_len(n.knots)) {
              # test.tr.i = 1
              
              rows.test <- c(
                paste0(test.cond.i, "InCon#", test.tr.i - 1),
                paste0(test.cond.i, "Con#", test.tr.i - 1)
              )
              
              betas.ii.test <- betas.i[rows.test, , ]
              
              ## NB: "train" dataset is always run 1, "test" is always run 2
              
              ## projection onto discriminant method
              
              b_len <- sqrt(colSums(m_dif^2))
              
              proj.incon1 <- m_dif[, "run2"] %*% W2 %*% (betas.ii.test[1, , "run1"] - m_bar[, "run2"]) / b_len["run2"] 
              proj.incon2 <- m_dif[, "run1"] %*% W2 %*% (betas.ii.test[1, , "run2"] - m_bar[, "run1"]) / b_len["run1"]
              
              proj.congr1 <- m_dif[, "run2"] %*% W2 %*% (betas.ii.test[2, , "run1"] - m_bar[, "run2"]) / b_len["run2"]
              proj.congr2 <- m_dif[, "run1"] %*% W2 %*% (betas.ii.test[2, , "run2"] - m_bar[, "run1"]) / b_len["run1"]
              
              a.subj.i.roi.i[train.tr.i, test.tr.i, train.cond.i, test.cond.i, "discriminant"] <- 
                mean(c(proj.incon1, proj.incon2, -proj.congr1, -proj.congr2))
              
              ## correlation with template method
              
              cor.incon1.incon2 <- whitecor(betas.ii.train[1, , "run1"], betas.ii.test[1, , "run2"], W2)
              cor.incon1.congr2 <- whitecor(betas.ii.train[1, , "run1"], betas.ii.test[2, , "run2"], W2)
              
              cor.congr1.incon2 <- whitecor(betas.ii.train[2, , "run1"], betas.ii.test[1, , "run2"], W2)
              cor.congr1.congr2 <- whitecor(betas.ii.train[2, , "run1"], betas.ii.test[2, , "run2"], W2)
              
              a.subj.i.roi.i[train.tr.i, test.tr.i, train.cond.i, test.cond.i, "template"] <- 
                mean(
                  c(
                    atanh(cor.incon1.incon2) - atanh(cor.congr1.incon2),
                    atanh(cor.congr1.congr2) - atanh(cor.incon1.congr2)
                  )
                )

            }  ## end 'test' loop
            
          }  ##
            
        }

      }
      
      a.subj.i.roi.i  ## return results from cluster
        
      }  ## roi loop end
      
      stopCluster(cl)
      
      pb$tick()  ## progress bar
      
      ## extract results and store in arrays
      
      
      for (roi.i in seq_along(results.subj.i)) a[, , , , , roi.i, subj.i] <- results.subj.i[[roi.i]]
      
      ## take out trash (overkill, but just in case...)
      
      rm(results.subj.i, betas, whitened, a.subj.i.roi.i)
      if (need.to.calc.invcov) rm(resid)
      gc()
      
    }  ## subj loop end
    
    
    ## save ----    
    
    for (projtype in projtypes) {
      # projtype = "discriminant"
      
      saveRDS(
        a[, , , , projtype, , ], 
        here(
          "out", "cctg", 
          paste0("matrix_", projtype, "_prw_shaefer400_", name.sess.i, "_", name.glm.i, ".rds")
          )
      )
      
    }
    
  
}  ## session loop end


time.run <- Sys.time() - time.begin
print(time.run)

