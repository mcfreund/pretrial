## about ----
## 
## 
## mike freund, 2020-02-19


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

## functions

whitening <- function(E, shrinkage = "ledoitwolf") {
  
  S <- cov(E)  ## sample cov
  # corrplot::corrplot(cov2cor(S), method = "color")
  H <- diag(nrow(S))  ## target to shrink towards (Ledoit-Wolf's 'F')

  if (shrinkage %in% c("lw", "ledoitwolf", "LW"))  {
    k <- tawny::shrinkage.intensity(E, H, S)
    lambda <- max(c(0, min(k / nrow(E), 1)))  ## shrinkage factor
  } else lambda <- shrinkage

  S_hat <- lambda * H + (1 - lambda) * S  ## shrunken matrix
  
  W2 <- solve(S_hat)  ## mahalanobis whitening matrix^2
  
  list(W2 = W2, lambda = lambda)
  
}

is_equal <- function(x, y, tol = .Machine$double.eps^0.5) {
  ## https://stackoverflow.com/questions/35097815/vectorized-equality-testing
  abs(x - y) < tol
}


## variables

parcellation <- read_atlas("schaefer400")

dir.ccp.hcp <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS"
dir.subsubj <- "/data/nil-external/ccp/witzermanm/AFNI_ANALYSIS_SUBSUBJECT"
dir.results <- file.path(dir.subsubj, "RESULTS_RUNWISE")

subjs <- list.dirs(dir.results, recursive = FALSE, full.names = FALSE)
# sessi <- c("baseline", "proactive", "reactive")
# sessi.short <- c("Bas", "Pro", "Rea")
sessi <- "baseline"
sessi.short <- "Bas"
n.knots <- 6
glms <- "Congruency_EVENTS_censored"
regressors <- c("PC50InCon", "biasInCon", "PC50Con", "biasCon")
measures <- c("corr", "eucl")
normalizations <- c("raw", "prw")
rsatypes <- c("vanilla", "crossval")


## initialize arrays

r.vn <- array(
  NA,
  dim = c(
    .row   = length(regressors) * 2,
    .col   = length(regressors) * 2,
    subj   = length(subjs), 
    sess   = length(sessi),
    roi    = length(parcellation$key), 
    measu  = length(measures),
    norma  = length(normalizations),
    knot   = n.knots,
    glm    = length(glms)
  ),
  dimnames = list(
    .row   = combo_paste(regressors, c("run1", "run2")),
    .col   = combo_paste(regressors, c("run1", "run2")),
    subj   = subjs, 
    sess   = sessi,
    roi    = parcellation$key, 
    measu  = measures,
    norma  = normalizations,
    knot   = paste0("knot", seq_len(n.knots)),
    glm    = glms
  )
)

r.cv <- array(
  NA,
  dim = c(
    .row   = length(regressors),
    .col   = length(regressors),
    subj   = length(subjs), 
    sess   = length(sessi),
    roi    = length(parcellation$key), 
    measu  = length(measures),
    norma  = length(normalizations),
    knot   = n.knots,
    glm    = length(glms)
  ),
  dimnames = list(
    .row   = regressors,
    .col   = regressors,
    subj   = subjs, 
    sess   = sessi,
    roi    = parcellation$key, 
    measu  = measures,
    norma  = normalizations,
    knot   = paste0("knot", seq_len(n.knots)),
    glm    = glms
  )
)

## for tallying unresponsive voxels
counts.silent <- as.data.table(expand.grid(subj = subjs, session = sessi, glm = glms, roi = parcellation$key))
counts.silent$n <- as.numeric(NA)

shrinkages <- as.data.table(expand.grid(subj = subjs, session = sessi, glm = glms, roi = parcellation$key))
shrinkages$lambda.run1 <- as.numeric(NA)
shrinkages$lambda.run2 <- as.numeric(NA)

## loop ----

n.iter <- length(sessi) * length(subjs) * length(glms) * length(parcellation$key)
pb <- progress_bar$new(
  format = " downloading [:bar] :percent eta: :eta (elapsed :elapsed)",
  total = n.iter, clear = FALSE, width = 120
)

time.begin <- Sys.time()
for (sess.i in seq_along(sessi)) {
  # sess.i = 1
  
  name.sess.i <- sessi[sess.i]
  
  for (subj.i in seq_along(subjs)) {
    ## subj.i = 1
    
    name.subj.i <- subjs[subj.i]
    
    for (name.glm.i in glms) {
      # name.glm.i <- glms[1]
      
      dirs <- file.path(dir.results, name.subj.i, "Stroop", name.sess.i, paste0(name.sess.i, "_", name.glm.i))
      
      betas <- vector("list", 4) %>% setNames(combo_paste(c("run1", "run2"), c("L", "R")))
      resid <- betas
      
      for (name.run.i in c("run1", "run2")) {
        # name.run.i = "run1"
        for (name.hemi.i in c("L", "R")) {
          # name.hemi.i = "L"
          
          f.betas <- file.path(
            paste0(dirs, "_", name.run.i), paste0("stats_", subjs[subj.i], "_", name.hemi.i, ".func.gii")
            )
          f.resid <- file.path(
            paste0(dirs, "_", name.run.i), paste0("wherr_", subjs[subj.i], "_", name.hemi.i, ".func.gii")
            )
          
          file.is.missing <- !file.exists(f.betas) | !file.exists(f.resid)
          if (file.is.missing) next  ## go to next task!
          
          name.run.hemi.i <- paste0(name.run.i, "_", name.hemi.i)
          
          betas[[name.run.hemi.i]] <- collate_surface_params(
            f.betas, 
            pattern = combopaste(regressors, c("#[0-9]*_Coef")) %>% paste0(collapse = "|")
          )
          
          resid[[name.run.hemi.i]] <- read_gifti2matrix(f.resid)
          
        }
      }
      
      betas <- abind(
        run1 = cbind(betas[["run1_L"]], betas[["run1_R"]]),
        run2 = cbind(betas[["run2_L"]], betas[["run2_L"]]),
        rev.along = 0
      )
      resid <- abind(
        run1 = cbind(resid[["run1_L"]], resid[["run1_R"]]),
        run2 = cbind(resid[["run2_L"]], resid[["run2_L"]]),
        rev.along = 0
      )
      
      
      ## estimation ----
      
      for (roi.i in seq_along(parcellation$key)) {
        # roi.i = 1
# profvis({
        ## extract roi
        
        name.roi.i <- parcellation$key[roi.i]
        betas.i <- betas[, parcellation$atlas == roi.i, ]
        resid.i <- resid[, parcellation$atlas == roi.i, ]
        
        ## remove unresponsive vertices (no variance across conditions in any run)
        
        var.vert <- apply(resid.i, c(2, 3), var)
        is.silent <- rowSums(is_equal(var.vert, 0)) > 0
        betas.i <- betas.i[, !is.silent, ]
        resid.i <- resid.i[, !is.silent, ]
        
        counts.silent[
          subj == name.subj.i & session == name.sess.i & glm == name.glm.i & roi == name.roi.i,
          "n"
          ] <- sum(is.silent)  ## tally
        
        n.vert <- ncol(betas.i)  ## number vertices
        
        for (knot.i in seq_len(n.knots)) {
          # knot.i = 1
          
          ## extract knot
          
          rows.knot.i <- grep(paste0("#", knot.i - 1, "_Coef"), rownames(betas.i))  ## minus one b/c 0-based ind.
          betas.ii <- betas.i[rows.knot.i, , ]
          
          ## reshape betas to matrix & rename/rearrange to match dims of storage array
          
          betas.ii.mat <- t(rbind(betas.ii[, , 1], betas.ii[, , 2]))
          colnames(betas.ii.mat) <- paste0(colnames(betas.ii.mat), rep(c("_run1", "_run2"), each = length(regressors)))
          colnames(betas.ii.mat) <- gsub("#.*_Coef", "", colnames(betas.ii.mat))  ## remove knot info
          betas.ii.mat <- betas.ii.mat[, rownames(r.vn)]  ## rearrange col order

          ## prewhiten patterns (bottleneck)
          
          whitened.run1 <- whitening(resid.i[, , 1])
          whitened.run2 <- whitening(resid.i[, , 2])
          
          saveRDS(
            list(run1 = whitened.run1, run2 = whitened.run2), 
            here(
              "out", "rsa", "whitening_matrices", 
              paste0("glm-", name.glm.i, "_schaefer400-", roi.i, "_subj-", name.subj.i)
              )
            )
          
          (whitened.run1$W1 + whitened.run1$W2) / 2  ## mean of inverse cov matrices
          
          W2 <- (whitening(resid.i[, , 1]) + whitening(resid.i[, , 2])) / 2  ## mean of inverse cov matrices
          W <- expm::sqrtm(W2)  ## square root (mahalanobis whitening matrix)
          betas.ii.mat.w <- W %*% betas.ii.mat
          
          # betas.ii <- plyr::aaply(betas.ii, 3, function(b) b %*% W)  ## apply
          # betas.ii <- aperm(betas.ii1, c(2, 3, 1))  ## permute array back to (condition * vertex * run )
                    
          ## vanilla RSA (includes both cross-run and within-run similarity matrices)

          ## pearson's and euclideans on 'raw' (unnormalized) patterns
          
          r.vn[, , subj.i, sess.i, roi.i, "corr", "raw", knot.i, name.glm.i] <- cor(betas.ii.mat)
          r.vn[, , subj.i, sess.i, roi.i, "eucl", "raw", knot.i, name.glm.i] <- dist2mat(betas.ii.mat) / n.vert
          
          ## now on prewhitened patterns
          
          r.vn[, , subj.i, sess.i, roi.i, "corr", "prw", knot.i, name.glm.i] <- cor(betas.ii.mat.w)
          r.vn[, , subj.i, sess.i, roi.i, "eucl", "prw", knot.i, name.glm.i] <- dist2mat(betas.ii.mat.w) / n.vert
          
          # qcor(cor(betas.ii.mat.w))
          # qcor(cor(betas.ii.mat))
          
          
          ## TODO: cross-validated RSA
          ## each measure: corrleation, euclidean
          ## each normalization: raw, prew

        }  ## knot loop end
        
        pb$tick()  ## track progress
# })
      }  ## roi loop end
      
    }  ## glm loop end
    
  }  ## subj loop end
  
}  ## end session loop

time.run <- Sys.time() - time.begin
print(time.run)

## save ----

saveRDS(r.vn, here("out", "rsa", "rmatrix_vanilla_shaefer400.rds"))
# saveRDS(r.cv, here("out", "rsa", "rmatrix_crossva_shaefer400.rds"))
