## about ----
## 
## calculates condition-by-condition similarity matrices (per parcel, subject, knot) using betas from runwise GLMs.
## - two 'types' of similarity matrices are estimated: "vanilla" and "cross-validated".
## - for each 'type', three similarity measurments are used: linear correlation ("corr"), euclidean distance ("eucl"), 
##   and normalized euclidean distance ("neuc"; i.e., euclidean's on betas z-scored within parcel)
## - each type*measure combination is performed on non-prewhitened ("raw") and prewhitened patterns.
## 
## these matrices are collated into array objects and saved as .RDS for easy I/O.
## 
## additionally, this script reads/writes from ./out/rsa/whitening_matrices.
## in this directory, vertex-by-vertex whitening matrices are saved as .rds files (one for each parcel).
## computing these matrices is necessary for prewhitening, but is a computationally expensive operation.
## prior to estimating them, however, the directory is checked for their existance for a given subject.
## if all Nparcel matrices exist for a given subject, the matrices are read instead of computed.
## else, they are computed and written.
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
library(foreach)
library(doParallel)

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

contrast_matrix <- function(n, condition.names) {
  # n <- 10
  # condition.names <- letters[1:n]
  
  if (n < 2) stop("you need more than 1 condition")
  
  W <- matrix(0, nrow = n^2, ncol = n)
  
  if (missing(condition.names)) {
    dimnames(W) <- list(contrast = NULL, condition = NULL)
  } else {
    dimnames(W) <- list(
      contrast = paste0(rep(condition.names, each = n), "_", rep(condition.names, n)),
      condition = condition.names
    )
  }

  for (condition.i in seq_len(n)) {
    # condition.i = 1
    
    row.beg <- (condition.i - 1) * n + 1
    row.end <- (condition.i - 1) * n + n
    W.i <- W[row.beg:row.end, ]  ## square matrix
    
    W.i[, condition.i] <- 1
    diag(W.i) <- diag(W.i) - 1
    
    W[row.beg:row.end, ] <- W.i
    
  }
  
  W
  
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
glms <- c("Congruency_EVENTS_tentzero01210_censored", "Congruency_EVENTS_censored")
n.knots <- setNames(c(6, 8), glms)
regressors <- c("PC50InCon", "biasInCon", "PC50Con", "biasCon")
measures <- c("corr", "eucl", "neuc")
normalizations <- c("raw", "prw")
# rsatypes <- c("vanilla", "crossval")


## for tallying unresponsive voxels
# counts.silent <- as.data.table(expand.grid(subj = subjs, session = sessi, glm = glms, roi = parcellation$key))
# counts.silent$n <- as.numeric(NA)

## for saving estimated shrinkage factors
# shrinkages <- as.data.table(expand.grid(subj = subjs, session = sessi, glm = glms, roi = parcellation$key))
# shrinkages$lambda.run1 <- as.numeric(NA)
# shrinkages$lambda.run2 <- as.numeric(NA)


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

  for (name.glm.i in glms) {
    # name.glm.i <- glms[1]
    
    ## initialize arrays
    
    r.vn <- array(
      NA,
      dim = c(
        .row   = length(regressors) * 2,
        .col   = length(regressors) * 2,
        measu  = length(measures),
        norma  = length(normalizations),
        knot   = n.knots[name.glm.i],
        roi    = length(parcellation$key), 
        subj   = length(subjs)
      ),
      dimnames = list(
        .row   = combo_paste(regressors, c("run1", "run2")),
        .col   = combo_paste(regressors, c("run1", "run2")),
        measu  = measures,
        norma  = normalizations,
        knot   = paste0("knot", seq_len(n.knots[name.glm.i])),
        roi    = parcellation$key, 
        subj   = subjs
      )
    )
    
    r.cv <- array(
      NA,
      dim = c(
        .row   = length(regressors),
        .col   = length(regressors),
        measu  = length(measures),
        norma  = length(normalizations),
        knot   = n.knots[name.glm.i],
        roi    = length(parcellation$key), 
        subj   = length(subjs)
      ),
      dimnames = list(
        .row   = regressors,
        .col   = regressors,
        measu  = measures,
        norma  = normalizations,
        knot   = paste0("knot", seq_len(n.knots[name.glm.i])),
        roi    = parcellation$key, 
        subj   = subjs
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
              paste0("glm-", name.glm.i, "_schaefer400-", 1:400, "_subj-", name.subj.i, ".RDS")
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
        
        ## compute similarities

        cl <- makeCluster(n.cores - 1, type = "FORK")
        registerDoParallel(cl)
        
        results.subj.i <- foreach(roi.i = seq_along(parcellation$key)) %dopar% {
          # roi.i = 1
          
          name.roi.i <- parcellation$key[roi.i]
          betas.i <- betas[, parcellation$atlas == roi.i, ]
          if (need.to.calc.invcov) resid.i <- resid[, parcellation$atlas == roi.i, ]
          
          ## remove unresponsive vertices (no variance across conditions in any run)
          ## use betas instead of resids, bc resids might not be in environment
          
          # var.vert <- apply(resid.i, c(2, 3), var)
          # is.silent <- rowSums(is_equal(var.vert, 0)) > 0
          is.silent <- rowSums(colSums(is_equal(betas.i, 0))) > 0
          betas.i <- betas.i[, !is.silent, ]
          if (need.to.calc.invcov) resid.i <- resid.i[, !is.silent, ]
          
          # counts.silent[
          #   subj == name.subj.i & session == name.sess.i & glm == name.glm.i & roi == name.roi.i,
          #   "n"
          #   ] <- sum(is.silent)  ## tally
          
          n.vert <- ncol(betas.i)  ## number vertices
          
          r.vn.subj.i.roi.i <- r.vn[, , , , , 1, 1]  ## initialize array slice
          r.vn.subj.i.roi.i[] <- NA
          
          r.cv.subj.i.roi.i <- r.cv[, , , , , 1, 1]  ## initialize array slice
          r.cv.subj.i.roi.i[] <- NA
          
          
          ## get prewhitening matrices (a bottleneck, so save/read results)
          
          fname.mahal <- here(
            "out", "rsa", "whitening_matrices", 
            paste0("glm-", name.glm.i, "_schaefer400-", roi.i, "_subj-", name.subj.i, ".RDS")
          )
          
          if (file.exists(fname.mahal)) {
            
            whitened <- readRDS(fname.mahal)
            
          } else {
            
            whitened <- list(
              run1 = whitening(resid.i[, , 1], shrinkage = 0.4),
              run2 = whitening(resid.i[, , 2], shrinkage = 0.4)
            )
            
            saveRDS(whitened, fname.mahal)
            
          }

          for (knot.i in seq_len(n.knots[name.glm.i])) {
            # knot.i = 1
            
            ## extract knot
            
            rows.knot.i <- grep(paste0("#", knot.i - 1), rownames(betas.i))  ## minus one b/c 0-based ind.
            betas.ii <- betas.i[rows.knot.i, , ]
            
            ## reshape betas to matrix & rename/rearrange to match dims of storage array
            
            betas.ii.mat <- t(rbind(betas.ii[, , 1], betas.ii[, , 2]))
            colnames(betas.ii.mat) <- paste0(colnames(betas.ii.mat), rep(c("_run1", "_run2"), each = length(regressors)))
            colnames(betas.ii.mat) <- gsub("#[0-9]", "", colnames(betas.ii.mat))  ## remove knot info
            betas.ii.mat <- betas.ii.mat[, rownames(r.vn)]  ## rearrange col order
            
            W2 <- (whitened$run1$W2 + whitened$run1$W2) / 2  ## mean of inverse cov matrices
            W <- expm::sqrtm(W2)  ## square root (mahalanobis whitening matrix)
            betas.ii.mat.w <- W %*% betas.ii.mat
            
            # betas.ii <- plyr::aaply(betas.ii, 3, function(b) b %*% W)  ## apply
            # betas.ii <- aperm(betas.ii1, c(2, 3, 1))  ## permute array back to (condition * vertex * run )
            
            ## vanilla RSA (includes both cross-run and within-run similarity matrices)
            
            r.vn.subj.i.roi.i[, , "corr", "raw", knot.i] <- cor(betas.ii.mat)  ## pearson
            r.vn.subj.i.roi.i[, , "corr", "prw", knot.i] <- cor(betas.ii.mat.w)
            
            r.vn.subj.i.roi.i[, , "eucl", "raw", knot.i] <- dist2mat(betas.ii.mat) / n.vert  ## euclidean
            r.vn.subj.i.roi.i[, , "eucl", "prw", knot.i] <- dist2mat(betas.ii.mat.w) / n.vert
            
            r.vn.subj.i.roi.i[, , "neuc", "raw", knot.i] <- dist2mat(scale(betas.ii.mat)) / n.vert  ## norm. euclidean
            r.vn.subj.i.roi.i[, , "neuc", "prw", knot.i] <- dist2mat(scale(betas.ii.mat.w)) / n.vert
            
            ## cross-validated RSA
            
            cmat <- contrast_matrix(length(regressors), regressors)

            B <- t(betas.ii.mat)
            B.run1 <- B[grep("run1", rownames(B)), ]
            B.run2 <- B[grep("run2", rownames(B)), ]
            
            
            B.run1.c <- sweep(B.run1, 2, colMeans(B.run1))  ## center
            B.run2.c <- sweep(B.run2, 2, colMeans(B.run2))
            ## cross-validated covariance matrix (save this, not the correlation matrix!):
            r.cv.subj.i.roi.i[, , "corr", "raw", knot.i] <- B.run1.c %*% t(B.run2.c)
            r.cv.subj.i.roi.i[, , "corr", "prw", knot.i] <- B.run1.c %*% W2 %*% t(B.run2.c)
            
            cveucl.v <- colSums(t(cmat %*% B.run1) * t(B.run2) %*% t(cmat))  ## euclidean
            cvmaha.v <- colSums(t(cmat %*% B.run1 %*% W2) * t(B.run2) %*% t(cmat))  ## mahalanobis
            r.cv.subj.i.roi.i[, , "eucl", "raw", knot.i] <- matrix(cveucl.v, ncol = length(regressors))
            r.cv.subj.i.roi.i[, , "eucl", "prw", knot.i] <- matrix(cvmaha.v, ncol = length(regressors))

            cvneuc.v <- colSums(t(cmat %*% scale(B.run1)) * t(scale(B.run2)) %*% t(cmat))  ## norm. euclidean
            cvnmah.v <- colSums(t(cmat %*% scale(B.run1) %*% W2) * t(scale(B.run2)) %*% t(cmat))  ## norm. mahalanobis
            r.cv.subj.i.roi.i[, , "neuc", "raw", knot.i] <- matrix(cvneuc.v, ncol = length(regressors))
            r.cv.subj.i.roi.i[, , "neuc", "prw", knot.i] <- matrix(cvnmah.v, ncol = length(regressors))

          }  ## knot loop end
          
          list(vanilla = r.vn.subj.i.roi.i, crossvalidated = r.cv.subj.i.roi.i)  ## return results from cluster
          
        }  ## roi loop end
        
        stopCluster(cl)
        
        pb$tick()  ## progress bar
        
        ## extract results and store in arrays
        
        for (roi.i in seq_along(results.subj.i)) {
          r.vn[, , , , , roi.i, subj.i] <- results.subj.i[[roi.i]]$vanilla
          r.cv[, , , , , roi.i, subj.i] <- results.subj.i[[roi.i]]$crossvalidated
        }
        
        rm(results.subj.i, betas)
        if (need.to.calc.invcov) rm(resid)
        gc()
        
    }  ## subj loop end
    
    
    ## save ----    
    
    for (measure.i in measures) {
      
      saveRDS(
        r.vn[, , measure.i, , , , ], 
        here("out", "rsa", paste0("rmatrix_vanilla_", measure.i, "_shaefer400_", name.sess.i, "_", name.glm.i, ".rds"))
      )
      
      saveRDS(
        r.cv[, , measure.i, , , , ], 
        here("out", "rsa", paste0("rmatrix_crossva_", measure.i, "_shaefer400_", name.sess.i, "_", name.glm.i, ".rds"))
      )

    }
    
    
  }  ## glm loop end
  
}  ## session loop end

time.run <- Sys.time() - time.begin
print(time.run)

