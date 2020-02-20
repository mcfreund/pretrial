## about ----
## 
## 
## mike freund, 2020-01-26

## setup ----

library(here)
library(magrittr)
library(gifti)
library(cifti)
library(abind)
library(data.table)
library(mikeutils)

## functions

whitening <- function(E, shrinkage = "ledoitwolf") {
  
  S <- cov(E)  ## sample cov
  # corrplot::corrplot(cov2cor(S), method = "color")
  H <- diag(nrow(S))  ## target to shrink towards (Ledoit-Wolf's 'F')
  
  if (shrink %in% c("lw", "ledoitwolf", "LW"))  {
    k <- tawny::shrinkage.intensity(E, H, S)
    lambda <- max(c(0, min(k / nrow(E), 1)))  ## shrinkage factor
  }
  
  S_hat <- lambda * H + (1 - lambda) * S  ## shrunken matrix
  
  solve(S_hat)  ## mahalanobis whitening matrix^2
  
}

## variables

parcellation <- read_atlas("schaefer400")

dir.ccp.hcp <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS"
dir.subsubj <- "/data/nil-external/ccp/witzermanm/AFNI_ANALYSIS_SUBSUBJECT"
dir.results <- file.path(dir.subsubj, "RESULTS_RUNWISE")

subjs <- list.dirs(dir.results, recursive = FALSE, full.names = FALSE)
# tasks <- c("Axcpt", "Cuedts", "Stern", "Stroop")
tasks <- "Stroop"
# sessi <- c("baseline", "proactive", "reactive")
# sessi.short <- c("Bas", "Pro", "Rea")
sessi <- "baseline"
sessi.short <- "Bas"
n.knots <- 8





## create "label table" (lt) for loop indices

lt <- data.table(
  variable = c("PC50Con", "PC50InCon", "biasCon", "biasInCon"),
  glm.name = c(rep( "ListLength_EVENTS_censored", 4))
)

# contrast <- c("InCon_Con_bias", "InCon_Con_PC50", "InCon_Con_PC50bias", "error_correct")

m <- list(
  Axcpt   = rbind(
    ##    ax  ay  ang  bx  by  bng  error
    cbind( 1/2,  1/2,    0, -1/2, -1/2,    0,    0),
    cbind(-1/2,  1/2,    0,  1/2, -1/2,    0,    0),
    cbind(-1/4, -1/4,  1/2, -1/4, -1/4,  1/2,    0),
    cbind(-1/4, -1/4,    0, -1/4, -1/4,    0,  1/2)
  ),
  Cuedts  = rbind(
    ##    ax  ay  ang  bx  by  bng  error
    cbind( 1/2,  1/2,    0, -1/2, -1/2,    0,    0),
    cbind(-1/2,  1/2,    0,  1/2, -1/2,    0,    0),
    cbind(-1/4, -1/4,  1/2, -1/4, -1/4,  1/2,    0),
    cbind(-1/4, -1/4,    0, -1/4, -1/4,    0,  1/2)
  ),
  Stern   = rbind(
  ),
  Stroop  = rbind(
    
  ),
)


for (ii in seq_along(m)) dimnames(m[[ii]]) <- list(contrast = contrast[[task.i]], param = )


statistics <- combo.paste(c("cor", "euc"), c("raw", "unn", "mnn"))

## lists of arrays

r <- vector("list", length(tasks)) %>% setNames(tasks)  ## d within
u <- r
for (task.i in seq_along(tasks)) {
  
  r[[task.i]] <- array(
    NA,
    dim = c(
      subj  = length(subjs), 
      sess  = length(sessi),
      roi   = length(parcellation$key), 
      cont  = length(contrast[[task.i]]),
      stat  = length(statistics),
      knot  = n.knots,
      run   = 2
    ),
    dimnames = list(
      subj  = subjs, 
      sess  = sessi,
      roi   = parcellation$key, 
      cont  = contrast[[task.i]],
      stat  = statistics,
      knot  = paste0("knot", seq_len(n.knots)),
      run   = c("run1", "run2")
    )
  )
  
  u[[task.i]] <- array(
    NA,
    dim = c(
      subj  = length(subjs), 
      sess  = length(sessi),
      roi   = length(parcellation$key), 
      cont  = length(contrast[[task.i]]),
      knot  = n.knots,
      run   = 2
    ),
    dimnames = list(
      subj  = subjs, 
      sess  = sessi,
      roi   = parcellation$key, 
      cont  = contrast[[task.i]],
      knot  = paste0("knot", seq_len(n.knots)),
      run   = c("run1", "run2")
    )
  )
  
}


## loop ----


for (task.i in seq_along(tasks)) {
  # task.i = 1
  
  name.task.i <- tasks[task.i]
  
  for (sess.i in seq_along(sessi)) {
    # sess.i = 1
    
    name.sess.i <- sessi[sess.i]
    
    for (subj.i in seq_along(subjs)) {
      ## subj.i = 1
      
      name.subj.i <- subjs[subj.i]
      
      lt.i <- lt[task == name.task.i, c("variable", "glm.name")]
      
      for (name.glm.i in unique(lt.i$glm.name)) {
        # name.glm.i = unique(lt.i$glm.name)[1]
        
        vars.i <- lt.i[glm.name == name.glm.i]$variable
        
        dirs <- combopaste(
          file.path(dir.results, name.subj.i, name.task.i),
          paste0("/", name.sess.i, "/", name.sess.i, "_", name.glm.i)
        )
        
        betas <- vector("list", 4) %>% setNames(combo_paste(c("run1", "run2"), c("L", "R")))
        resid <- betas
        
        for (name.run.i in c("run1", "run2")) {
          # name.run.i = "run1"
          for (name.hemi.i in c("L", "R")) {
            # name.hemi.i = "L"
            
            f.betas <- file.path(paste0(dirs, "_", name.run.i), paste0("stats_", subjs[subj.i], "_", name.hemi.i, ".func.gii"))
            f.resid <- file.path(paste0(dirs, "_", name.run.i), paste0("wherr_", subjs[subj.i], "_", name.hemi.i, ".func.gii"))
            
            file.is.missing <- !file.exists(f.betas) | !file.exists(f.resid)
            if (file.is.missing) next  ## go to next task!
            
            name.run.hemi.i <- paste0(name.run.i, "_", name.hemi.i)
            
            betas[[name.run.hemi.i]] <- collate_surface_params(
              f.betas, 
              pattern = combopaste(vars.i, c("#[0-9]_Coef")) %>% paste0(collapse = "|")
            )
            
            resid[[name.run.hemi.i]] <- read_gifti2matrix(f.resid)
            
          }  ## end hemi loop
        }  ## end run loop
        
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
          name.roi.i <- parcellation$key[roi.i]
          betas.i <- betas[, parcellation$atlas == roi.i, ]
          
          for (knot.i in seq_len(n.knots)) {
            # knot.i = 1
            
            ## mask (for knot and region)
            
            rows.knot.i <- grep(paste0("#", knot.i - 1, "_Coef"), rownames(betas.i))  ## minus one b/c 0-based ind.
            betas.ii <- betas.i[rows.knot.i, , ]
            
            for (contrast.i in contrast[[task.i]]) {
              
              ## check order of rownames
              
              ## get contrasts
              
              
              contrasts.i <- abind(
                run1 = m %*% betas.ii[, , "run1"],
                run2 = m %*% betas.ii[, , "run2"],
                rev.along = 0
              )
              
              ## estimate
              
              means <- rowMeans(contrasts.i)
              relia <- apply(contrasts.i, 1, function(.) cor(.)[1, 2])
              
              
              
            }
            
            
            
            
          }
          
          
          
        }
        
        
      }  ## glm loop end
      
      
      
    }  ## subj loop end
    
    
  }  ## end session loop
  
  
  ## save task array
}




## save ----

saveRDS(r, here("out", "runwise", "reliability-contrasts_runwise_schaefer400.rds"))
saveRDS(u, here("out", "runwise", "mean-contrasts_runwise_schaefer400.rds"))











































































#!/usr/bin/env Rscript

## about ----
## 
## reads in afni images (beta estimates from GLM) into a list.
## this list contains one element per subject per parcel per hemisphere.
## given a parcellation atlas or mask, similarity measures are then calculated from this list, and saved as
##  .RDS files to stroop-rsa/out/rsa/.
## additionaly saved are mean values for conducting a univariate analysis.
## 
## this script is configured to run as an executable on a *nix system.
## however, it can also be run locally (on mike's lenovo) in an interactive session.
## 
## mike freund, 2019-12-24

## TODO
## function for matching brick string

doc <- 
  "Usage:
   2_estimate_rsms.R [-a <do_atlases> -m <do_masks> -u <univariate>]
Options:
   -a Conduct analysis using atlases (Glasser's Multi Modal Parcellation, and Gordon's RSFC communities)? [default: 0]
   -m Conduct analysis using user-specified masks? [default: 0]
   -u Estimate univariate statistics? [default: 0]
 ]"

opts <- docopt::docopt(doc)

do.atlas <- as.logical(as.integer(opts$a))
do.masks <- as.logical(as.integer(opts$m))
do.univa <- as.logical(as.integer(opts$u))

## defaults for interactive use (e.g., debugging, ...):
if (interactive()) {
  do.atlas <- TRUE
  do.masks <- TRUE
  do.univa <- TRUE
}

if (!any(do.atlas, do.masks, do.univa)) stop(paste0("you must do something!"))

## setup ----

source(here::here("code", "strings.R"))
if (do.atlas) source(here::here("code", "read_atlases.R"))
if (do.masks) source(here::here("code", "read_masks.R"))

## paths, vars

dir.analysis <- here::here("glms")
glm.name <- "pro_bias_acc-only"
files.dir.analysis <- list.files(dir.analysis, pattern = "stats_var", recursive = TRUE)  ## get fit.subjs
files.dir.analysis <- files.dir.analysis[grep(glm.name, files.dir.analysis)]
fit.subjs <- unique(gsub("/results/.*", "", files.dir.analysis))

## regs will be used to pull out (via string match) the statistic from the afni brick;
## thus must have the reg.suffix
regs <- c(bias.items, "pc50_c", "pc50_i", "sustained", "transient", "nuisance")
n.regs <- length(regs)
n.subj <- length(fit.subjs)
n.bias.items <- length(bias.items)

## this variable defines the outermost loop.
## each iteration collates RSMs into a single array, and saves it as a single .rds file.
sets.of.rois <- character(0)
if (do.atlas) sets.of.rois <- c(sets.of.rois, names(atlas))
if (do.masks) sets.of.rois <- c(sets.of.rois, "masks")


## loop over sets of ROIs ----

for (set.i in sets.of.rois) {
  # set.i = "masks"
  
  ## get numbers and create storage objects
  
  if (set.i != "masks") {
    n.roi <- nrow(atlas.key[[set.i]])
    roi.names <- atlas.key[[set.i]]$roi
  } else {
    n.roi <- length(masks)
    roi.names <- names(masks)
  }
  
  rsarray.pearson <- array(  ## for representational similarity matrices
    NA,
    dim = c(n.bias.items, n.bias.items, n.subj, n.roi),
    dimnames = list(
      .row = bias.items,
      .col = bias.items,
      subj = fit.subjs,
      roi  = roi.names
    )
  )
  rsarray.euclidean <- rsarray.pearson
  
  ## for tallying voxels:
  voxels.silent <- matrix(NA, nrow = n.roi, ncol = n.subj)  ## for num unresponsive / roi
  dimnames(voxels.silent) <- list(roi = roi.names, subj = fit.subjs)
  voxels.number <- numeric(n.roi)  ## for total number of voxels / roi
  
  if (do.univa) {  ### for saving unvariate stats (means)
    roi.means <- matrix(NA, nrow = n.roi, ncol = n.regs)
    dimnames(roi.means) <- list(roi = roi.names, reg = regs)
  }
  
  ## loop over subjs and rois ----
  
  for (subj.i in seq_along(fit.subjs)) {
    # subj.i <- 1
    
    ## get afni images
    
    dir.glms <- file.path(dir.analysis, fit.subjs[subj.i], "results", glm.name)
    fname.nii <- file.path(dir.glms, paste0("stats_", fit.subjs[subj.i], ".nii.gz"))
    
    if (file.exists(fname.nii)) {
      ## dims of image.run.full [i, j, k, ???, regressor]
      image.full <- oro.nifti::readNIfTI(fname.nii, reorient = FALSE)
    } else stop("file nonexistant! ", paste0(fname.nii))
    
    ## get brick numbers for regressors
    
    brick.nums <- rep(NA, length(regs))  ## + 1 for sustained
    for (reg.i in regs) {
      # reg.i <- regs[1]
      
      image.label <- paste0(reg.i, "#0_Coef")
      
      brick.str <- mikeutils::afni("3dinfo", paste("-label2index", image.label, fname.nii))
      
      ## for error checking
      
      has.error <- grepl("error", brick.str, ignore.case = TRUE)
      if (any(has.error)) stop("error loading brick nums: ", paste0(fit.subjs[subj.i], " ", reg.i))
      
      ## to remove function call that is included in output (when is.local.session)
      
      brick.str <- brick.str[!grepl("3dinfo: AFNI version", brick.str)]
      brick.num <- as.numeric(brick.str)
      brick.nums[which(reg.i == regs)] <- brick.num
      
    }
    
    if (any(is.na(brick.nums))) stop("brick nums equal zero! ", paste0(fit.subjs[subj.i]))
    
    ## put subset of images (only the relevant regressors) into one array
    
    ## image.betas with dims [i, j, k, regs]:
    image.betas <- image.full[, , , 1, brick.nums + 1]
    rm(image.full)
    gc()  ## take out the garbage
    dimnames(image.betas) <- list(i = NULL, j = NULL, k = NULL, reg = regs)
    
    ## generate rsm for each roi
    
    for (roi.i in seq_len(n.roi)) {  ## careful: roi.i is an index for as.numeric(rois)
      # roi.i <- 1
      
      ## get and apply mask for roi.i
      
      if (set.i != "masks") mask.i <- atlas[[set.i]] == roi.i else mask.i <- masks[[roi.i]] == 1
      
      roi.betas <- apply(image.betas, "reg", function(.) .[mask.i])
      
      ## tally number of unresponsive voxels
      
      is.all.zero <- vapply(rowSums(roi.betas), function(.) isTRUE(all.equal(., 0)), logical(1))
      voxels.silent[n.roi, subj.i] <- sum(is.all.zero)
      
      ## tally number of voxels (but only do once; same for all subjects)
      
      n.voxels <- nrow(roi.betas)
      if (subj.i == 1) voxels.number[n.roi] <- n.voxels
      
      ## get rsm (pearson and euclidean)
      
      rsarray.pearson[, , subj.i, roi.i] <- cor(roi.betas[, bias.items])
      rsarray.euclidean[, , subj.i, roi.i] <- mikeutils::dist2mat(roi.betas[, bias.items]) / n.voxels
      
      ## get univariate stats (across-voxel means)
      
      if (do.univa) roi.means[roi.i, ] <- apply(roi.betas, "reg", mean)
      
    }
    
    print(paste0(subj.i, ": subj ", fit.subjs[subj.i], " done!"), quote = FALSE)
    
  }  ## end subject loop
  
  
  ## store ----
  
  ## RSA results
  
  saveRDS(
    rsarray.pearson, 
    here::here(
      "out", "rsa", "obsv", 
      paste0("rsarray_", glm.name, "_", set.i, "_pearson.rds")
    )
  )
  
  saveRDS(
    rsarray.euclidean, 
    here::here(
      "out", "rsa", "obsv", 
      paste0("rsarray_", glm.name, "_", set.i, "_euclidean.rds")
    )
  )
  
  ## univariate results
  
  if (do.univa) {
    
    saveRDS(
      roi.means, 
      here::here(
        "out", "rsa", "obsv",  ## not an RSA, but save in ./out/rsa/ for consistency...
        paste0("roi-means_", glm.name, "_", set.i, ".rds")
      )
    )
    
  }
  
  ## voxel information
  
  voxels.silent <- data.table::as.data.table(voxels.silent)
  voxels.silent$roi <- roi.names
  voxels.silent <- data.table::melt(voxels.silent, id.vars = "roi", variable.name = "subj", value.name = "n.silent")
  data.table::fwrite(
    voxels.silent,
    here::here(
      "out", "summaries", paste0("voxel-counts_unresponsive_", set.i, ".csv")
    )
  )
  
  data.table::fwrite(
    data.table::data.table(roi = roi.names, n.total = voxels.number),
    here::here(
      "out", "summaries", paste0("voxel-counts_total_", set.i, ".csv")
    )
  )
  
  
}  ## end atlas loop