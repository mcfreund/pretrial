## about ----
## 
## fits general linear models on each subject's RSMs from each parcel then writes the results to files.
## 
## mike freund, created 2020-05-01

## TODO
## loop over session, glm, measure, subj
## lapply over knot, norma, roi
## 
## ANALYSES:
##  - vanilla: for both between and within RSA

## setup ----

library(here)
library(mikeutils)
library(dplyr)
library(purrr)
library(magrittr)
library(data.table)

source(here("src", "setup.R"))

## read models

models <- lapply(
  modname,
  function(.) data.table::fread(here("out", "rsa", "models", paste0("rsmodel_", ., ".csv")), data.table = FALSE)
)
names(models) <- modname


## create design matrices

models.d <- lapply(models, reshape2::melt)
Xcross.d <- models.d[grep("cross", modname)] %>% 
  bind_rows(.id = "model") %>% 
  tidyr::pivot_wider(values_from = "value", names_from = "model")
names(Xcross.d) <- gsub("cross_", "", names(Xcross.d))
Xcross <- Xcross.d %>% select(-V1, -variable) %>% as.matrix
Xcross <- cbind(intercept = 1, Xcross)

Xvanil.d <- models.d[grep("vanil", modname)] %>% 
  bind_rows(.id = "model") %>% 
  tidyr::pivot_wider(values_from = "value", names_from = "model")
names(Xvanil.d) <- gsub("vanil_", "", names(Xvanil.d))
Xvanil <- Xvanil.d %>% select(-V1, -variable) %>% as.matrix
Xvanil <- cbind(intercept = 1, Xvanil)

subjs <- subjs.analysis

## loop ----

## normal: prew, raw
## method: vanil, cross 
## roi
## model


for (session.i in seq_along(session)) {
  # session.i = 1
  
  name.session.i <- session[session.i]
  
  ## get knots
  
  for (glm.i in seq_along(glmname)) {
    # glm.i = 1
  
    name.glm.i <- glmname[glm.i]
    
    for (measure.i in seq_along(measure)) {
      # measure.i = 1
      
      name.measure.i <- measure[measure.i]
      
      ## read data
      
      vanilla <- readRDS(
        here(
        "out", "rsa", "observed",
        paste0("rmatrix_vanilla_", name.measure.i, "_shaefer400_", name.session.i, "_", name.glm.i, ".rds"))
      )
      
      crossva <- readRDS(
        here(
          "out", "rsa", "observed",
          paste0("rmatrix_crossva_", name.measure.i, "_shaefer400_", name.session.i, "_", name.glm.i, ".rds"))
      )
      
      ## transformations
      
      if (name.measure.i %in% c("neuc", "eucl")) {
        
        crossva <- -crossva
        vanilla <- -vanilla
        
      } else if (name.measure.i %in% "corr") {
        
        vanilla <- atanh(vanilla)
        crossva <- atanh(crossva)
        
      }
      
      for (subj.i in seq_along(subjs)) {
        # subj.i = 1
        
        name.subj.i <- subjs[subj.i]
        
        ## vanilla: between-run
        
        a.vanil <- reshape2::melt(vanilla[, , , , , subj.i])
        g.vanil <- interaction(a.vanil$norma, a.vanil$knot, a.vanil$roi)
        l.vanil <- split(a.vanil, g.vanil)
        f.vanil <- lapply(l.vanil, function(., x) coef(.lm.fit(x, .$value)), x = Xvanil)
        
        ## vanilla: within-run
        
        ## crossvalidated
        
        a.cross <- reshape2::melt(crossva[, , , , , subj.i])
        g.cross <- interaction(a.cross$norma, a.cross$knot, a.cross$roi)
        l.cross <- split(a.cross, g.cross)
        f.cross <- lapply(l.cross, function(., x) coef(.lm.fit(x, .$value)), x = Xcross)
        
        
        ## collate
        
        
        
      }
      
      
    }


  }
  
}


## loop over models
## loop over ROIs
## loop over subjs
## fit glms
## format
## write







## models will be variants upon these main effects:
variables <- as.matrix(rsv.models.ltri[sapply(rsv.models.ltri, is.numeric)])
variables <- cbind(variables, conclust = variables[, "incongruency"] | variables[, "congruency"])

variables.std <- scale(variables)

## models will contain variables with these names:
structures <- list(
  tdic    = c("target", "distractor", "incongruency", "congruency"),
  tdi     = c("target", "distractor", "incongruency"),
  tdclust = c("target", "distractor", "conclust"),
  all     = c("target", "cielab", "tphono", "distractor", "silhou", "orthog", "dphono", "incongruency", "congruency"),
  continu = c("target", "cielab", "tphono", "distractor", "silhou", "dphono", "incongruency", "congruency")
)

m <- list(
  tdic    = cbind(b0 = 1, variables[, structures$tdic]),
  tdi     = cbind(b0 = 1, variables[, structures$tdi]),
  tdclust = cbind(b0 = 1, variables[, structures$tdclust]),
  all     = cbind(b0 = 1, variables[, structures$all]),  ## continuous and categorical
  continu = cbind(b0 = 1, variables[, structures$continu])  ## perceptual and phonological for target and distractor
)
m.std <- lapply(m, function(x) scale(x[, -1]))  ## drop the intercept

## funs

get.partr <- function(x, colname, ...) {
  rsv <- x[, colname]
  d <- cbind(rsv, X[, -1])
  psych::partial.r(d, ...)[1, -1]
}


## loop over sets of ROIs ----

for (set.i in sets.of.rois) {
  # set.i = "mmp"
  
  ## read observed similarity matrices (arrays)
  
  rsarray.rank <- readRDS(
    here(
      "out", "rsa", "obsv",
      paste0("rsarray_", glmname, "_", set.i, "_pearson_residual-rank.rds")
    )
  )
  rsarray.line <- readRDS(
    here(
      "out", "rsa", "obsv",
      paste0("rsarray_", glmname, "_", set.i, "_pearson_residual-linear.rds")
    )
  )
  
  
  ## prepare similarity matrices for regression ----
  
  ## check if rows and col names are equal (should be, but just to be sure...)
  
  are.rowcol.equal <- isTRUE(
    all.equal(dimnames(rsarray.rank)[[1]], dimnames(rsarray.rank)[[2]]) & 
      all.equal(dimnames(rsarray.line)[[1]], dimnames(rsarray.line)[[2]]) &
      all.equal(dimnames(rsarray.line)[[1]], dimnames(rsarray.rank)[[2]])
  )
  if(!are.rowcol.equal) stop("rsarray isn't rsarray!")
  
  ## get indices and values
  
  n.dim <- dim(rsarray.line)[1]
  is.lower.tri <- lower.tri(diag(n.dim))
  subjs <- dimnames(rsarray.line)$subj
  rois <- dimnames(rsarray.line)$roi
  n.mods <- length(subjs) * length(rois)
  
  ## unwrap into lower-triangle vector
  
  rsvectors <- vector("list", n.mods)
  names(rsvectors) <- combo_paste(subjs, rois)
  
  for (subj.i in seq_along(subjs)) {
    for (roi.j in seq_along(rois)) {
      # subj.i = 1; roi.j = 1
      
      rsm.line <- rsarray.line[, , subj.i, roi.j]  ## get slice
      rsm.rank <- rsarray.rank[, , subj.i, roi.j]
      
      r <- rsm.line[is.lower.tri]  ## get lower.triangle vector
      rank <- rsm.rank[is.lower.tri]
      
      name.ij <- paste0(subjs[subj.i], "_", rois[roi.j])  ## to match name
      rsvectors[[name.ij]] <- cbind(r, rank)
      
    }
  }
  
  ## check numbers
  
  n.rows.rsvectors <- map_dbl(rsvectors, nrow)
  if(sum(n.rows.rsvectors != 120) > 0) stop("missing row somewhere in rsvectors!")
  
  ## scale
  
  rsvectors.std <- lapply(rsvectors, scale)
  
  
  ## fit glms ----
  
  stats.subjs <- setNames(vector("list", length(m)), names(m))
  
  for (m.i in seq_along(m)) {
    # m.i = 3
    
    X <- m[[m.i]]
    X.std <- m.std[[m.i]]
    
    ## get unstandardized coefficients
    
    fits <- rsvectors %>% map(~ .lm.fit( x = X, y = .))  ## fit models (two-column y)
    coefs <- as.data.frame(do.call(rbind, lapply(fits, coef))) ## to single data.frame
    
    names(coefs) <- c("r", "rank")  ## same order as in models
    coefs$param <- rep(colnames(X), n.mods)  ## same orders as in models
    coefs$id <- rep(names(rsvectors), each = ncol(X))
    coefs <- coefs[coefs$param != "b0", ]  ## ditch intercept
    ## put rank and linear (r) coefs in single long-form column:
    coefs <- melt(as.data.table(coefs), id.vars = c("id", "param"), variable = "y", value.name = "coef")
    
    ## get betas
    ## NB: refitting model (with standardized coefficients), as opposed to scaling estimated coefficients,
    ## allows interactions to be specified (not applicable, however)
    
    fits.std <- rsvectors.std %>% map(~ .lm.fit( x = X.std, y = .))
    betas <- as.data.frame(do.call(rbind, lapply(fits.std, coef)))
    
    names(betas) <- c("r", "rank")
    betas$param <- rep(colnames(X.std), n.mods)
    betas$id <- rep(names(rsvectors.std), each = ncol(X.std))
    betas <- melt(as.data.table(betas), id.vars = c("id", "param"), variable = "y", value.name = "beta")
    
    ## get partial cors
    
    partr.r <- do.call(rbind, lapply(rsvectors, get.partr, colname = "r", method = "pearson"))
    partr.r <- tibble::rownames_to_column(as.data.frame(partr.r), "id")
    partr.r <- reshape2::melt(partr.r, variable = "param", value.name = "partr", id.vars = "id")
    
    partr.rank <- do.call(rbind, lapply(rsvectors, get.partr, "rank", method = "spearman"))
    partr.rank <- tibble::rownames_to_column(as.data.frame(partr.rank), "id")
    partr.rank <- reshape2::melt(partr.rank, variable = "param", value.name = "partr", id.vars = "id")
    
    partr.r$y <- "r"
    partr.rank$y <- "rank"
    partr <- rbind(partr.r, partr.rank)
    
    ## bind and save
    
    d <- full_join(coefs, betas, by = c("y", "id", "param"))
    stats.subjs[[m.i]] <- full_join(d, partr, by = c("y", "id", "param"))
    
  }
  
  
  ## format ----
  
  ## collate into data.frame
  
  stats.subjs <- bind_rows(stats.subjs, .id = "model")
  
  ## create subj, roi, and hemi cols from id col
  
  stats.subjs <- bind_cols(
    stats.subjs,
    reshape2::colsplit(stats.subjs$id, pattern = "_", names = c("subj", "roi"))
  )
  
  ## add is.analysis.group col
  
  stroop <- fread(here("data", "behavior-and-events_group201902.csv"))
  sample.analysis <- unique(filter(stroop, is.analysis.group)$subj)
  stats.subjs$is.analysis.group <- stats.subjs$subj %in% sample.analysis
  
  ## rearrange cols (and drop id col)
  
  stats.subjs %<>% select(subj, is.analysis.group, roi, model, y, param, coef, beta, partr)
  
  
  ## write ----
  ## NB. break into smaller files to keep within github limit of 100 MB / file
  ## break by model structure
  
  for (ii in seq_along(structures)) {
    
    model.i <- names(structures)[ii]
    stats.i <- filter(stats.subjs, model == model.i)
    
    fwrite(
      stats.i,
      here(
        "out", "rsa", "stats", 
        paste0("subjs_", glmname, "_", set.i, "_pearson_residual_glm-", model.i, ".csv")
      )
    )
    
  }
  
}
