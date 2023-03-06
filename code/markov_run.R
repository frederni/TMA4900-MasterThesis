#' ---
#' title: "Main notebook"
#' output: pdf_document
#' date: '2023-02-01'
#' ---
#' 
## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------

#' 
#' # Acknowledgements
#' 
#' The data and a large portion of data preprocessing is provided by Jane Reid. The re-implementation into INLA is also largely based on the work from Stefanie Muff.
#' 
#' # Data loading
#' 
## ----data loading-------------------------------------------------------------------------------------------------------------------------------------
# library(MCMCglmm)
library(MASS)
library(nadiv)
library(bdsmatrix)
library(INLA)
# library(SMisc)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("GeneticsPed", quietly = TRUE))
    BiocManager::install("GeneticsPed")
# if (!require("MCMCglmm", quietly = TRUE))
#     install.packages("../MCMCglmm-rbv-patch.tar.gz")

library("GeneticsPed") 
library(QGglmm)
# library("MCMCglmm")


qg.data.gg.inds <- read.table("../data/qg.data.gg.inds.steffi.txt", header=T)
d.ped <- ped.prune.inds <- read.table("../data/ped.prune.inds.steffi.txt", header=T)
d.Q <-  read.table("../data/Q.data.steffi.txt", header=T)

qg.data.gg.inds$natalyr.id <- qg.data.gg.inds$natalyr.no


#' 
#' Scaling continuous variances can be more stable: (I think?)
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------
qg.data.gg.inds$f.coef.sc <- scale(qg.data.gg.inds$f.coef,scale=FALSE)
qg.data.gg.inds$g1.sc <- scale(qg.data.gg.inds$g1,scale=FALSE)
qg.data.gg.inds$natalyr.no.sc <- scale(qg.data.gg.inds$natalyr.no,scale=FALSE)
qg.data.gg.inds$brood.date.sc <- scale(qg.data.gg.inds$brood.date,scale=FALSE)

#' 
#' The sex covariate is either `1` or `2`, so we binarize this covariate:
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------
qg.data.gg.inds$sex <- qg.data.gg.inds$sex.use.x1 - 1 

#' 
#' ## Deriving *A*
#' 
#' In order to derive the A matrix, we need to work a bit
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------
d.ped <- nadiv::prepPed(d.ped)

#' 
#' In particular, for INLA we need ids that run from 1 to the number of individuals
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------
d.ped$id <- 1:(nrow(d.ped))

#' 
#' Need a map file to keep track of the ids
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------
d.map <- d.ped[,c("ninecode","id")]
d.map$g1 <- d.Q[match(d.map$ninecode,d.Q$ninecode),"g1"]
d.map$foc0 <- d.Q[match(d.map$ninecode,d.Q$ninecode),"foc0"]

#' 
#' Give mother and father the id
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------
d.ped$mother.id <- d.map[match(d.ped$gendam, d.map$ninecode),"id"]
d.ped$father.id <- d.map[match(d.ped$gensire, d.map$ninecode),"id"]

#' 
#' Make the inverse A matrix using the `nadiv` package:
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------
Cmatrix <- nadiv::makeAinv(d.ped[,c("id","mother.id","father.id")])$Ainv

#' 
#' Store the id twice: Once for the breeding value, and once for the independent residuals u with variance 1 (the latter are not going to be included in the end, but we checked what happened when they were there)
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------
qg.data.gg.inds$id <- d.map[match(qg.data.gg.inds$ninecode, d.map$ninecode), "id"]
qg.data.gg.inds$u <- 1:nrow(qg.data.gg.inds)

#' 
#' # INLA
#' 
#' The general INLA formula is provided below, where `f()` encode the random effect:
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------

formula.inla.scaled = surv.ind.to.ad ~ f.coef.sc + g1.sc + natalyr.no.sc + brood.date.sc + sex +
  f(nestrec, model="iid",hyper=list(
    prec=list(initial=log(1/0.05), prior="pc.prec",param=c(1,0.05)) # PC priors
  )) +
  f(natalyr.id, model="iid",hyper=list(
    prec=list(initial=log(1/0.25), prior="pc.prec",param=c(1,0.05)) # PC priors
  )) +
  f(id,model="generic0", # Here we need to specify the covariance matrix 
    Cmatrix=Cmatrix,     #    via the inverse (Cmatrix)
    constr = F, # Doesn't really matter
    hyper=list(
      prec=list(initial=log(1/10), prior="pc.prec",param=c(1,0.05)) # PC priors
     ))  # + # this last part is only needed if the u~N(0,1) residuals are included
# f(u, model="iid",
#   constr=TRUE,
#   hyper=list(
#     prec=list(initial=log(1), fixed=TRUE) # Fixed variance to 1
#   )) # This last component adds the independent residuals with variance 1;

#' 
#' Now we call INLA models, this is the slowest part:
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------

fit.inla.probit = inla(formula=formula.inla.scaled, family="binomial",
                       data=qg.data.gg.inds,
                       control.compute=list(dic=T, return.marginals.predictor=T),
                       control.family = list(link = "probit"),
)

fit.inla.logit = inla(formula=formula.inla.scaled, family="binomial",
                      data=qg.data.gg.inds,
                      control.compute=list(dic=T, return.marginals.predictor=T),
                      control.family = list(link = "logit"),
                      control.predictor = list(compute=T)
)

fit.inla.gaussian = inla(formula=formula.inla.scaled, family="gaussian",
                             data=qg.data.gg.inds,
                             control.compute=list(dic=T, cpo=T) 
)


inla.posterior.marginal.latent.mode <- function(fit){
  modes = c()
  iter = 1
  for(predictor in names(fit$marginals.linear.predictor)){
    xy = get(predictor, fit$marginals.linear.predictor)
    modes[iter] = xy[, "x"][which.max(xy[, "y"])]
    iter = iter + 1
  }
  modes
}

marginal.latent.samples <- function(fit, nsamples){
  #' TODO : Work in progress!!
  #' TODO : improve runtime, seems like inla.rmarginal takes the longest time.
  #' What if we sample.marginal `nsamples` samples for each predictor
  #' Output is list(c(...), c(...), ..., c(...)) w/ nsamples list elems,
  #' and 2000ish elements in each c() (the different predictor samples)
  out_transpose = matrix(nrow=nsamples, data=rep(0,nsamples))
  cat("DEBUG: Entering sampling for loop\n")
  iter = 0
  for(predictor in names(fit$marginals.linear.predictor)){
    iter = iter + 1
    xy = get(predictor, fit$marginals.linear.predictor)
    predictor_samples = inla.rmarginal(nsamples, xy)
    out_transpose = cbind(out_transpose, predictor_samples)
  }
  # Need to slice off the first 0-column:
  out_transpose = out_transpose[, 2:ncol(out_transpose)]
  cat("DEBUG: Entering transposing loop\n")
  # We want list where each list element is one sample, i.e. the transposed
  # Start naive method, might be optimized:
  out = list()
  for(i in 1:nsamples){
    out[[i]] = out_transpose[i, ]
  }
  return(out)
}



get.h2 <- function(inla.fit, n, use.scale=F, model=NA){
  #' Get n samples of heritability from INLA fit
  #' h^2 computed as reciprocal of precision for id, over the sum of 
  #' (the reciprocal of) all random effects / hyperparameters 
  samples <- inla.hyperpar.sample(n=n,inla.fit)
  denominator = 0
  for(cname in colnames(samples)){ denominator = denominator + 1/samples[, cname]}
  
  if(use.scale){
    # We need model specification to use scale
    stopifnot(model %in% c("binom1.probit", "binom1.logit"))
    scale.param = ifelse(model == "binom1.probit", 1, pi^2/3)
    denominator = denominator + scale.param
  }
  
  h2.inla <- (1/samples[,"Precision for id"]) / denominator
  return(h2.inla)
}

get.h2.from.qgparams <- function(inla.fit, modelname, n, n.obs=nrow(d.ped)){
  #' Computes a posterior of data-scale heritability (h2) using QGParams
  #'  
  #' Params:
  #' inla.fit    the fitted INLA object
  #' modelname  a string specifying model
  #'    ("Gaussian", "binomN1.probit", "binom1.logit")
  #' n.obs      keyword argument if binomial model is has N != 1 trials
  #'    NB: Should explicitly be set to NULL if not relevant (Gaussian model)
  #'    Currently not in use as we're always dealing with binary (1 trial)
  stopifnot(modelname %in% c("Gaussian", "binom1.probit", "binom1.logit"))
  samples.posterior <- inla.hyperpar.sample(n=n,inla.fit)
  vp.samples = 0
  for(cname in colnames(samples.posterior)){
    vp.samples = vp.samples + 1/samples.posterior[, cname]
    }
  mu = inla.fit$summary.fixed$mean[1]
  va.samples = 1/samples.posterior[,"Precision for id"]
  kwargs = list(verbose=F)
  # out = mapply(QGparams, mu, va.samples, vp.samples, modelname, MoreArgs = kwargs)
  h2.getter = function(...){
    get("h2.obs", suppressWarnings(QGparams(...)))
    }
  out = mapply(h2.getter, mu, va.samples, vp.samples, modelname, MoreArgs=kwargs)
  return(out)
}


new.h2.transf <- function(fit, modelname, nsamples, marginal.mode=NA){
  # This is very much WIP but works I think now
  posterior.samples = inla.hyperpar.sample(nsamples, fit)
  
  vp.samples = 0
  for(cname in colnames(posterior.samples)){
    vp.samples = vp.samples + 1/posterior.samples[, cname]
    }
  
  df <- data.frame(va = as.vector(1/posterior.samples[, "Precision for id"]),
                   vp = as.vector(vp.samples)
                   )
  if(is.na(marginal.mode[1])){
    df$predict = marginal.latent.samples(fit, nsamples)
    posterior = do.call("rbind", apply(df, 1, function(row){
      QGparams(predict=row[["predict"]], var.a=row[["va"]], var.p=row[["vp"]],
               model=modelname, verbose=F)
    }))
  }
  else{
    posterior = do.call("rbind", apply(df, 1, function(row){
      QGparams(predict=marginal.mode, var.a=row[["va"]], var.p=row[["vp"]],
               model=modelname, verbose=F)
    }))
  }
  posterior
}

compute.observation.h2.SSH <- function(fit, nsamples, modelname){
  cat("Called function at:", format(Sys.time(), "%a %b %d %X %Y"), "\n")
  t0 = Sys.time()
  latent.heritability = get.h2(fit, nsamples)
  scaled.heritability = get.h2(fit, nsamples, use.scale=T, model=modelname)
  
  ti = Sys.time()
  cat("Entering `new.h2.transf` (Bayesian sampling)\n")
  heritability_averaged = new.h2.transf(fit, modelname, nsamples)
  cat("Done after", difftime(Sys.time(),ti, units="mins"), "minutes. Entering same function without the Bayesian stuff.\n")
  
  ti = Sys.time()
  marginal.mode = inla.posterior.marginal.latent.mode(fit)
  heritability_frequentist_avged = new.h2.transf(fit, modelname, nsamples, marginal.mode=marginal.mode)
  cat("Done after", difftime(Sys.time(),ti, units="mins"), "minutes. Now without averaging:")
  
  ti = Sys.time()
  heritability_notavged = get.h2.from.qgparams(fit, modelname, nsamples)
  cat("Done after", difftime(Sys.time(),ti, units="mins"), "minutes.\n")
  
  cat("Saving to disk...\n")
  save(latent.heritability, scaled.heritability, heritability_averaged,
       heritability_frequentist_avged, heritability_notavged,
       file=paste0("heritabilities_SSH_", modelname, ".Rdata"))
  return(0)
}

compute.observation.h2.SSH(fit.inla.logit,  10000, "binom1.logit")
compute.observation.h2.SSH(fit.inla.probit, 10000, "binom1.probit")



