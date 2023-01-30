#'---
#'title: Analysis of juvenile survival from independence to adulthood
#'author: "S. Muff"
#'date: "`r format(Sys.time(), '%d %B, %Y')`"
#'output:
#'  html_document:
#'    toc: yes
#'---

#'  Analysis of juvenile survival from independence to adulthood
#'  Original analysis carried out by Jane Reid
#'  Here we implement the same model using INLA


#+ include = FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)

library(MCMCglmm)
library(MASS)
library(nadiv)
library(bdsmatrix)
library(INLA)

library(SMisc)

#' # Preparing the data

qg.data.gg.inds <- read.table("../data/qg.data.gg.inds.steffi.txt", header=T)
d.ped <- ped.prune.inds <- read.table("../data/ped.prune.inds.steffi.txt", header=T)
d.Q <-  read.table("../data/Q.data.steffi.txt", header=T)



qg.data.gg.inds$natalyr.id <- qg.data.gg.inds$natalyr.no

#' I like to scale the continuous variances, especially in a Bayesian approach
qg.data.gg.inds$f.coef.sc <- scale(qg.data.gg.inds$f.coef,scale=FALSE)
qg.data.gg.inds$g1.sc <- scale(qg.data.gg.inds$g1,scale=FALSE)
qg.data.gg.inds$natalyr.no.sc <- scale(qg.data.gg.inds$natalyr.no,scale=FALSE)
qg.data.gg.inds$brood.date.sc <- scale(qg.data.gg.inds$brood.date,scale=FALSE)

#' Make sex binary (0-1 rather than 1-2):
qg.data.gg.inds$sex <- qg.data.gg.inds$sex.use.x1 - 1 

#' Look at the data once more
#+ eval = FALSE
str(qg.data.gg.inds)


#' # MCMCglmm
#' 


run.MCMC_chain <- function(){
  
  #' Some data preparation that seems obviously needed to run the MCMCglmm code:
  qg.data.gg.inds$animal <- qg.data.gg.inds$ninecode
  qg.data.gg.inds$sex.factor <- as.factor(qg.data.gg.inds$sex.use.x1)
  qg.data.gg.inds$natalyr.factor <- as.factor(qg.data.gg.inds$natalyr)
  qg.data.gg.inds$nestrec.factor <- as.factor(qg.data.gg.inds$nestrec)
  
  #' It's possible to play a bit with different variances for the residuals; My observation is that the chains have a mixing problem when the choice is too small. V=1 seems to mix well, V=10 even better... 
  VR <- 1
  prior = list(R = list(V = VR, fix=T), 
               G = list(
                 G1 = list(V = 1, nu = 0, alpha.mu = 0, alpha.V = 1000),
                 G2 = list(V = 1, nu = 0, alpha.mu = 0, alpha.V = 1000),
                 G3 = list(V = 1, nu = 0, alpha.mu = 0, alpha.V = 1000))
               )
  model.surv.ind.ad.gg1 <- MCMCglmm(
    surv.ind.to.ad ~ f.coef.sc + g1.sc + natalyr.no.sc + brood.date.sc + sex,
    random = ~ animal + nestrec.factor + natalyr.factor, 
    # family = "ordinal", # probit regression, to compare to INLA
    family = "categorical", # logistic regression (what Jane does on paper)
    # family = "multinomial2",
    data = qg.data.gg.inds, prior = prior,
    pedigree = ped.prune.inds, 
    verbose = F,
    pr = TRUE,
    burnin = 1000, nitt = 11000, thin=5
    )
  
  summary(model.surv.ind.ad.gg1)
  
  plot(model.surv.ind.ad.gg1$VCV[,"animal"])
  plot(model.surv.ind.ad.gg1$VCV[,"nestrec.factor"])
  plot(model.surv.ind.ad.gg1$VCV[,"natalyr.factor"])
  
  h2.mcmc <- model.surv.ind.ad.gg1$VCV[,"animal"]/(
    model.surv.ind.ad.gg1$VCV[,"animal"] +
      model.surv.ind.ad.gg1$VCV[,"nestrec.factor"] +
      model.surv.ind.ad.gg1$VCV[,"natalyr.factor"] + 
      VR
    )
  
  truehist(h2.mcmc)
  
  mean(h2.mcmc)  
} 

#' # INLA
#' ## Homogeneous variances

#' In order to derive the A matrix, we need to work a bit
d.ped <- nadiv::prepPed(d.ped)

#' In particular, for INLA we need ids that run from 1 to the number of individuals
d.ped$id <- 1:(nrow(d.ped))

#' Need a map file to keep track of the ids
d.map <- d.ped[,c("ninecode","id")]
d.map$g1 <- d.Q[match(d.map$ninecode,d.Q$ninecode),"g1"]
d.map$foc0 <- d.Q[match(d.map$ninecode,d.Q$ninecode),"foc0"]

#' give mother and father the id
d.ped$mother.id <- d.map[match(d.ped$gendam, d.map$ninecode),"id"]
d.ped$father.id <- d.map[match(d.ped$gensire, d.map$ninecode),"id"]


#' Make the inverse A matrix using the nadiv package:
Cmatrix <- nadiv::makeAinv(d.ped[,c("id","mother.id","father.id")])$Ainv

#' Store the id twice: Once for the breeding value, and once for the independent residuals u with variance 1 (the latter are not going to be included in the end, but we checked what happened when they were there)
qg.data.gg.inds$id <- d.map[match(qg.data.gg.inds$ninecode, d.map$ninecode), "id"]
qg.data.gg.inds$u <- 1:nrow(qg.data.gg.inds)

#' inla formula, where f() encode the random effect:
formula.surv.ind.to.ad = surv.ind.to.ad ~ f.coef.sc + g1.sc + natalyr.no.sc + brood.date.sc + sex +
  f(nestrec,model="iid",hyper=list(
    prec=list(initial=log(1/0.05), prior="pc.prec",param=c(1,0.05)) # PC priors
  )) +
  f(natalyr.id, model="iid",hyper=list(
    prec=list(initial=log(1/0.25), prior="pc.prec",param=c(1,0.05)) # PC priors
  )) +
  f(id,model="generic0", # Here we need to specify the covariance matrix via the inverse (Cmatrix)
    Cmatrix=Cmatrix,
    # constr = TRUE,
    hyper=list(
      prec=list(initial=log(1/10), prior="pc.prec",param=c(1,0.05)) # PC priors
    ))  #+ # this last part is only needed if the u~N(0,1) residuals are included
# f(u, model="iid",
#   constr=TRUE,
#   hyper=list(
#     prec=list(initial=log(1), fixed=TRUE) # Fixed variance to 1
#   )) # This last component adds the independent residuals with variance 1;

#' Now we call INLA, this is the slowest part:
r.inla.surv.ind.to.ad = inla(formula=formula.surv.ind.to.ad, family="binomial",
                             data=qg.data.gg.inds,
                             control.compute=list(dic=T), # config = TRUE
                             control.family = list(link = "logit")
)

#' Check if there is a problem (ok of =0) and dic
r.inla.surv.ind.to.ad$mode$mode.status
r.inla.surv.ind.to.ad$dic$dic

summary(r.inla.surv.ind.to.ad)

r.inla.surv.ind.to.ad$summary.hyperpar
r.inla.surv.ind.to.ad$summary.fixed

inla_mmarginal(r.inla.surv.ind.to.ad)
inla_emarginal(r.inla.surv.ind.to.ad)

#' Plotting the posterior marginals for the variances; This is a bit cumbersome, because INLA works with precisions, so we need a transformation. The code below does it for us:
par(mfrow=c(1,3))
plot(inla.tmarginal(function(x) 1/x,r.inla.surv.ind.to.ad$marginals.hyperpar$`Precision for nestrec`),type="l",main = "next")
plot(inla.tmarginal(function(x) 1/x,r.inla.surv.ind.to.ad$marginals.hyperpar$`Precision for natalyr.id`),type="l",main="natalyr")
plot(inla.tmarginal(function(x) 1/x,r.inla.surv.ind.to.ad$marginals.hyperpar$`Precision for id`),type="l",main="animal")

#' To obtain the posterior marginal of heritability, we need to resample from the posterior of the hyperparameters (this is very efficient):

nsamples <- 100000
sample.posterior <- inla.hyperpar.sample(n=nsamples,r.inla.surv.ind.to.ad)

h2.inla <- 1/sample.posterior[,"Precision for id"] / ((
  1/sample.posterior[,"Precision for id"]) + (1/sample.posterior[,"Precision for natalyr.id"]) +
    (1/sample.posterior[,"Precision for nestrec"]) + 0)

par(mfrow=c(1,1))
truehist(h2.inla)


#### Start QGglmm?
install.packages("QGglmm")
library(QGglmm)
QGparams(
  mu= r.inla.surv.ind.to.ad$summary.fixed$mean[1],
  var.a = mean(1/sample.posterior[,"Precision for id"]),
  var.p = mean((1/sample.posterior[,"Precision for id"]) +
    (1/sample.posterior[,"Precision for natalyr.id"]) +
    (1/sample.posterior[,"Precision for nestrec"])),
  model = "binom1.logit",
  n.obs = nrow(d.ped)
)
# ___ QUESTIONS ___
# * It's correct that we have a binomial case (binary trait) with N trials, each trial corresponding to an individual,
#       so that n.obs = 2722 (nrow(d.ped)) ?
# * Does it matter if I choose to look at the mean or mode of the re-sampled variances?
# * In Gaussian INLA, it makes more sense to me to just not specify any link functions, considering that
#     in gaussian regression, the response is eta so g(\eta)=\eta \iff g = identity (right?)
# Back-transformed INLA (binomial-probit) yields h^2 0.1781276 on data scale.

#### End QGglmm

#### Start gaussian INLA
r.inla.gaussian.surv.ind.to.ad = inla(formula=formula.surv.ind.to.ad, family="gaussian",
                             data=qg.data.gg.inds,
                             control.compute=list(dic=T), # config = TRUE
                             #control.family = list(link = "probit")
)

sample.gaussian.posterior <- inla.hyperpar.sample(n=nsamples,r.inla.gaussian.surv.ind.to.ad)
h2.gaussian.inla <- 1/sample.gaussian.posterior[,"Precision for id"] / (
  (1/sample.gaussian.posterior[,"Precision for id"]) +
    (1/sample.gaussian.posterior[,"Precision for natalyr.id"]) +
    (1/sample.gaussian.posterior[,"Precision for nestrec"]) #+1
  )
mean(h2.gaussian.inla) # Gives us 0.2025554 no?
truehist(h2.gaussian.inla)
quantile(h2.gaussian.inla, probs=c(0.05, 0.95))
#### End gaussian INLA


################### Not in use at the moment

#' ## Heterogeneous variances

#' make A (I was a bit lazy, because I have the heterogeneous VA code for the A matrix, but would be easy to make it for Ainv directly)
A <- nadiv::makeA(d.ped[,c("id","mother.id","father.id")])

#' Cholesky decomposition of the A matrix
tmp.full <- gchol(as.matrix(A)) 

D <- diag(tmp.full)

TT <- round(as(as.matrix(tmp.full), "dtCMatrix"),8)
TT[1:10,1:10]


#' Do the scaling on the full A matrix and only then reduce to AA0 (same for AA1 and AA2)
scaling0 <- ifelse(d.map$foc0>0,(d.map$foc0),1e-6)
scaling1 <- ifelse(d.map$g1>0,(d.map$g1),1e-6)

DD0 <- Diagonal(x=1-scaling0*(1-D)) 
TT0 <- TT %*%  Diagonal(x=scaling0)
A0 <- crossprod(t(TT0 %*% sqrt(DD0))) 
Cmatrix0 <- round(solve(A0),5)
Cmatrix0 <- as(Cmatrix0,"dgCMatrix")
Cmatrix0[Cmatrix0==0] <- 0

DD1 <- Diagonal(x=1-scaling1*(1-D)) 
TT1 <- TT %*%  Diagonal(x=scaling1)
A1 <- crossprod(t(TT1 %*% sqrt(DD1))) 
Cmatrix1 <- round(solve(A1),5)
Cmatrix1 <- as(Cmatrix1,"dgCMatrix")
Cmatrix1[Cmatrix1==0] <- 0

# Need to make a new column to use same id twice:
qg.data.gg.inds$id2 <- qg.data.gg.inds$id

qg.data.gg.inds$IndexA0 <- qg.data.gg.inds$id # ifelse(d.map[match(qg.data.gg.inds$ninecode,d.map$ninecode),"foc0"]>0, qg.data.gg.inds$id, NA)  
qg.data.gg.inds$IndexA1 <- qg.data.gg.inds$id # ifelse(d.map[match(qg.data.gg.inds$ninecode,d.map$ninecode),"g1"]>0, qg.data.gg.inds$id, NA)  

formula.2groups = surv.ind.to.ad ~ f.coef.sc + g1.sc + natalyr.no.sc + brood.date.sc + sex +
  f(nestrec,model="iid",hyper=list(
    prec=list(initial=log(1/0.05), prior="pc.prec",param=c(1,0.05))
  )) +
  f(natalyr.id, model="iid",hyper=list(
    prec=list(initial=log(1/0.25), prior="pc.prec",param=c(1,0.05))
  )) +
  f(IndexA0,model="generic0",
    Cmatrix=Cmatrix0,
    # constr = TRUE,
    hyper=list(
      prec=list(initial=log(10), prior="pc.prec",param=c(1,0.05))
    )) +
  f(IndexA1,model="generic0",
    Cmatrix=Cmatrix1,
    # constr = TRUE,
    hyper=list(
      prec=list(initial=log(10), prior="pc.prec",param=c(1,0.05))
    )) 

r.inla.2groups = inla(formula=formula.2groups, family="binomial",
                      data=qg.data.gg.inds,
                      control.compute=list(dic=T),
                      control.family = list(link = "probit")
                      # control.compute=list(config = TRUE)
)

r.inla.2groups$mode$mode.status
r.inla.2groups$dic$dic

summary(r.inla.2groups)

r.inla.2groups$summary.hyperpar
r.inla.2groups$summary.fixed

inla_mmarginal(r.inla.2groups)
inla_emarginal(r.inla.2groups)

par(mfrow=c(1,3))
plot(inla.tmarginal(function(x) 1/x,r.inla.2groups$marginals.hyperpar$`Precision for nestrec`),type="l",main="nest")
plot(inla.tmarginal(function(x) 1/x,r.inla.2groups$marginals.hyperpar$`Precision for natalyr.id`),type="l",main="natalyr")
plot(inla.tmarginal(function(x) 1/x,r.inla.2groups$marginals.hyperpar$`Precision for IndexA0`),ylim=c(0,10),xlim=c(0,1),lwd=2,type="l",main="animal")
lines(inla.tmarginal(function(x) 1/x,r.inla.2groups$marginals.hyperpar$`Precision for IndexA1`),lwd=2,col=2)
legend("topright",legend=c("Resident VA","Immigrant VA"),col=1:2,lty=1)



