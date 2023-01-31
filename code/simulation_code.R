###############################################
########### R code for simulation #############
###############################################

# load libraries
library("synbreed")
library("pedigree")
library("MCMCglmm")
library("mvtnorm")
library("GeneticsPed") 
library("pedigreemm")


fun.simulation <- function(nsim){
  
  h_err <- c()
  VA <- c()
  VE <- c()
  VPE <- c()
  h_noerr <- c()
  err.est <- c()
  VERR <- c()
  VA_noerr <- c()
  VE_noerr <- c()
  VPE_noerr <- c()
  sigmaE_mean <- c()
  sigma_err_mean <- c()
  h_true<- c()
  VA_true<- c()
  VE_true<- c()
  VPE_true<- c()
  
  for (jj in 1: nsim){
    
    # set a value for the Ne/Nc ratio
    # here is 0.5 but it can be set to 0.15 or 0.05
    NeNc <- 0.5
    # set other parameters 
    # (idgen = id per generations, nGen = number of generations)
    idgen <- 100
    nGen <- 9
    
    # generate the pedigree using the function generatePedigree from GeneticsPed
    # nFather and nMother are the number of fathers and mothers per generation
    # they are generated according to the selected Ne/Nc ratio
    ped0 <- generatePedigree(nId = idgen, nGeneration = nGen, 
                             nFather = idgen * NeNc, nMother = idgen * NeNc)  
    
    # N individuals
    Nid <- length(ped0$id)
    
    # Vector with animal id
    id <- 1 : Nid
    
    # fixed effects are mu and sex
    # mu is the baseline, here is chosen as 10 but it can be any value
    mu <- 10
    X <- rep(1, Nid)
    sex <- ped0$sex
    
    # set correct format for pedigree
    pedigree <- ped0[ , c(1,3,2)]
    names(pedigree) <- c("id", "dam", "sire")
    
    # generate breeding values (always constant)
    # using the function rbv from MCMCglmm
    # fix the additive variance value, here is 0.3 but it can be any value
    Z <- diag(Nid)
    sigmaA <- 0.4
    pedigree <- ped0[ ,1:3]
    u <- rbv(pedigree, sigmaA)
    
    # add repeated measurements to each individuals
    # number of observations
    sessions <- 10
    sigma_PE <- 0.25 # permanent environmental variance, costant across sessions
    # generate permanent environment effect (constant)
    id_eff <- rnorm(Nid, 0, sqrt(sigma_PE))
    
    # this quantities change across sessions
    session <- c()
    Y_rep <- c()
    Y_rep0 <- c()
    id_rep <- c()
    sex_rep <- c()
    for (s in 1:sessions){
      # number of repeats for esch individuals 
      repeats <- rbinom(Nid, 6, 0.5)
      repeats[repeats == 0] <- 1
      # E changes across session 
      sigmaE <- 0.3
      e <- rnorm(Nid, 0, sqrt(sigmaE))
      
      # generate repeated measurements for each individual
      # sex and id are always the same
      # Y is generated as many times per individual as its number of repeats
      # each time Y is generated from the covariates 
      # beta sex is 2 
      
      trait <- data.frame(sex_rep, id_rep, Y_rep, Y_rep0, session)
      for ( i in 1:Nid){
        Y1 <- c()
        Y0 <- c()
        session <- c(session, rep(s, repeats[i]))
        sex_rep <- c(sex_rep, rep(sex[i], repeats[i]))
        id_rep <- c(id_rep, rep(id[i], repeats[i]))
        
        sigma_err <- 0.5
        #sigma_err <- runif(1, 0.05, 0.1)
        for (k in 1: repeats[i]){
          Y1 <- mu+ 2*sex[i]+u[i]+e[i]+id_eff[i]+rnorm(1, 0, sqrt(sigma_err))
          Y0 <- mu+ 2*sex[i]+u[i]+e[i]+id_eff[i]
          Y_rep <- c(Y_rep, Y1)
          Y_rep0 <- c(Y_rep0, Y0)
        }
      }
      # store the dataset with all repeats
      trait <- rbind(trait, data.frame(sex_rep, id_rep, Y_rep, Y_rep0, session))
    }
    
    # additional id covariate to model permanent environment 
    trait$animal <- trait$id_rep
    trait$session <- as.factor(trait$session)
    trait$sex <- as.factor(trait$sex_rep)
    trait <- trait[ ,-1]
    names(trait) <- c('id', 'Y', 'Y0', 'session' ,'animal', 'sex')
    M <- c()
    for (m in 1:sessions) {
      M[m] <- length(trait[trait$session == m, 1])
    }
    
    trait <- trait[order(trait$session), ]
    trait$ids <- c(1:M[1], 1:M[2], 1:M[3], 1:M[4], 1:M[5], 
                  1:M[6], 1:M[7], 1:M[8], 1:M[9], 1:M[10])
    # MCMCglmm to compute heritability taking into account error!
    # prior assumes that the variance is equally divided into the four components
    p.var <- var(trait$Y,na.rm = TRUE)
    prior <- list(G = list(G1 = list(V = matrix(0.4), nu = 1),
                           G2 = list(V = matrix(0.25), nu = 1), 
                           G3 = list(V = matrix(0.5), nu = 1)),
                  R = list(V = matrix(0.3), nu =1))
    
    
    model_err <- MCMCglmm(Y~1+sex, random = ~animal+id+session:ids,
                          pedigree = pedigree, data = trait, prior = prior,
                          nitt = 11000, thin = 100, burnin = 1000, verbose = T)
    
    # heritability
    # it should not have the error component
    posterior.heritability_err <- model_err$VCV[,"animal"]/
      (model_err$VCV[,"animal"]+model_err$VCV[,"units"]+model_err$VCV[,"id"])
    h_err[jj] <- mean(posterior.heritability_err)
    
    # variance components
    VA[jj] <- mean(model_err$VCV[,"animal"])
    VE[jj] <- mean(model_err$VCV[,"units"])
    VPE[jj] <- mean(model_err$VCV[, "id"])
    VERR[jj] <- mean(model_err$VCV[ , "session:ids"])
    
    # true value on Y0 (without error) 
    # prior assumes that the variance is equally divided into the three components
    prior1 <- list(G = list(G1 = list(V = matrix(0.4), nu = 1),
                            G2 = list(V = matrix(0.25), nu = 1)),
                   R = list(V = matrix(0.8), nu = 1))
    
    model_true <- MCMCglmm(Y0~1+sex, random = ~animal+id,
                           pedigree = pedigree, data = trait, prior = prior1,
                           nitt = 11000, thin = 100, burnin = 1000, verbose = F)
    
    # heritability
    posterior.heritability_true <- model_true$VCV[,"animal"]/
      (model_true$VCV[,"animal"]+model_true$VCV[,"units"]+model_true$VCV[,"id"])
    h_true[jj] <- mean(posterior.heritability_true)
    
    # variance components
    VA_true[jj] <- mean(model_true$VCV[,"animal"])
    VE_true[jj] <- mean(model_true$VCV[,"units"])
    VPE_true[jj] <- mean(model_true$VCV[, "id"])
    
    
    # MCMCglmm to compute heritability without taking into account error!
    # prior assumes that the variance is equally divided into the three components
    
    prior2 <- list(G = list(G1 = list(V = matrix(0.4), nu = 1),
                            G2 = list(V = matrix(0.25), nu = 1)),
                   R = list(V = matrix(0.8), nu = 1))
    
    model_noerr <- MCMCglmm(Y~1+sex, random = ~animal+id,
                            pedigree = pedigree, data = trait, prior = prior2,
                            nitt = 11000, thin = 100, burnin = 1000, verbose = F)
    
    # heritability is given as VA/(VA+VPE+VE)
    # but it contains the error variance in the denominator! 
    posterior.heritability_noerr <- model_noerr$VCV[,"animal"]/
      (model_noerr$VCV[,"animal"]+model_noerr$VCV[,"units"]+model_noerr$VCV[,"id"])
    
    h_noerr[jj] <- mean(posterior.heritability_noerr)
    VA_noerr[jj] <- mean(model_noerr$VCV[,"animal"])
    VE_noerr[jj] <- mean(model_noerr$VCV[,"units"])
    VPE_noerr[jj] <- mean(model_noerr$VCV[, "id"])
    
    mm <- c()
    den <- c()
 
    for ( i in 1: Nid){
      for ( k in 1:sessions){
        for ( l in 1:length(trait[trait$session == k & trait$id == i, 1])){
          mm <- c(mm, (trait[trait$session == k & trait$id == i, ]$Y[l] - 
                         mean(trait[trait$session == k & trait$id == i, ]$Y, na.rm = TRUE))^2)
          
          
        }
        den<- c(den, length(trait[trait$session == k & trait$id == i, 1])-1)
        
      }
    }
    err.est[jj] <- sum(mm)/sum(den)
   
  }
  results <- data.frame(h_err, VA, VE, VPE, VERR,
                        h_noerr, VA_noerr, VE_noerr, VPE_noerr, 
                        h_true, VA_true, VE_true, VPE_true, err.est)
  
  # h if we just correct by subtracting the error from the denominator
  
  results$crude_h <- results$VA_noerr/
    (results$VA_noerr+results$VPE_noerr+results$VE_noerr-results$err.est)
  
  return(results)
}

# 100 simulated pedigrees
nsim <- 100
library(parallel)
cl = makeCluster(detectCores()-1)
clusterExport(cl,ls())
clusterExport(cl,c("orderPed", "generatePedigree", "calcInbreeding", "glmer", "MCMCglmm", "rbv"))
system.time(results <- parLapply(cl,nsim,fun.simulation))
stopCluster(cl)


write.table(results, 'results_simul_err0.5.txt')
results <- read.table('results_simul_1205.txt')


h <- c(results$h_err, results$h_noerr, results$h_true, results$crude_h)
VA <- c(results$VA, results$VA_noerr, results$VA_true)
VE <- c(results$VE, results$VE_noerr, results$VE_true)
VPE <- c(results$VPE, results$VPE_noerr, results$VPE_true)
#VERR <- c(results$h_err, results$h_noerr, results$h_true)
model <- c(rep('corrected',100), rep('biased',100), rep('true', 100))
model0 <- c(rep('corrected',100), rep('biased',100), rep('true', 100), rep('crude', 100))
f <- data.frame(VA, VE, VPE, model)
f0 <-  data.frame(h, model0)
postscript(file="results_simul.eps", width = 400, height = 600)
par(las = 1, mar = c(5,9,5,5))
boxplot(h~ model0, data = f0, ylab = "", 
        cex = 1.3, cex.axis = 2)
title(ylab = expression(hat(h^2)), line = 5, cex.lab = 2)
dev.off()

postscript(file="results_simulVA.eps", width = 400, height = 600)
par(las = 1, mar = c(5,9,5,5))
boxplot(VA~ model, data = f, ylab = "",
        cex = 1.3, cex.axis = 2.5)
title(ylab = expression(hat(sigma[A]^2)), line = 5, cex.lab = 2.5)
dev.off()

postscript(file="results_simulVE.eps", width = 400, height = 600)
par(las = 1, mar = c(5,9,5,5))
boxplot(VE~ model, data = f, ylab = "",
        cex = 1.3, cex.axis = 2.5)
title(ylab = expression(hat(sigma[E]^2)), line = 5, cex.lab = 2.5)
dev.off()


postscript(file="results_simulVPE.eps", width = 400, height = 600)
par(las = 1, mar = c(5,9,5,5))
boxplot(VPE~ model, data = f, ylab = "",
        cex = 1.3, cex.axis = 2.5)
title(ylab = expression(hat(sigma[PE]^2)), line = 5, cex.lab = 2.5)
dev.off()

VERR <- c(results$VERR, rep(0,100))
model1 <- c(rep('corrected', 100), rep('biased', 100))
f1 <- data.frame(VERR, model1)
postscript(file="results_simulVERR.eps", width = 400, height = 600)
par(las = 1, mar = c(5,9,5,5))
boxplot(VERR~ model1, data = f1, ylab = "",
        cex = 1.3, cex.axis = 2.5)
title(ylab = expression(hat(sigma[e]^2)), line = 5, cex.lab = 2.5)
dev.off()

