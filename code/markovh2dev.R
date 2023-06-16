# R script to re-run plot h2 dev and get errorbars:)

args = commandArgs(trailingOnly = TRUE)
if(length(args)==0) stop("No 'ntimes' parameter provided. Exiting.")

.libPaths("../rlib/")
req.packages <- c("BiocManager", "GeneticsPed", "ggplot2",
                  "latex2exp", "nadiv", "logger")
for (pack in req.packages){
  if(!require(pack, quietly = TRUE)) install.packages(pack,
                                                      lib="../rlib/")
}
if (!require("MCMCglmm", quietly = TRUE)) {
  install.packages("devtools", lib="../rlib/")
  devtools::install_local("MCMCglmm-rbv-patch.tar.gz", lib="../rlib/")
}
if(!require("INLA", quietly = TRUE)){
  install.packages(
    "INLA",
    repos=c(getOption("repos"),
            INLA="https://inla.r-inla-download.org/R/stable"),
    dep=TRUE)
}
library(ggplot2)
library(INLA)
library(GeneticsPed)
library(MCMCglmm)
library(latex2exp)
library(nadiv)
library(logger)
#################

get.h2 <- function(inla.fit, n, use.scale = FALSE, model = NA,
                   include.fixed = FALSE) {
  #' Get heritability
  #'
  #' Get n samples of heritability (h^2) from INLA object
  #' @param inla.fit fitted model
  #' @param n number of samples
  #' @param use.scale flag for adding link variance to denominator
  #' @param model string representation of model type
  #' @param include.fixed flag for including fixed effects variance
  #' @return n-sized vector of heritability samples
  samples <- inla.hyperpar.sample(n = n, inla.fit)
  denominator <- 0
  for (cname in colnames(samples)) {
    denominator <- denominator + 1 / samples[, cname]
  }
  if (include.fixed) {
    # Currently it just takes the SD instead of sampling,
    # didn't find an easy way to sample fixed effect estimates..
    denominator <- denominator + sum(inla.fit$summary.fixed[, "sd"]^2)
  }
  
  if (use.scale) {
    scales.dictionary <- list(
      binom1.probit = 1, binom1.logit = pi^2 / 3,
      round = 0.25
    )
    scale.param <- get(model, scales.dictionary)
    denominator <- denominator + scale.param
  }
  
  h2.inla <- (1 / samples[, "Precision for id"]) / denominator
  return(h2.inla)
}

threshold.scaling.param <- function(p) {
  # h^2_l = threshold.scaling.param * h^2_obs
  p * (1 - p) / (dnorm(qnorm(p)))^2
}

simulated.heritability <- function(NeNc = 0.5, idgen = 100, nGen = 9,
                                   sigmaA = 0.8, linear.predictor = NA,
                                   simulated.formula = NA,
                                   dichotomize = "round",
                                   pc.prior = NA, probit.model = FALSE,
                                   simulated.formula.probit = NA,
                                   DIC = FALSE) {
  #' Simulate and fit animal model
  #'
  #' Generate pedigree, fit Gaussian (INLA) model and provide
  #' heritability estimate.
  #' @param NeNc Effective/Census population mean, used to determine
  #' number of fathers and mothers per generation,
  #' @param idgen Number of individuals per generation in pedigree
  #' @param nGen Number of generations in pedigree
  #' @param sigmaA:  Additive genetic variance
  #' @param linear.predictor Callable function of two parameters, the
  #' first 'u' (breeding values), the second for the data
  #' @param simulated.formula Formula expression using the response
  #' name `simulated.response`, param `id` and `Cmatrix`, both of
  #' which are defined locally in this method.
  #' @param dichotomize Dichotomization method (round, binomial, etc)
  #' @param pc.prior (optional) parameters for PC prior
  #' @param probit.model (optional) flag to fit binomial probit model
  #' in addition to the Gaussian model.
  #' @param simulated.formula.probit (opt) Formula for probit model
  #' @param DIC (optional) Flag for computing DIC for the models
  #'
  #' @return A list with the following items:
  #' heritability: Posterior latent heritability samples
  #' summary: List of mean, standard deviation and quantiles of h^2
  #' p: Portion of `TRUE` observations in simulated response
  #' simulated.response:  The observed values in simulation dataset
  #' fit: Gaussian fitted model
  #' fit.probit: Probit fitted model, if `probit.model=F`, is `NULL`.
  ped0 <- generatePedigree(
    nId = idgen, nGeneration = nGen,
    nFather = idgen * NeNc, nMother = idgen * NeNc
  )
  # Set correct format for pedigree
  pedigree <- ped0[, c(1, 3, 2, 5)]
  names(pedigree) <- c("id", "dam", "sire", "sex")
  
  # Generate random breeding values
  # The following will CRASH if you don't
  # use the patched MCMCglmm package!
  u <- rbv(pedigree[, c(1, 2, 3)], sigmaA)
  
  simulated.d.ped <- nadiv::prepPed(pedigree, gender = "sex")
  # Binarize from (1,2) to (0,1)
  simulated.d.ped$sex <- simulated.d.ped$sex - 1 
  simulated.Cmatrix <- nadiv::makeAinv(pedigree[, c(1, 2, 3)])$Ainv
  # Make index to allow iid noise random effect
  simulated.d.ped$ind <- seq_len(nrow(simulated.d.ped))
  
  # Generating "true" y_i
  
  if (dichotomize == "binom1.logit") {
    simulated.response <- rbinom(length(u),
                                 size = 1,
                                 prob = pnorm(linear.predictor(u, simulated.d.ped))
    )
  } else if (dichotomize == "round") {
    # This assumes mean of \eta_i is 0
    simulated.response <- ifelse(
      linear.predictor(u, simulated.d.ped) <= 0, 0, 1
    )
  } else if (dichotomize == "round_balanced") {
    # Get balanced residuals for unbalanced linear predictor
    cutoff <- mean(linear.predictor(u, simulated.d.ped))
    simulated.response <- ifelse(
      linear.predictor(u, simulated.d.ped) <= cutoff, 0, 1
    )
  } else if (is.numeric(dichotomize)) {
    stopifnot(dichotomize >= 0 & dichotomize <= 1)
    eta_values <- linear.predictor(u, simulated.d.ped)
    cutoff <- quantile(eta_values, 1 - dichotomize)
    # e.g. dichotomize=0.1, p should be about 0.1
    simulated.response <- ifelse(eta_values <= cutoff, 0, 1)
  } else {
    stop(paste0(
      "Unknown dichotomization method '", dichotomize, "'. ",
      "Consider using 'binom1.logit' or 'round'."
    ))
  }
  p <- mean(simulated.response) # portion of true responses
  
  
  # Model fitting LMM for binary trait
  # First reload formula environment to access local variables
  environment(simulated.formula) <- environment()
  environment(simulated.formula.probit) <- environment()
  
  simulated.fit.inla <- inla(
    formula = simulated.formula, family = "gaussian",
    data = simulated.d.ped, control.compute = list(dic = DIC)
  )
  
  # Checks for error status in INLA fit,
  if (simulated.fit.inla$mode$mode.status != 0) {
    log_warn(("INLA status {simulated.fit.inla$mode$mode.status}"))
  }
  heritability <- get.h2(simulated.fit.inla, 10000) 
  
  # Also fit probit if specified
  if (probit.model) {
    fit.probit <- inla(
      formula = simulated.formula.probit, family = "binomial",
      data = simulated.d.ped,
      control.compute = list(return.marginals.predictor = TRUE,
                             dic = DIC)
    )
  } else {
    fit.probit <- NULL
  }
  list(
    heritability = heritability,
    summary = list(
      mean = mean(heritability),
      standard.deviation = sd(heritability),
      quantiles = quantile(heritability, probs = c(0.025, 0.5, 0.975))
    ),
    p = p,
    simulated.response = simulated.response,
    fit = simulated.fit.inla,
    fit.probit = fit.probit
  )
}

plot.h2.deviation <- function(
    dichotomize = "round",
    title = "Simulation heritability",
    SAVE.PLOT = TRUE, plot.fn = NA, sigma.scale = "log",
    lin.pred = NULL,
    dynamic.priors = FALSE, simulated.formula = NULL, Ve = NULL,
    fixedeffects = FALSE) {
  #' Plot h^2 estimate, alongside true value, for a series of V_A
  #'
  #' For each V_A, generate simulation and fit a Gaussian model.
  #' Then, plot the obtained h^2 for observation and liability scale,
  #' alongside the true value.
  #' @param dichotomize Dichotomization method, used for simulation
  #' @param title ggplot legend title used as title for all (sub)plots
  #' @param SAVE.PLOT flag for saving plot to disk
  #' @param plot.fn (optional) String to add to the end of the 
  #' filename, before file extension, when saving the plot.
  #' @param sigma.scale either "log" or "small", deciding what values,
  #' and the spacing between values of V_A to be iterated over.
  #' @param lin.pred linear predictor for simulation
  #' @param dynamic.priors Flag for changing model priors based on V_A
  #' @param simulated.formula Formula for simulatin
  #' @param Ve Residual variance, or fixed effects variance
  #' for that model
  #' @param fixedeffects Flag to determine if model has fixed effects
  #' @return List with item "p" for the ggplot object.
  
  if (sigma.scale == "log") {
    sigmaA.list <- c(1:10 %o% 10^(-3:3)) # Log scale[10^-3,  10^3]
  } else if (sigma.scale == "small") { # Linear scale [10^-3, 0.259]
    sigmaA.list <- seq(0.001, 0.26, by = 0.01)
  } else {
    stop("Unrecognized scale for sigmaA.")
  }
  pc.U.list <- c(rep(10, 10) %o% 10^(-3:3))
  error.flag <- FALSE
  estimates <- c()
  latent <- c()
  true.vals <- c()
  est.CI.u <- c()
  est.CI.l <- c()
  plist <- c()
  iter.num <- 0
  Ve0 <- Ve
  for (sigmaA in sigmaA.list) {
    cat(">")
    if (dynamic.priors) {
      iter.num <- iter.num + 1
      if (sigmaA < 1) {
        pc.prior <- c(1, 0.05)
      } else {
        pc.prior <- c(pc.U.list[iter.num], 0.05)
      }
    } else {
      pc.prior <- c(1, 0.05)
    }
    if (is.null(simulated.formula)) { # Defaults eta = a_i + e
      simulated.formula <- simulated.response ~  f(id,
                                                   model = "generic0",
                                                   Cmatrix = simulated.Cmatrix,
                                                   constr = FALSE,
                                                   hyper = list(
                                                     prec = list(initial = log(sigmaA), prior = "pc.prec",
                                                                 param = pc.prior)
                                                   )
      )
    }
    if (is.null(lin.pred)) { # The 'usual' linear predictor
      lin.pred <- function(u, .) u + rnorm(length(u))
    }
    result <- tryCatch(
      simulated.heritability(
        NeNc = 0.5, idgen = 100, nGen = 9, sigmaA = sigmaA,
        linear.predictor = lin.pred,
        simulated.formula = simulated.formula,
        dichotomize = dichotomize,
        pc.prior = pc.prior
        ),
      error=function(cond){
        message("INLA was not able to fit model! Traceback:\n", cond,
                "\nReturning NAs")
        error.flag <- TRUE
        return(
          list(heritability=NA,
               summary=list(mean=NA, standard.deviation=NA, quantiles=NA),
               p=NA, fit=NA, fit.probit=NA)
        )
      }
    )
    if(error.flag){
      latent <- c(latent, NA)
      true.vals <- c(true.vals, NA)
      estimates <- c(estimates, NA)
      est.CI.l <- c(est.CI.l, NA)
      est.CI.u <- c(est.CI.u, NA)
      plist <- c(plist, NA)
    }
    else {
      if (is.null(Ve)) {
        # Fallback residual variance
        Ve <- 1
      }
      if (fixedeffects) {
        # beta^2 * Var(x_fixedeffect):
        Ve <- Ve * var(result$simulated.response) 
        posterior <- get.h2(result$fit, 10000, include.fixed = TRUE)
      } else {
        posterior <- result$heritability
      }
      
      simulation.h2.true <- sigmaA / (sigmaA + Ve)
      
      latent <- c(latent, mean(posterior)) # Observatoin-level
      threshold.scaled.h2 <- threshold.scaling.param(result$p)*posterior
      
      true.vals <- c(true.vals, simulation.h2.true)
      estimates.CI <- quantile(threshold.scaled.h2,
                               probs = c(0.025, 0.975))
      estimates <- c(estimates, mean(threshold.scaled.h2))
      est.CI.l <- c(est.CI.l, estimates.CI[1])
      est.CI.u <- c(est.CI.u, estimates.CI[2])
      plist <- c(plist, result$p)
      # Reset Ve
      Ve <- Ve0
    }
  }
  res <- data.frame(
    estimates = estimates,
    true.vals = true.vals, latent = latent,
    est.CI.l = est.CI.l, est.CI.u = est.CI.u,
    sigmaA = sigmaA.list,
    plist = plist
  )
  # Plotting
  p <- ggplot(data = res, aes(x = sigmaA)) +
    geom_ribbon(aes(ymin = est.CI.l, ymax = est.CI.u), alpha = 0.1) +
    geom_line(aes(y = true.vals, color = "atrue"), size = 1.5) +
    geom_line(aes(y = estimates, color = "bliab"), size = 1.5) +
    geom_line(aes(y = latent, color = "obs"), size = 1.5) +
    xlab(TeX("$\\sigma_A^2$")) +
    ylab(TeX("$h^2$")) +
    scale_x_log10() +
    scale_color_manual(
      name = title,
      values = c(
        "atrue" = "darkred", "bliab" = "steelblue",
        "obs" = "chartreuse3"
      ),
      labels = c(
        expression("True " * h["liab"]^2),
        expression("Fitted " * h["liab"]^2),
        expression("Fitted " * h["obs"]^2)
      )
    ) +
    theme(text = element_text(size = 18), legend.text.align = 0)
  if (SAVE.PLOT) {
    ggsave(
      paste0(
        "../figures/simulation_deviance_",
        if (!is.na(plot.fn)) plot.fn, ".pdf"
      ), p + theme(legend.position = "none"),
      width = 20, height = 20,
      units = "cm"
    )
    # Save legend as separate plot
    p.legend <- cowplot::get_legend(p)
    pdf(paste0(
      "../figures/simulation_deviance",
      if (fixedeffects) "_fixedeffects", "_legend.pdf"
    ), width = 7.87402, height = 7.87402)
    grid.newpage()
    grid.draw(p.legend)
    dev.off()
  }
  return(list(p = p))
}


######################

multiple.h2.dev <- function(sigma.scale, ntimes,
                            title = "Simulation heritability", ...){
  if (sigma.scale == "log") {
    sigmaA.list <- c(1:10 %o% 10^(-3:3))
  } else if (sigma.scale == "small") {
    sigmaA.list <- seq(0.001, 0.26, by = 0.01)
  } else {
    stop("Unrecognized scale for sigmaA.")
  }
  n.sigma <- length(sigmaA.list)
  # Initialize containers: each col is one run
  all.h2obs <- matrix(ncol = ntimes, nrow = n.sigma)
  all.h2liab <- matrix(ncol = ntimes, nrow = n.sigma)
  all.truevals <- matrix(ncol = ntimes, nrow = n.sigma)
  for(i in 1:ntimes){
    cat(paste0("\n [Run ", i ,"/", ntimes, "]\n"))
    res <- plot.h2.deviation(SAVE.PLOT = FALSE, sigma.scale=sigma.scale, ...)
    all.h2obs[, i] <- res$p$data$latent
    all.h2liab[, i] <- res$p$data$estimate
    all.truevals[, i] <- res$p$data$true.vals
  }
  plot.data <- data.frame(
    model = rep(c("h2obs", "h2liab", "true"),
                each = n.sigma, times = 2),
    top = c(rowMeans(all.h2obs) + apply(all.h2obs, 1, sd),
            rowMeans(all.h2liab) + apply(all.h2liab, 1, sd),
            rowMeans(all.truevals) + apply(all.truevals, 1, sd)
    ),
    mid = c(rowMeans(all.h2obs), rowMeans(all.h2liab),
            rowMeans(all.truevals)
    ),
    btm = c(rowMeans(all.h2obs) - apply(all.h2obs, 1, sd),
            rowMeans(all.h2liab) - apply(all.h2liab, 1, sd),
            rowMeans(all.truevals) - apply(all.truevals, 1, sd)
    ),
    xax = rep(sigmaA.list, times=3)
  )
  return(plot.data)
}


# run
log_threshold(TRACE)
log_info("Running R file 'markovh2dev.R' with arguments:")
log_info(args)
nruns <- as.numeric(args[1])

log_info("Plotting h2 deviation for smaller values")
markov.result2<- multiple.h2.dev("small", nruns, dynamic.priors=TRUE)
log_info("Plotting h2 deviation, round, with dynamic priors")
markov.result1<- multiple.h2.dev("log", nruns, dynamic.priors=T)
log_info("Plotting h2 deviation, binom, with dynamic priors")
markov.result3<- multiple.h2.dev("log", nruns, dichotomize="binom1.logit", dynamic.priors=TRUE)

log_info("Done! Saving dataframes to Rdata")
save(markov.result1, markov.result2, markov.result3, file=paste0("markovh2dev_", nruns, "_runs.Rdata"))

## For fixed effects system;
linear_predictor_fixedeffects <- function(u, simulated.d.ped) {
  #' \Tilde{\eta} = a + N(0, varE) + betaSex x_{sex}
  varE <- 1
  betaSex <- 100
  out <- c()
  intercept <- 0
  residuals <- rnorm(length(u), mean = 0, sd = sqrt(varE))
  for (idx in seq_along(u)) {
    out <- c(
      out,
      intercept + betaSex * simulated.d.ped$sex[idx] + u[idx] +
        residuals[idx]
    )
  }
  out
}
simulated.formula.fixedeffects <- simulated.response ~ sex +
  f(id,
    model = "generic0",
    Cmatrix = simulated.Cmatrix,
    constr = FALSE,
    hyper = list(
      prec = list(
        initial = log(1 / 10), prior = "pc.prec",
        param = c(1, 0.05)
      ) # PC priors
    )
  )

log_info("Fixed effects for beta=100, p balanced")
res.fixed1 <- multiple.h2.dev(
  "log", ntimes, dichotomize = "round_balanced",
  lin.pred = linear_predictor_fixedeffects,
  simulated.formula = simulated.formula.fixedeffects,
  Ve = 100^2, fixedeffects = TRUE, dynamic.priors=TRUE)

# Redfining linear predictor to smaller beta:
linear_predictor_fixedeffects <- function(u, simulated.d.ped) {
  varE <- 1
  betaSex <- 10
  out <- c()
  intercept <- 0 #-4.5
  residuals <- rnorm(length(u), mean = 0, sd = sqrt(varE))
  for (idx in seq_along(u)) {
    out <- c(
      out,
      intercept + betaSex * simulated.d.ped$sex[idx] + u[idx] +
        residuals[idx]
    )
  }
  out
}

log_info("Fixed effects for beta=10, p balanced")
res.fixed2 <- multiple.h2.dev(
  "log", ntimes, dichotomize = "round_balanced",
  lin.pred = linear_predictor_fixedeffects,
  simulated.formula = simulated.formula.fixedeffects,
  Ve = 10^2, fixedeffects = TRUE, dynamic.priors=TRUE)

log_info("Fixed effects for beta=10, p unbalanced 0.1")
res.fixed3 <- multiple.h2.dev(
  "log", ntimes, dichotomize = 0.1,
  lin.pred = linear_predictor_fixedeffects,
  simulated.formula = simulated.formula.fixedeffects,
  Ve = 10^2, fixedeffects = TRUE, dynamic.priors=TRUE)

log_info("Done! Saving results from fixed effects")
save(res.fixed1, res.fixed2, res.fixed3, file=paste0("markovfixed_", nruns, "_runs.Rdata"))

log_info("Done! Exiting")

## Plotter to use locally
plotter <- function(df){
  ggplot(df, aes(x=xax)) +
    geom_pointrange(aes(ymax=top, ymin=btm, y=mid, color=model)) +
    xlab(TeX("$\\sigma_A^2$")) +
    ylab(TeX("$h^2$")) +
    scale_x_log10() +
    scale_color_manual(
      name = "",
      values = c(
        "true" = "darkred", "h2liab" = "steelblue",
        "h2obs" = "chartreuse3"
      ),
      labels = c(
        expression("True " * h["liab"]^2),
        expression("Fitted " * h["liab"]^2),
        expression("Fitted " * h["obs"]^2)
      )
    ) +
    theme(text = element_text(size = 18), legend.text.align = 0)  
}

