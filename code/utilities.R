# Utilities

## Analyze things ran on Markov server
plot.markov.result <- function(rdata.file, ssh.samples, plot.title = NA, saveplot = F) {
  #' Analyze heritability from Rdata file. This is useful for binomial logit
  #' models or other models without closed form solutions
  load(rdata.file)
  df <- data.frame(
    Heritability = c(
      heritability_averaged$h2.obs,
      heritability_frequentist_avged$h2.obs,
      heritability_notavged,
      scaled.heritability
    ),
    Method = c(
      rep("h^2_Psi, Bayesian", ssh.samples),
      rep("h^2_Psi, Frequentist", ssh.samples),
      rep("h^2_Psi, No averaging", ssh.samples),
      rep("h^2_Phi", ssh.samples)
    )
  )
  color_palette <- c("#FFA600", "#FA2A00", "#7C43A6", "#0096FF")
  p <- ggplot(df, aes(x = Heritability, fill = Method)) +
    geom_density(alpha = 0.5) +
    ylab("Density") +
    scale_x_log10() +
    scale_fill_manual(
      values = color_palette,
      labels = c(
        expression(h[Phi]^2),
        expression(h[Psi]^2 * ", Bayesian"),
        expression(h[Psi]^2 * ", Frequentist"),
        expression(h[Psi]^2 * ", No averaging")
      )
    ) +
    theme(legend.text.align = 0)
  if (!is.na(plot.title)) p <- p + labs(title = plot.title)
  if (saveplot) {
    ggsave(
      filename = paste0(
        "../figures/",
        substring(rdata.file, 1, nchar(rdata.file) - 6),
        ".pdf"
      ),
      plot = p
    )
  }
  p
}


compute.observation.h2.SSH <- function(fit, nsamples, modelname, fn.prepend=""){
  #' Storing this momentarily to use when implementing runtime computations
  #' in main notebook. Otherwise not anything worth saving.
  cat("Called function at:", format(Sys.time(), "%a %b %d %X %Y"), "\n")
  t0 = Sys.time()
  latent.heritability = get.h2(fit, nsamples)
  scaled.heritability = get.h2(fit, nsamples, use.scale=T, model=modelname)
  
  ti = Sys.time()
  cat("Entering `new.h2.transf` (Bayesian sampling)\n")
  heritability_averaged = new.h2.transf(fit, modelname, nsamples)
  cat("Done after", difftime(Sys.time(),ti, units="mins"), "minutes.\nEntering `new.h2.transf` (Frequentist)\n")
  
  ti = Sys.time()
  marginal.mode = inla.posterior.marginal.latent.mode(fit)
  heritability_frequentist_avged = new.h2.transf(fit, modelname, nsamples, marginal.mode=marginal.mode)
  cat("Done after", difftime(Sys.time(),ti, units="mins"), "minutes. \nNow without averaging:\n")
  
  ti = Sys.time()
  heritability_notavged = get.h2.from.qgparams(fit, modelname, nsamples)
  cat("Done after", difftime(Sys.time(),ti, units="mins"), "minutes.\n")
  
  cat("Saving to disk...\n")
  save(latent.heritability, scaled.heritability, heritability_averaged,
       heritability_frequentist_avged, heritability_notavged,
       file=paste0(fn.prepend, "heritabilities_SSH_", modelname, ".Rdata"))
  return(0)
}


# Hyperparamter exploration for simulation case
simulation.hyperparameters <- function(NeNc_list = 0.1, sigmaA_list = 0.4) {
  #' Utility function for exploring hyperparameters in pedigree generation
  #' (Ne/Nc and sigmaA)
  p <- list()
  if (length(NeNc_list) > 1) {
    cat("Running simulations for NeNc values, sigmaA is fixed at", sigmaA_list[1])
    for (NeNc in NeNc_list) {
      simulation_results <- simulated.heritability(
        NeNc = NeNc, idgen = 100, nGen = 10, sigmaA = sigmaA_list,
        linear.predictor = function(u, .) u + rnorm(length(u)),
        simulated.formula = simulated.formula
      )
      p[[paste0(NeNc)]] <- simulation_results$heritability
    }
  } else {
    cat("Simulations for sigmaA, NeNc is fixed at", NeNc_list[1])
    for (sigmaA in sigmaA_list) {
      simulation_results <- simulated.heritability(
        NeNc = NeNc_list, idgen = 100, nGen = 10, sigmaA = sigmaA,
        linear.predictor = function(u, .) u + rnorm(length(u)),
        simulated.formula = simulated.formula
      )
      p[[paste0(sigmaA)]] <- simulation_results$heritability
    }
  }
  return(p)
}

hyperparameters_wrapper <- function() {
  # To avoid being ran every time...
  p <- simulation.hyperparameters(NeNc_list = c(0.05, 0.1, 0.25, 0.5))
  p2 <- simulation.hyperparameters(sigmaA_list = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75) * 100)
  
  
  par(mfrow = c(2, 2))
  for (i in names(p)) {
    truehist(get(i, p), xlab = paste("NeNc =", i))
    abline(v = 0.4)
  }
  par(mfrow = c(3, 2))
  for (i in names(p2)) {
    truehist(get(i, p2), xlab = paste("sigmaA =", i))
    cat("True h2", as.numeric(i) / (as.numeric(i) + 1), "\n")
    abline(v = as.numeric(i) / (as.numeric(i) + 1)) # a+e_i, e_i~N
  }
}

# Gaussian analysis - unsued plot
# <Plot 5> Recreation of QQ plot using residuals computed above instead of PIT values
#           Unused
plot.inla.qqplot <- function() {
  n.obs <- length(qg.data.gg.inds$u)
  qqplot(
    qnorm(ppoints(n.obs),
          mean = mean(qg.data.gg.inds$surv.ind.to.ad),
          sd = sd(qg.data.gg.inds$surv.ind.to.ad)
    ),
    resids,
    xlab = "Theoretical quantiles", ylab = "Sample Quantiles",
    main = "Q-Q plot from 'residuals' above"
  )
}