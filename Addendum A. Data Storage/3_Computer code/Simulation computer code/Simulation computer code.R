#################################################
########## Code Master Thesis Project  ##########
#################################################

### V number of the software. 
R.Version()

### Operation system of the machine.
#Windows 10

####################
##### Packages #####
####################

## "network" package.
# Check whether "network" package is installed and install if necessary.
if(!require(network)) install.packages("network")
# Push "network" package to library.
library(network)

## "ergm" package.
# Check whether "ergm" package is installed and install if necessary.
if(!require(ergm)) install.packages("ergm")
# Push "ergm" package to library.
library(ergm)

## "MASS" package.
# Check whether "MASS" package is installed and install if necessary.
if(!require(MASS)) install.packages("MASS")
# Push "MASS" package to library.
library(MASS)

## "metafor" package.
# Check whether "metafor" package is installed and install if necessary.
if(!require(metafor)) install.packages("metafor")
# Push "metafor" package to library.
library(metafor)

######################
### Pre-requisites ###
######################

### Load pre-imported dependencies.
load("myEnvironment.RData")

############################
##### Define functions #####
############################

timeR <- function(input){
  start.time <- Sys.time()
  input
  end.time <- Sys.time()
  (time.taken <- end.time - start.time)
}

resaveR <- function(..., list = character(), file) {
  previous  <- load(file)
  var.names <- c(list, as.character(substitute(list(...)))[-1L])
  for (var in var.names) assign(var, get(var, envir = parent.frame()))
  save(list = unique(c(previous, var.names)), file = file)
}

cov_sim_het0 <- function(nsim){
  
  ### Define objects / lists for storing function output. 
  ## Next counter.
  next_counter <- 0
  ## ERGM results. 
  ergm_results <- vector(mode = "list", length = nsim)
  ## ERGM fits. 
  ergm_fits <- vector(mode = "list", length = nsim)
  ## Lists for storing model statistics.
  t <- vector(mode = "list", length = nsim)
  mahalm <- vector(mode = "list", length = nsim)
  ## List for storing in-degree GOF. 
  tideg <- vector(mode = "list", length = nsim)
  ideg <- vector(mode = "list", length = nsim)
  ## List for storing out-degree GOF. 
  todeg <- vector(mode = "list", length = nsim)
  odeg <- vector(mode = "list", length = nsim)
  ## List for storing geodesic distance GOF. 
  tdist <- vector(mode = "list", length = nsim)
  dist <- vector(mode = "list", length = nsim)
  ## List for storing edgewise shared partners GOF. 
  tesp <- vector(mode = "list", length = nsim)
  esp <- vector(mode = "list", length = nsim)
  ## List for storing dyadwise shared partners GOF. 
  tdsp <- vector(mode = "list", length = nsim)
  dsp <- vector(mode = "list", length = nsim)
  ## List for storing triadic census GOF. 
  ttriad <- vector(mode = "list", length = nsim)
  triad <- vector(mode = "list", length = nsim)
  ## List for model parameter estimate differences relative to true model. 
  diffs <- vector(mode = "list", length = nsim)
  ## List for storing AIC and BIC.
  aic <- vector(mode = "list", length = nsim)
  bic <- vector(mode = "list", length = nsim)
  ## Seeds for reproducibility.
  seed1 <- as.integer(Sys.time()) %% 1000
  seed2 <- as.integer(Sys.time()) %% 10000 
  
  for (i in 1 : nsim){
    
    cat(paste0("\nCondition: param 'cov' @ size ", nsim, " @ het_0. Iteration ", i, ".\n")) 
    
    tryCatch({
      
      ### Fitting ERGM.
      
      ## Define model.
      model <- ergm_sims_het0[[i]] ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex")
      ## Fit model to "i"th simulated network.
      ergm_results[[i]] <- ergm(formula = model, control = control.ergm(MCMC.interval = 1 * 1024, MCMC.burnin = (1 * 1024) * 16, MCMC.samplesize = 1 * 1024, 
                                                                        seed = seed1))
      ## Calculate model and auxillary statistics fit of model. 
      ergm_fits[[i]] <- gof(object = model, GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus, 
                            coef = ergm_results[[i]]$coef, control.gof.formula(seed = seed2, nsim = 1000))
      
      ### Extracting GOF. 
      
      ## Model.
      # t-statistic.
      meant <- apply(ergm_fits[[i]]$sim.model, 2, mean)
      sdt <- apply(ergm_fits[[i]]$sim.model, 2, sd)
      obst <- summary(model)
      t[[i]] <- (meant-obst) / sdt
      # Mahalanobis.
      range <- which(sdt != 0)
      covm <- cov(as.matrix(ergm_fits[[i]]$sim.model[, range]))
      rangem <- range[which(eigen(covm)$values > 0.001)]
      covmodel <- cov(as.matrix(ergm_fits[[i]]$sim.model[, rangem]))
      mahalmodel <- mahalanobis(meant[rangem], ergm_fits[[i]]$obs.model[rangem], covmodel, tol = 1e-20)
      dfm <- length(rangem)
      pm <- 1 - pchisq(mahalmodel, dfm)
      mahalm[[i]] <- c(mahalmodel, dfm, pm)
      
      ## In-degree.
      # t-statistic.
      meantideg <- apply(ergm_fits[[i]]$sim.ideg, 2, mean)
      sdtideg <- apply(ergm_fits[[i]]$sim.ideg, 2, sd)
      tideg[[i]] <- (meantideg - ergm_fits[[i]]$obs.ideg) / sdtideg
      # Mahalanobis.
      range <- which(sdtideg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[,range])) 
      rangei <- range[which(eigen(covc)$values > 0.001)] 
      covideg <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[, rangei]))
      mahalideg <- mahalanobis(meantideg[rangei], ergm_fits[[i]]$obs.ideg[rangei], covideg, tol = 1e-20)
      dfi <- length(rangei)
      pi <- 1 - pchisq(mahalideg, dfi)
      ideg[[i]] <- c(mahalideg, dfi, pi)
      
      ## Out-degree.
      # t-statistic.
      meantodeg <- apply(ergm_fits[[i]]$sim.odeg,2,mean)
      sdtodeg <- apply(ergm_fits[[i]]$sim.odeg,2,sd)
      todeg[[i]] <- (meantodeg-ergm_fits[[i]]$obs.odeg)/sdtodeg
      # Mahalanobis.
      range <- which(sdtodeg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, range]))
      rangeo <- range[which(eigen(covc)$values > 0.001)] 
      covodeg <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, rangeo]))
      mahalodeg <- mahalanobis(meantodeg[rangeo], ergm_fits[[i]]$obs.odeg[rangeo], covodeg, tol = 1e-20)
      dfo <- length(rangeo)
      po <- 1 - pchisq(mahalodeg, dfo)
      odeg[[i]] <- c(mahalodeg, dfo, po)
      
      ## Geodesic distance.
      # t-statistic.
      meantdist <- apply(ergm_fits[[i]]$sim.dist, 2, mean)
      sdtdist <- apply(ergm_fits[[i]]$sim.dist, 2, sd)
      tdist[[i]] <- (meantdist - ergm_fits[[i]]$obs.dist) / sdtdist
      # Mahalanobis.
      range <- which(sdtdist != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, range]))
      ranged <- range[which(eigen(covc)$values > 0.001)] 
      covdist <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, ranged]))
      mahaldist <- mahalanobis(meantdist[ranged], ergm_fits[[i]]$obs.dist[ranged], covdist, tol = 1e-20)
      dfdist <- length(ranged)
      pdist <- 1 - pchisq(mahaldist, dfdist)
      dist[[i]] <- c(mahaldist,dfdist,pdist)
      
      ## Edgewise shared partners.
      # t-statistic.
      meantesp <- apply(ergm_fits[[i]]$sim.esp, 2, mean)
      sdtesp <- apply(ergm_fits[[i]]$sim.esp, 2, sd)
      tesp[[i]] <- (meantesp - ergm_fits[[i]]$obs.esp) / sdtesp
      # Mahalanobis.
      range <- which(sdtesp != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, range]))
      rangee <- range[which(eigen(covc)$values > 0.001)] 
      covesp <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, rangee]))
      mahalesp <- mahalanobis(meantesp[rangee], ergm_fits[[i]]$obs.esp[rangee], covesp, tol = 1e-20)
      dfesp <- length(rangee)
      pesp <- 1 - pchisq(mahalesp, dfesp)
      esp[[i]] <- c(mahalesp, dfesp, pesp)
      
      ## Dyadwise shared partners.
      # t-statistic.
      meantdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, mean)
      sdtdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, sd)
      tdsp[[i]] <- c(meantdsp - ergm_fits[[i]]$obs.dsp) / sdtdsp
      # Mahalanobis.
      range <- which(sdtdsp != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, range]))
      rangeds <- range[which(eigen(covc)$values > 0.001)] 
      covdsp <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, rangeds]))
      mahaldsp <- mahalanobis(meantdsp[rangeds], ergm_fits[[i]]$obs.dsp[rangeds], covdsp, tol = 1e-20)
      dfdsp <- length(rangeds)
      pdsp <- 1 - pchisq(mahaldsp, dfdsp)
      dsp[[i]] <- c(mahaldsp, dfdsp, pdsp)
      
      ## Triadcensus.
      # t-statistic.
      meanttriad <- apply(ergm_fits[[i]]$sim.triad, 2, mean)
      sdttriad <- apply(ergm_fits[[i]]$sim.triad, 2, sd)
      ttriad[[i]] <- (meanttriad - ergm_fits[[i]]$obs.triad) / sdttriad
      # Mahalanobis.
      range <- which(sdttriad != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, range]))
      ranget <- range[which(eigen(covc)$values > 0.001)] 
      covtriad <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, ranget]))
      mahaltriad <- mahalanobis(meanttriad[ranget], ergm_fits[[i]]$obs.triad[ranget], covtriad, tol = 1e-20)
      dftriad <- length(ranget)
      ptriad <- 1 - pchisq(mahaltriad, dftriad)
      triad[[i]] <- c(mahaltriad, dftriad, ptriad)
      
      ## Bias and rmse.
      # Differences.
      diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het0$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")]   
      
      ## AIC / BIC.
      # AIC.
      aic[i] <- summary(ergm_results[[i]])[17]  
      # BIC.
      bic[i] <- summary(ergm_results[[i]])[18]
      
      retry <<- FALSE
      
    }, error = function(e){
      
      cat(paste0("Condition: param 'cov' @ size ", nsim, " @ het_0. Error in iteration ", i, ".1. Re-simulating and re-fitting.\n"))
      next_counter <<- next_counter + 1
      retry <<- TRUE
      
    })
    
    if(retry){
      
      k <- 1
      err <<- TRUE
      
      while(isTRUE(err)){
        
        tryCatch({
          
          k <- k + 1
          
          ### Resimulate.
          ergm_sims_het0[[i]] <<- simulate(true_ergm_results_het0, nsim = 1, set.seed = seed1)
          ### Define model.
          model <- ergm_sims_het0[[i]] ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex")
          ### Fit model to "i"th simulated network.
          ergm_results[[i]] <- ergm(formula = model, control = control.ergm(MCMC.interval = 1 * 1024, MCMC.burnin = (1 * 1024) * 16, MCMC.samplesize = 1 * 1024, 
                                                                            seed = seed1))
          ### Calculate model and auxillary statistics fit of model. 
          ergm_fits[[i]] <- gof(object = model, GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus, 
                                coef = ergm_results[[i]]$coef, control.gof.formula(seed = seed2, nsim = 1000))
          
          ### GOF. 
          
          ## Model.
          # t-statistic.
          meant <- apply(ergm_fits[[i]]$sim.model, 2, mean)
          sdt <- apply(ergm_fits[[i]]$sim.model, 2, sd)
          obst <- summary(model)
          t[[i]] <- (meant-obst) / sdt
          # Mahalanobis.
          range <- which(sdt != 0)
          covm <- cov(as.matrix(ergm_fits[[i]]$sim.model[, range]))
          rangem <- range[which(eigen(covm)$values > 0.001)]
          covmodel <- cov(as.matrix(ergm_fits[[i]]$sim.model[, rangem]))
          mahalmodel <- mahalanobis(meant[rangem], ergm_fits[[i]]$obs.model[rangem], covmodel, tol = 1e-20)
          dfm <- length(rangem)
          pm <- 1 - pchisq(mahalmodel, dfm)
          mahalm[[i]] <- c(mahalmodel, dfm, pm)
          
          ## In-degree.
          # t-statistic.
          meantideg <- apply(ergm_fits[[i]]$sim.ideg, 2, mean)
          sdtideg <- apply(ergm_fits[[i]]$sim.ideg, 2, sd)
          tideg[[i]] <- (meantideg - ergm_fits[[i]]$obs.ideg) / sdtideg
          # Mahalanobis.
          range <- which(sdtideg != 0)
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[,range])) 
          rangei <- range[which(eigen(covc)$values > 0.001)] 
          covideg <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[, rangei]))
          mahalideg <- mahalanobis(meantideg[rangei], ergm_fits[[i]]$obs.ideg[rangei], covideg, tol = 1e-20)
          dfi <- length(rangei)
          pi <- 1 - pchisq(mahalideg, dfi)
          ideg[[i]] <- c(mahalideg, dfi, pi)
          
          ## Out-degree.
          # t-statistic.
          meantodeg <- apply(ergm_fits[[i]]$sim.odeg,2,mean)
          sdtodeg <- apply(ergm_fits[[i]]$sim.odeg,2,sd)
          todeg[[i]] <- (meantodeg-ergm_fits[[i]]$obs.odeg)/sdtodeg
          # Mahalanobis.
          range <- which(sdtodeg != 0)
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, range]))
          rangeo <- range[which(eigen(covc)$values > 0.001)] 
          covodeg <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, rangeo]))
          mahalodeg <- mahalanobis(meantodeg[rangeo], ergm_fits[[i]]$obs.odeg[rangeo], covodeg, tol = 1e-20)
          dfo <- length(rangeo)
          po <- 1 - pchisq(mahalodeg, dfo)
          odeg[[i]] <- c(mahalodeg, dfo, po)
          
          ## Geodesic distance.
          # t-statistic.
          meantdist <- apply(ergm_fits[[i]]$sim.dist, 2, mean)
          sdtdist <- apply(ergm_fits[[i]]$sim.dist, 2, sd)
          tdist[[i]] <- (meantdist - ergm_fits[[i]]$obs.dist) / sdtdist
          # Mahalanobis.
          range <- which(sdtdist != 0)
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, range]))
          ranged <- range[which(eigen(covc)$values > 0.001)] 
          covdist <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, ranged]))
          mahaldist <- mahalanobis(meantdist[ranged], ergm_fits[[i]]$obs.dist[ranged], covdist, tol = 1e-20)
          dfdist <- length(ranged)
          pdist <- 1 - pchisq(mahaldist, dfdist)
          dist[[i]] <- c(mahaldist,dfdist,pdist)
          
          ## Edgewise shared partners.
          # t-statistic.
          meantesp <- apply(ergm_fits[[i]]$sim.esp, 2, mean)
          sdtesp <- apply(ergm_fits[[i]]$sim.esp, 2, sd)
          tesp[[i]] <- (meantesp - ergm_fits[[i]]$obs.esp) / sdtesp
          # Mahalanobis.
          range <- which(sdtesp != 0) 
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, range]))
          rangee <- range[which(eigen(covc)$values > 0.001)] 
          covesp <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, rangee]))
          mahalesp <- mahalanobis(meantesp[rangee], ergm_fits[[i]]$obs.esp[rangee], covesp, tol = 1e-20)
          dfesp <- length(rangee)
          pesp <- 1 - pchisq(mahalesp, dfesp)
          esp[[i]] <- c(mahalesp, dfesp, pesp)
          
          ## Dyadwise shared partners.
          # t-statistic.
          meantdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, mean)
          sdtdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, sd)
          tdsp[[i]] <- c(meantdsp - ergm_fits[[i]]$obs.dsp) / sdtdsp
          # Mahalanobis.
          range <- which(sdtdsp != 0)
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, range]))
          rangeds <- range[which(eigen(covc)$values > 0.001)] 
          covdsp <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, rangeds]))
          mahaldsp <- mahalanobis(meantdsp[rangeds], ergm_fits[[i]]$obs.dsp[rangeds], covdsp, tol = 1e-20)
          dfdsp <- length(rangeds)
          pdsp <- 1 - pchisq(mahaldsp, dfdsp)
          dsp[[i]] <- c(mahaldsp, dfdsp, pdsp)
          
          ## Triadcensus.
          # t-statistic.
          meanttriad <- apply(ergm_fits[[i]]$sim.triad, 2, mean)
          sdttriad <- apply(ergm_fits[[i]]$sim.triad, 2, sd)
          ttriad[[i]] <- (meanttriad - ergm_fits[[i]]$obs.triad) / sdttriad
          # Mahalanobis.
          range <- which(sdttriad != 0) 
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, range]))
          ranget <- range[which(eigen(covc)$values > 0.001)] 
          covtriad <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, ranget]))
          mahaltriad <- mahalanobis(meanttriad[ranget], ergm_fits[[i]]$obs.triad[ranget], covtriad, tol = 1e-20)
          dftriad <- length(ranget)
          ptriad <- 1 - pchisq(mahaltriad, dftriad)
          triad[[i]] <- c(mahaltriad, dftriad, ptriad)
          
          ## Bias and rmse.
          # Differences.
          diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het0$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")]   
          
          ## AIC / BIC.
          # AIC.
          aic[i] <- summary(ergm_results[[i]])[17]  
          # BIC.
          bic[i] <- summary(ergm_results[[i]])[18]
          
          err <<- FALSE
          
        }, error = function(e){
          
          cat(paste0("Condition: param 'cov' @ size ", nsim, " @ het_0. Error in iteration ", i, ".", k, ". Re-simulating and re-fitting.\n"))
          err <<- TRUE
          
        })
        
      }
      
    }
    
  }
  
  ### Meta-analysis.
  ma_coeffs <- matrix(nrow = nsim, ncol = length(diffs[[1]]))
  ma_ses <- matrix(nrow = nsim, ncol = length(diffs[[1]]))
  ma_model <- list()
  for (i in 1 : nsim){
    for (j in 1 : length(diffs[[1]])){
      ma_coeffs[i, ] <- ergm_results[[i]]$coef[1 : length(diffs[[1]])]
      ma_ses[i, ] <- sqrt(diag(ergm_results[[i]]$covar))
      suppressWarnings(ma_model[[j]] <- rma(yi = ma_coeffs[, j], sei = ma_ses[, j], control = list(stepadj = 0.5, maxiter = 1000)))
    }
  }
  
  ### Bias and rmse. 
  ## Bias.
  bias_matrix <- matrix(data = NA, nrow = nsim, ncol = length(diffs[[1]]))
  for (i in 1 : nsim){
    bias_matrix[i, ] <- diffs[[i]]
  }
  bias <- colSums(bias_matrix) / nsim
  
  ## Rmse.
  rmse <- sqrt(colSums(bias_matrix ** 2) / nsim)
  
  ### Mean AIC and BIC.
  mean_aic <- mean(unlist(aic)) # Mean "meta-analysis" AIC.
  mean_bic <- mean(unlist(bic)) # Mean "meta-analysis" BIC.
  
  return(list(next_counter = next_counter,
              model_GOF = list(t = t, mahalm = mahalm), indegree_aux = list(tideg = tideg, ideg = ideg), 
              outdegree_aux = list(todeg = todeg, odeg = odeg), dist_aux = list(tdist = tdist, dist = dist), 
              esp_aux = list(tesp = tesp, esp = esp), dsp_aux = list(tdsp = tdsp, dsp = dsp), 
              triad_aux = list(ttriad = ttriad, triad = triad), aic = list(aic = unlist(aic), mean_aic = mean_aic), 
              bic = list(bic = unlist(bic), mean_bic = mean_bic), ma_model = ma_model, 
              bias = bias, rmse = rmse, seed1 = seed1, seed2 = seed2))
  
}

gwodegree_sim_het0 <- function(nsim){
  
  ### Define objects / lists for storing function output. 
  ## Flag.
  flag <- TRUE
  ## Next counter.
  next_counter <- 0
  ## ERGM results. 
  ergm_results <- vector(mode = "list", length = nsim)
  ## ERGM fits. 
  ergm_fits <- vector(mode = "list", length = nsim)
  ## List for storing model statistics. 
  t <- vector(mode = "list", length = nsim)
  mahalm <- vector(mode = "list", length = nsim)
  ## List for storing in-degree GOF. 
  tideg <- vector(mode = "list", length = nsim)
  ideg <- vector(mode = "list", length = nsim)
  ## List for storing out-degree GOF. 
  todeg <- vector(mode = "list", length = nsim)
  odeg <- vector(mode = "list", length = nsim)
  ## List for storing geodesic distance GOF. 
  tdist <- vector(mode = "list", length = nsim)
  dist <- vector(mode = "list", length = nsim)
  ## List for storing edgewise shared partners GOF. 
  tesp <- vector(mode = "list", length = nsim)
  esp <- vector(mode = "list", length = nsim)
  ## List for storing dyadwise shared partners GOF. 
  tdsp <- vector(mode = "list", length = nsim)
  dsp <- vector(mode = "list", length = nsim)
  ## List for storing triadic census GOF. 
  ttriad <- vector(mode = "list", length = nsim)
  triad <- vector(mode = "list", length = nsim)
  ## List for model parameter estimate differences relative to true model. 
  diffs <- vector(mode = "list", length = nsim)
  ## List for storing AIC and BIC.
  aic <- vector(mode = "list", length = nsim)
  bic <- vector(mode = "list", length = nsim)
  ## Seeds for reproducibility.
  seed1 <- as.integer(Sys.time()) %% 1000
  seed2 <- as.integer(Sys.time()) %% 10000 
  
  for (i in 1 : nsim){
    
    cat(paste0("\nCondition: param 'gwodeg' @ size ", nsim, " @ het_0. Iteration ", i, ".\n"))    
    
    tryCatch({
      
      ### Define model.
      model <- ergm_sims_het0[[i]] ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T)
      ### Fit model to "i"th simulated network.
      ergm_results[[i]] <- ergm(formula = model, control = control.ergm(MCMC.interval = 1 * 1024, MCMC.burnin = (1 * 1024) * 16, MCMC.samplesize = 1 * 1024, 
                                                                        seed = seed1))
      ### Calculate model and auxillary statistics fit of model. 
      ergm_fits[[i]] <- gof(object = model, GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus, 
                            coef = ergm_results[[i]]$coef, control.gof.formula(seed = seed2, nsim = 1000))
      
      ### GOF. 
      
      ## Model.
      # t-statistic.
      meant <- apply(ergm_fits[[i]]$sim.model, 2, mean)
      sdt <- apply(ergm_fits[[i]]$sim.model, 2, sd)
      obst <- summary(model)
      t[[i]] <- (meant-obst) / sdt
      # Mahalanobis.
      range <- which(sdt != 0)
      covm <- cov(as.matrix(ergm_fits[[i]]$sim.model[, range]))
      rangem <- range[which(eigen(covm)$values > 0.001)]
      covmodel <- cov(as.matrix(ergm_fits[[i]]$sim.model[, rangem]))
      mahalmodel <- mahalanobis(meant[rangem], ergm_fits[[i]]$obs.model[rangem], covmodel, tol = 1e-20)
      dfm <- length(rangem)
      pm <- 1 - pchisq(mahalmodel, dfm)
      mahalm[[i]] <- c(mahalmodel, dfm, pm)
      
      ## In-degree.
      # t-statistic.
      meantideg <- apply(ergm_fits[[i]]$sim.ideg, 2, mean)
      sdtideg <- apply(ergm_fits[[i]]$sim.ideg, 2, sd)
      tideg[[i]] <- (meantideg - ergm_fits[[i]]$obs.ideg) / sdtideg
      # Mahalanobis.
      range <- which(sdtideg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[,range])) 
      rangei <- range[which(eigen(covc)$values > 0.001)] 
      covideg <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[, rangei]))
      mahalideg <- mahalanobis(meantideg[rangei], ergm_fits[[i]]$obs.ideg[rangei], covideg, tol = 1e-20)
      dfi <- length(rangei)
      pi <- 1 - pchisq(mahalideg, dfi)
      ideg[[i]] <- c(mahalideg, dfi, pi)
      
      ## Out-degree.
      # t-statistic.
      meantodeg <- apply(ergm_fits[[i]]$sim.odeg,2,mean)
      sdtodeg <- apply(ergm_fits[[i]]$sim.odeg,2,sd)
      todeg[[i]] <- (meantodeg-ergm_fits[[i]]$obs.odeg)/sdtodeg
      # Mahalanobis.
      range <- which(sdtodeg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, range]))
      rangeo <- range[which(eigen(covc)$values > 0.001)] 
      covodeg <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, rangeo]))
      mahalodeg <- mahalanobis(meantodeg[rangeo], ergm_fits[[i]]$obs.odeg[rangeo], covodeg, tol = 1e-20)
      dfo <- length(rangeo)
      po <- 1 - pchisq(mahalodeg, dfo)
      odeg[[i]] <- c(mahalodeg, dfo, po)
      
      ## Geodesic distance.
      # t-statistic.
      meantdist <- apply(ergm_fits[[i]]$sim.dist, 2, mean)
      sdtdist <- apply(ergm_fits[[i]]$sim.dist, 2, sd)
      tdist[[i]] <- (meantdist - ergm_fits[[i]]$obs.dist) / sdtdist
      # Mahalanobis.
      range <- which(sdtdist != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, range]))
      ranged <- range[which(eigen(covc)$values > 0.001)] 
      covdist <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, ranged]))
      mahaldist <- mahalanobis(meantdist[ranged], ergm_fits[[i]]$obs.dist[ranged], covdist, tol = 1e-20)
      dfdist <- length(ranged)
      pdist <- 1 - pchisq(mahaldist, dfdist)
      dist[[i]] <- c(mahaldist,dfdist,pdist)
      
      ## Edgewise shared partners.
      # t-statistic.
      meantesp <- apply(ergm_fits[[i]]$sim.esp, 2, mean)
      sdtesp <- apply(ergm_fits[[i]]$sim.esp, 2, sd)
      tesp[[i]] <- (meantesp - ergm_fits[[i]]$obs.esp) / sdtesp
      # Mahalanobis.
      range <- which(sdtesp != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, range]))
      rangee <- range[which(eigen(covc)$values > 0.001)] 
      covesp <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, rangee]))
      mahalesp <- mahalanobis(meantesp[rangee], ergm_fits[[i]]$obs.esp[rangee], covesp, tol = 1e-20)
      dfesp <- length(rangee)
      pesp <- 1 - pchisq(mahalesp, dfesp)
      esp[[i]] <- c(mahalesp, dfesp, pesp)
      
      ## Dyadwise shared partners.
      # t-statistic.
      meantdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, mean)
      sdtdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, sd)
      tdsp[[i]] <- c(meantdsp - ergm_fits[[i]]$obs.dsp) / sdtdsp
      # Mahalanobis.
      range <- which(sdtdsp != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, range]))
      rangeds <- range[which(eigen(covc)$values > 0.001)] 
      covdsp <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, rangeds]))
      mahaldsp <- mahalanobis(meantdsp[rangeds], ergm_fits[[i]]$obs.dsp[rangeds], covdsp, tol = 1e-20)
      dfdsp <- length(rangeds)
      pdsp <- 1 - pchisq(mahaldsp, dfdsp)
      dsp[[i]] <- c(mahaldsp, dfdsp, pdsp)
      
      ## Triadcensus.
      # t-statistic.
      meanttriad <- apply(ergm_fits[[i]]$sim.triad, 2, mean)
      sdttriad <- apply(ergm_fits[[i]]$sim.triad, 2, sd)
      ttriad[[i]] <- (meanttriad - ergm_fits[[i]]$obs.triad) / sdttriad
      # Mahalanobis.
      range <- which(sdttriad != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, range]))
      ranget <- range[which(eigen(covc)$values > 0.001)] 
      covtriad <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, ranget]))
      mahaltriad <- mahalanobis(meanttriad[ranget], ergm_fits[[i]]$obs.triad[ranget], covtriad, tol = 1e-20)
      dftriad <- length(ranget)
      ptriad <- 1 - pchisq(mahaltriad, dftriad)
      triad[[i]] <- c(mahaltriad, dftriad, ptriad)
      
      ## Bias and rmse.
      # Differences.
      diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1")] - true_ergm_results_het0$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1")]   
      ## AIC / BIC.
      # AIC.
      aic[i] <- summary(ergm_results[[i]])[17]  
      # BIC.
      bic[i] <- summary(ergm_results[[i]])[18]
      
    }, error = function(e){
      
      cat(paste0("Condition: param 'gwodeg' @ size ", nsim, " @ het_0. Error in iteration ", i, ". Skipping to next iteration.\n"))
      next_counter <<- next_counter + 1
      flag <<- FALSE
      
      ergm_results[i] <<- list(NULL) 
      ergm_fits[i] <<- list(NULL) 
      t[i] <<- list(NULL) 
      mahalm[i] <<- list(NULL) 
      tideg[i] <<- list(NULL) 
      ideg[i] <<- list(NULL) 
      todeg[i] <<- list(NULL) 
      odeg[i] <<- list(NULL) 
      tdist[i] <<- list(NULL) 
      dist[i] <<- list(NULL) 
      tesp[i] <<- list(NULL) 
      esp[i] <<- list(NULL) 
      tdsp[i] <<- list(NULL) 
      dsp[i] <<- list(NULL) 
      ttriad[i] <<- list(NULL) 
      triad[i] <<- list(NULL) 
      diffs[i] <<- list(NULL) 
      aic[i] <<- NA
      bic[i] <<- NA  
      
    })
    
    if (!flag) next
    
  }
  
  print(ergm_results)
  print(diffs)
  ergm_results <- ergm_results[!sapply(ergm_results, is.null)]
  ergm_fits <- ergm_fits[!sapply(ergm_fits, is.null)]
  t <- t[!sapply(t, is.null)]
  mahalm <- mahalm[!sapply(mahalm, is.null)]
  tideg <- tideg[!sapply(tideg, is.null)]
  ideg <- ideg[!sapply(ideg, is.null)]
  todeg <- todeg[!sapply(todeg, is.null)]
  odeg <- odeg[!sapply(odeg, is.null)]
  tdist <- tdist[!sapply(tdist, is.null)]
  dist <- dist[!sapply(dist, is.null)]
  tesp <- tesp[!sapply(tesp, is.null)]
  esp <- esp[!sapply(esp, is.null)]
  tdsp <- tdsp[!sapply(tdsp, is.null)]
  dsp <- dsp[!sapply(dsp, is.null)]
  ttriad <- ttriad[!sapply(ttriad, is.null)]
  triad <- triad[!sapply(triad, is.null)]
  diffs <- diffs[!sapply(diffs, is.null)]
  aic <- aic[!is.na(aic)]
  bic <- bic[!is.na(bic)]
  print(ergm_results)
  print(diffs)
  print(length(diffs))
  
  ### Meta-analysis.
  ma_coeffs <- matrix(nrow = length(diffs), ncol = length(diffs[[1]]))
  ma_ses <- matrix(nrow = length(diffs), ncol = length(diffs[[1]]))
  ma_model <- list()
  for (i in 1 : length(diffs)){
    for (j in 1 : length(diffs[[1]])){
      ma_coeffs[i, ] <- ergm_results[[i]]$coef[1 : length(diffs[[1]])]
      ma_ses[i, ] <- sqrt(diag(ergm_results[[i]]$covar[1 : length(diffs[[1]]), 1 : length(diffs[[1]])]))
      suppressWarnings(ma_model[[j]] <- rma(yi = ma_coeffs[, j], sei = ma_ses[, j], control = list(stepadj = 0.5, maxiter = 1000)))
    }
  }
  
  ### Bias and rmse. 
  ## Bias.
  bias_matrix <- matrix(data = NA, nrow = length(diffs), ncol = length(diffs[[1]]))
  for (i in 1 : length(diffs)){
    bias_matrix[i, ] <- diffs[[i]]
  }
  bias <- colSums(bias_matrix) / length(diffs)
  
  ## rmse.
  rmse <- sqrt(colSums(bias_matrix ** 2) / length(diffs))
  
  ### Mean AIC and BIC.
  mean_aic <- mean(unlist(aic)) # Mean "meta-analysis" AIC.
  mean_bic <- mean(unlist(bic)) # Mean "meta-analysis" BIC.
  
  return(list(next_counter = next_counter,
              model_GOF = list(t = t, mahalm = mahalm), indegree_aux = list(tideg = tideg, ideg = ideg), 
              outdegree_aux = list(todeg = todeg, odeg = odeg), dist_aux = list(tdist = tdist, dist = dist), 
              esp_aux = list(tesp = tesp, esp = esp), dsp_aux = list(tdsp = tdsp, dsp = dsp), 
              triad_aux = list(ttriad = ttriad, triad = triad), aic = list(aic = unlist(aic), mean_aic = mean_aic), 
              bic = list(bic = unlist(bic), mean_bic = mean_bic), ma_model = ma_model,
              bias = bias, rmse = rmse, seed1 = seed1, seed2 = seed2))
  
}

dtriad_sim_het0 <- function(nsim){
  
  ### Define objects / lists for storing function output. 
  ## Flag.
  flag <- TRUE
  ## Next counter.
  next_counter <- 0
  ## ERGM results. 
  ergm_results <- vector(mode = "list", length = nsim)
  ## ERGM fits. 
  ergm_fits <- vector(mode = "list", length = nsim)
  ## List for storing model statistics. 
  t <- vector(mode = "list", length = nsim)
  mahalm <- vector(mode = "list", length = nsim)
  ## List for storing in-degree GOF. 
  tideg <- vector(mode = "list", length = nsim)
  ideg <- vector(mode = "list", length = nsim)
  ## List for storing out-degree GOF. 
  todeg <- vector(mode = "list", length = nsim)
  odeg <- vector(mode = "list", length = nsim)
  ## List for storing geodesic distance GOF. 
  tdist <- vector(mode = "list", length = nsim)
  dist <- vector(mode = "list", length = nsim)
  ## List for storing edgewise shared partners GOF. 
  tesp <- vector(mode = "list", length = nsim)
  esp <- vector(mode = "list", length = nsim)
  ## List for storing dyadwise shared partners GOF. 
  tdsp <- vector(mode = "list", length = nsim)
  dsp <- vector(mode = "list", length = nsim)
  ## List for storing triadic census GOF. 
  ttriad <- vector(mode = "list", length = nsim)
  triad <- vector(mode = "list", length = nsim)
  ## List for model parameter estimate differences relative to true model. 
  diffs <- vector(mode = "list", length = nsim)
  ## List for storing AIC and BIC.
  aic <- vector(mode = "list", length = nsim)
  bic <- vector(mode = "list", length = nsim)
  ## Seeds for reproducibility.
  seed1 <- as.integer(Sys.time()) %% 1000
  seed2 <- as.integer(Sys.time()) %% 10000 
  
  for (i in 1 : nsim){
    
    cat(paste0("\nCondition: param 'dtriad' @ size ", nsim, " @ het_0. Iteration ", i, ".\n"))   
    
    tryCatch({
      
      ### Define model.
      model <- ergm_sims_het0[[i]] ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex") + triadcensus(14)
      ### Fit model to "i"th simulated network.
      ergm_results[[i]] <- ergm(formula = model, control = control.ergm(MCMC.interval = 1 * 1024, MCMC.burnin = (1 * 1024) * 16, MCMC.samplesize = 1 * 1024, 
                                                                        seed = seed1))
      ### Calculate model and auxillary statistics fit of model. 
      ergm_fits[[i]] <- gof(object = model, GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus, 
                            coef = ergm_results[[i]]$coef, control.gof.formula(seed = seed2, nsim = 1000))
      
      if (any(is.infinite(ergm_results[[i]]$coef))){
        
        ergm_results[i] <- list(NULL) 
        ergm_fits[i] <- list(NULL) 
        stop()
        
      } 
      
      ### GOF. 
      
      ## Model.
      # t-statistic.
      meant <- apply(ergm_fits[[i]]$sim.model, 2, mean)
      sdt <- apply(ergm_fits[[i]]$sim.model, 2, sd)
      obst <- summary(model)
      t[[i]] <- (meant-obst) / sdt
      # Mahalanobis.
      range <- which(sdt != 0)
      covm <- cov(as.matrix(ergm_fits[[i]]$sim.model[, range]))
      rangem <- range[which(eigen(covm)$values > 0.001)]
      covmodel <- cov(as.matrix(ergm_fits[[i]]$sim.model[, rangem]))
      mahalmodel <- mahalanobis(meant[rangem], ergm_fits[[i]]$obs.model[rangem], covmodel, tol = 1e-20)
      dfm <- length(rangem)
      pm <- 1 - pchisq(mahalmodel, dfm)
      mahalm[[i]] <- c(mahalmodel, dfm, pm)
      
      ## In-degree.
      # t-statistic.
      meantideg <- apply(ergm_fits[[i]]$sim.ideg, 2, mean)
      sdtideg <- apply(ergm_fits[[i]]$sim.ideg, 2, sd)
      tideg[[i]] <- (meantideg - ergm_fits[[i]]$obs.ideg) / sdtideg
      # Mahalanobis.
      range <- which(sdtideg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[,range])) 
      rangei <- range[which(eigen(covc)$values > 0.001)] 
      covideg <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[, rangei]))
      mahalideg <- mahalanobis(meantideg[rangei], ergm_fits[[i]]$obs.ideg[rangei], covideg, tol = 1e-20)
      dfi <- length(rangei)
      pi <- 1 - pchisq(mahalideg, dfi)
      ideg[[i]] <- c(mahalideg, dfi, pi)
      
      ## Out-degree.
      # t-statistic.
      meantodeg <- apply(ergm_fits[[i]]$sim.odeg,2,mean)
      sdtodeg <- apply(ergm_fits[[i]]$sim.odeg,2,sd)
      todeg[[i]] <- (meantodeg-ergm_fits[[i]]$obs.odeg)/sdtodeg
      # Mahalanobis.
      range <- which(sdtodeg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, range]))
      rangeo <- range[which(eigen(covc)$values > 0.001)] 
      covodeg <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, rangeo]))
      mahalodeg <- mahalanobis(meantodeg[rangeo], ergm_fits[[i]]$obs.odeg[rangeo], covodeg, tol = 1e-20)
      dfo <- length(rangeo)
      po <- 1 - pchisq(mahalodeg, dfo)
      odeg[[i]] <- c(mahalodeg, dfo, po)
      
      ## Geodesic distance.
      # t-statistic.
      meantdist <- apply(ergm_fits[[i]]$sim.dist, 2, mean)
      sdtdist <- apply(ergm_fits[[i]]$sim.dist, 2, sd)
      tdist[[i]] <- (meantdist - ergm_fits[[i]]$obs.dist) / sdtdist
      # Mahalanobis.
      range <- which(sdtdist != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, range]))
      ranged <- range[which(eigen(covc)$values > 0.001)] 
      covdist <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, ranged]))
      mahaldist <- mahalanobis(meantdist[ranged], ergm_fits[[i]]$obs.dist[ranged], covdist, tol = 1e-20)
      dfdist <- length(ranged)
      pdist <- 1 - pchisq(mahaldist, dfdist)
      dist[[i]] <- c(mahaldist,dfdist,pdist)
      
      ## Edgewise shared partners.
      # t-statistic.
      meantesp <- apply(ergm_fits[[i]]$sim.esp, 2, mean)
      sdtesp <- apply(ergm_fits[[i]]$sim.esp, 2, sd)
      tesp[[i]] <- (meantesp - ergm_fits[[i]]$obs.esp) / sdtesp
      # Mahalanobis.
      range <- which(sdtesp != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, range]))
      rangee <- range[which(eigen(covc)$values > 0.001)] 
      covesp <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, rangee]))
      mahalesp <- mahalanobis(meantesp[rangee], ergm_fits[[i]]$obs.esp[rangee], covesp, tol = 1e-20)
      dfesp <- length(rangee)
      pesp <- 1 - pchisq(mahalesp, dfesp)
      esp[[i]] <- c(mahalesp, dfesp, pesp)
      
      ## Dyadwise shared partners.
      # t-statistic.
      meantdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, mean)
      sdtdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, sd)
      tdsp[[i]] <- c(meantdsp - ergm_fits[[i]]$obs.dsp) / sdtdsp
      # Mahalanobis.
      range <- which(sdtdsp != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, range]))
      rangeds <- range[which(eigen(covc)$values > 0.001)] 
      covdsp <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, rangeds]))
      mahaldsp <- mahalanobis(meantdsp[rangeds], ergm_fits[[i]]$obs.dsp[rangeds], covdsp, tol = 1e-20)
      dfdsp <- length(rangeds)
      pdsp <- 1 - pchisq(mahaldsp, dfdsp)
      dsp[[i]] <- c(mahaldsp, dfdsp, pdsp)
      
      ## Triadcensus.
      # t-statistic.
      meanttriad <- apply(ergm_fits[[i]]$sim.triad, 2, mean)
      sdttriad <- apply(ergm_fits[[i]]$sim.triad, 2, sd)
      ttriad[[i]] <- (meanttriad - ergm_fits[[i]]$obs.triad) / sdttriad
      # Mahalanobis.
      range <- which(sdttriad != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, range]))
      ranget <- range[which(eigen(covc)$values > 0.001)] 
      covtriad <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, ranget]))
      mahaltriad <- mahalanobis(meanttriad[ranget], ergm_fits[[i]]$obs.triad[ranget], covtriad, tol = 1e-20)
      dftriad <- length(ranget)
      ptriad <- 1 - pchisq(mahaltriad, dftriad)
      triad[[i]] <- c(mahaltriad, dftriad, ptriad)
      ## Bias and rmse.
      # Differences.
      diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het0$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] 
      ## AIC / BIC.
      # AIC.
      aic[i] <- summary(ergm_results[[i]])[17]  
      # BIC.
      bic[i] <- summary(ergm_results[[i]])[18]
      
    }, error = function(e){
      
      if (e$message == ""){
        
        cat(paste0("Condition: param 'dtriad' @ size ", nsim, " @ het_0. Infinite coefficient in iteration ", i,". Skipping to next iteration.\n")) 
        next_counter <<- next_counter + 1
        flag <<- FALSE
        
      } else{
        
        cat(paste0("Condition: param 'dtriad' @ size ", nsim, " @ het_0. Error in iteration ", i, ". Skipping to next iteration.\n"))
        next_counter <<- next_counter + 1
        flag <<- FALSE
        
        ergm_results[i] <<- list(NULL) 
        ergm_fits[i] <<- list(NULL) 
        t[i] <<- list(NULL) 
        mahalm[i] <<- list(NULL) 
        tideg[i] <<- list(NULL) 
        ideg[i] <<- list(NULL) 
        todeg[i] <<- list(NULL) 
        odeg[i] <<- list(NULL) 
        tdist[i] <<- list(NULL) 
        dist[i] <<- list(NULL) 
        tesp[i] <<- list(NULL) 
        esp[i] <<- list(NULL) 
        tdsp[i] <<- list(NULL) 
        dsp[i] <<- list(NULL) 
        ttriad[i] <<- list(NULL) 
        triad[i] <<- list(NULL) 
        diffs[i] <<- list(NULL) 
        aic[i] <<- NA
        bic[i] <<- NA  
        
      }
      
    })
    
    if (!flag) next
    
  }
  
  print(ergm_results)
  print(diffs)
  ergm_results <- ergm_results[!sapply(ergm_results, is.null)]
  ergm_fits <- ergm_fits[!sapply(ergm_fits, is.null)]
  t <- t[!sapply(t, is.null)]
  mahalm <- mahalm[!sapply(mahalm, is.null)]
  tideg <- tideg[!sapply(tideg, is.null)]
  ideg <- ideg[!sapply(ideg, is.null)]
  todeg <- todeg[!sapply(todeg, is.null)]
  odeg <- odeg[!sapply(odeg, is.null)]
  tdist <- tdist[!sapply(tdist, is.null)]
  dist <- dist[!sapply(dist, is.null)]
  tesp <- tesp[!sapply(tesp, is.null)]
  esp <- esp[!sapply(esp, is.null)]
  tdsp <- tdsp[!sapply(tdsp, is.null)]
  dsp <- dsp[!sapply(dsp, is.null)]
  ttriad <- ttriad[!sapply(ttriad, is.null)]
  triad <- triad[!sapply(triad, is.null)]
  diffs <- diffs[!sapply(diffs, is.null)]
  aic <- aic[!is.na(aic)]
  bic <- bic[!is.na(bic)]
  print(ergm_results)
  print(diffs)
  print(length(diffs))
  
  ### Meta-analysis.
  ma_coeffs <- matrix(nrow = length(diffs), ncol = length(diffs[[1]]))
  ma_ses <- matrix(nrow = length(diffs), ncol = length(diffs[[1]]))
  ma_model <- list()
  for (i in 1 : length(diffs)){
    for (j in 1 : length(diffs[[1]])){
      ma_coeffs[i, ] <- ergm_results[[i]]$coef[1 : length(diffs[[1]])]
      ma_ses[i, ] <- sqrt(diag(ergm_results[[i]]$covar[1 : length(diffs[[1]]), 1 : length(diffs[[1]])]))
      suppressWarnings(ma_model[[j]] <- rma(yi = ma_coeffs[, j], sei = ma_ses[, j], control = list(stepadj = 0.5, maxiter = 1000)))
    }
  }
  
  ### Bias and rmse. 
  ## Bias.
  bias_matrix <- matrix(data = NA, nrow = length(diffs), ncol = length(diffs[[1]]))
  for (i in 1 : length(diffs)){
    bias_matrix[i, ] <- diffs[[i]]
  }
  bias <- colSums(bias_matrix) / length(diffs)
  
  ## rmse.
  rmse <- sqrt(colSums(bias_matrix ** 2) / length(diffs))
  
  ### Mean AIC and BIC.
  mean_aic <- mean(unlist(aic)) # Mean "meta-analysis" AIC.
  mean_bic <- mean(unlist(bic)) # Mean "meta-analysis" BIC.
  
  return(list(next_counter = next_counter,
              model_GOF = list(t = t, mahalm = mahalm), indegree_aux = list(tideg = tideg, ideg = ideg), 
              outdegree_aux = list(todeg = todeg, odeg = odeg), dist_aux = list(tdist = tdist, dist = dist), 
              esp_aux = list(tesp = tesp, esp = esp), dsp_aux = list(tdsp = tdsp, dsp = dsp), 
              triad_aux = list(ttriad = ttriad, triad = triad), aic = list(aic = unlist(aic), mean_aic = mean_aic), 
              bic = list(bic = unlist(bic), mean_bic = mean_bic), ma_model = ma_model,  
              bias = bias, rmse = rmse, seed1 = seed1, seed2 = seed2))
  
}

cov_sim_het1 <- function(nsim){
  
  ### Define objects / lists for storing function output. 
  ## Next counter.
  next_counter <- 0
  ## ERGM results. 
  ergm_results <- vector(mode = "list", length = nsim)
  ## ERGM fits. 
  ergm_fits <- vector(mode = "list", length = nsim)
  ## List for storing model statistics. 
  t <- vector(mode = "list", length = nsim)
  mahalm <- vector(mode = "list", length = nsim)
  ## List for storing in-degree GOF. 
  tideg <- vector(mode = "list", length = nsim)
  ideg <- vector(mode = "list", length = nsim)
  ## List for storing out-degree GOF. 
  todeg <- vector(mode = "list", length = nsim)
  odeg <- vector(mode = "list", length = nsim)
  ## List for storing geodesic distance GOF. 
  tdist <- vector(mode = "list", length = nsim)
  dist <- vector(mode = "list", length = nsim)
  ## List for storing edgewise shared partners GOF. 
  tesp <- vector(mode = "list", length = nsim)
  esp <- vector(mode = "list", length = nsim)
  ## List for storing dyadwise shared partners GOF. 
  tdsp <- vector(mode = "list", length = nsim)
  dsp <- vector(mode = "list", length = nsim)
  ## List for storing triadic census GOF. 
  ttriad <- vector(mode = "list", length = nsim)
  triad <- vector(mode = "list", length = nsim)
  ## List for model parameter estimate differences relative to true model. 
  diffs <- vector(mode = "list", length = nsim)
  ## List for storing AIC and BIC.
  aic <- vector(mode = "list", length = nsim)
  bic <- vector(mode = "list", length = nsim)
  ## Seeds for reproducibility.
  seed1 <- as.integer(Sys.time()) %% 1000
  seed2 <- as.integer(Sys.time()) %% 10000 
  
  for (i in 1 : (nsim * 0.80)){
    
    cat(paste0("\nCondition: param 'cov' @ size ", nsim, " @ het_1. Iteration ", i, ".\n")) 
    
    tryCatch({
      
      ### Define model.
      model <- ergm_sims_het1[[i]] ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex")
      ### Fit model to "i"th simulated network.
      ergm_results[[i]] <- ergm(formula = model, control = control.ergm(MCMC.interval = 1 * 1024, MCMC.burnin = (1 * 1024) * 16, MCMC.samplesize = 1 * 1024, 
                                                                        seed = seed1))
      ### Calculate model and auxillary statistics fit of model. 
      ergm_fits[[i]] <- gof(object = model, GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus, 
                            coef = ergm_results[[i]]$coef, control.gof.formula(seed = seed2, nsim = 1000))
      
      ### GOF. 
      
      ## Model.
      # t-statistic.
      meant <- apply(ergm_fits[[i]]$sim.model, 2, mean)
      sdt <- apply(ergm_fits[[i]]$sim.model, 2, sd)
      obst <- summary(model)
      t[[i]] <- (meant-obst) / sdt
      # Mahalanobis.
      range <- which(sdt != 0)
      covm <- cov(as.matrix(ergm_fits[[i]]$sim.model[, range]))
      rangem <- range[which(eigen(covm)$values > 0.001)]
      covmodel <- cov(as.matrix(ergm_fits[[i]]$sim.model[, rangem]))
      mahalmodel <- mahalanobis(meant[rangem], ergm_fits[[i]]$obs.model[rangem], covmodel, tol = 1e-20)
      dfm <- length(rangem)
      pm <- 1 - pchisq(mahalmodel, dfm)
      mahalm[[i]] <- c(mahalmodel, dfm, pm)
      
      ## In-degree.
      # t-statistic.
      meantideg <- apply(ergm_fits[[i]]$sim.ideg, 2, mean)
      sdtideg <- apply(ergm_fits[[i]]$sim.ideg, 2, sd)
      tideg[[i]] <- (meantideg - ergm_fits[[i]]$obs.ideg) / sdtideg
      # Mahalanobis.
      range <- which(sdtideg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[,range])) 
      rangei <- range[which(eigen(covc)$values > 0.001)] 
      covideg <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[, rangei]))
      mahalideg <- mahalanobis(meantideg[rangei], ergm_fits[[i]]$obs.ideg[rangei], covideg, tol = 1e-20)
      dfi <- length(rangei)
      pi <- 1 - pchisq(mahalideg, dfi)
      ideg[[i]] <- c(mahalideg, dfi, pi)
      
      ## Out-degree.
      # t-statistic.
      meantodeg <- apply(ergm_fits[[i]]$sim.odeg,2,mean)
      sdtodeg <- apply(ergm_fits[[i]]$sim.odeg,2,sd)
      todeg[[i]] <- (meantodeg-ergm_fits[[i]]$obs.odeg)/sdtodeg
      # Mahalanobis.
      range <- which(sdtodeg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, range]))
      rangeo <- range[which(eigen(covc)$values > 0.001)] 
      covodeg <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, rangeo]))
      mahalodeg <- mahalanobis(meantodeg[rangeo], ergm_fits[[i]]$obs.odeg[rangeo], covodeg, tol = 1e-20)
      dfo <- length(rangeo)
      po <- 1 - pchisq(mahalodeg, dfo)
      odeg[[i]] <- c(mahalodeg, dfo, po)
      
      ## Geodesic distance.
      # t-statistic.
      meantdist <- apply(ergm_fits[[i]]$sim.dist, 2, mean)
      sdtdist <- apply(ergm_fits[[i]]$sim.dist, 2, sd)
      tdist[[i]] <- (meantdist - ergm_fits[[i]]$obs.dist) / sdtdist
      # Mahalanobis.
      range <- which(sdtdist != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, range]))
      ranged <- range[which(eigen(covc)$values > 0.001)] 
      covdist <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, ranged]))
      mahaldist <- mahalanobis(meantdist[ranged], ergm_fits[[i]]$obs.dist[ranged], covdist, tol = 1e-20)
      dfdist <- length(ranged)
      pdist <- 1 - pchisq(mahaldist, dfdist)
      dist[[i]] <- c(mahaldist,dfdist,pdist)
      
      ## Edgewise shared partners.
      # t-statistic.
      meantesp <- apply(ergm_fits[[i]]$sim.esp, 2, mean)
      sdtesp <- apply(ergm_fits[[i]]$sim.esp, 2, sd)
      tesp[[i]] <- (meantesp - ergm_fits[[i]]$obs.esp) / sdtesp
      # Mahalanobis.
      range <- which(sdtesp != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, range]))
      rangee <- range[which(eigen(covc)$values > 0.001)] 
      covesp <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, rangee]))
      mahalesp <- mahalanobis(meantesp[rangee], ergm_fits[[i]]$obs.esp[rangee], covesp, tol = 1e-20)
      dfesp <- length(rangee)
      pesp <- 1 - pchisq(mahalesp, dfesp)
      esp[[i]] <- c(mahalesp, dfesp, pesp)
      
      ## Dyadwise shared partners.
      # t-statistic.
      meantdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, mean)
      sdtdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, sd)
      tdsp[[i]] <- c(meantdsp - ergm_fits[[i]]$obs.dsp) / sdtdsp
      # Mahalanobis.
      range <- which(sdtdsp != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, range]))
      rangeds <- range[which(eigen(covc)$values > 0.001)] 
      covdsp <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, rangeds]))
      mahaldsp <- mahalanobis(meantdsp[rangeds], ergm_fits[[i]]$obs.dsp[rangeds], covdsp, tol = 1e-20)
      dfdsp <- length(rangeds)
      pdsp <- 1 - pchisq(mahaldsp, dfdsp)
      dsp[[i]] <- c(mahaldsp, dfdsp, pdsp)
      
      ## Triadcensus.
      # t-statistic.
      meanttriad <- apply(ergm_fits[[i]]$sim.triad, 2, mean)
      sdttriad <- apply(ergm_fits[[i]]$sim.triad, 2, sd)
      ttriad[[i]] <- (meanttriad - ergm_fits[[i]]$obs.triad) / sdttriad
      # Mahalanobis.
      range <- which(sdttriad != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, range]))
      ranget <- range[which(eigen(covc)$values > 0.001)] 
      covtriad <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, ranget]))
      mahaltriad <- mahalanobis(meanttriad[ranget], ergm_fits[[i]]$obs.triad[ranget], covtriad, tol = 1e-20)
      dftriad <- length(ranget)
      ptriad <- 1 - pchisq(mahaltriad, dftriad)
      triad[[i]] <- c(mahaltriad, dftriad, ptriad)
      
      ## Bias and rmse.
      # Differences.
      if (i >= ((nsim * 0.80) + 1)){
        diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het1$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")]     
      } else{
        diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het0$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")]     
      }
      ## AIC / BIC.
      # AIC.
      aic[i] <- summary(ergm_results[[i]])[17]  
      # BIC.
      bic[i] <- summary(ergm_results[[i]])[18]
      
      retry <<- FALSE
      
    }, error = function(e){
      
      cat(paste0("Condition: param 'cov' @ size ", nsim, " @ het_1. Error in iteration ", i, ".1. Re-simulating and re-fitting.\n"))
      next_counter <<- next_counter + 1
      retry <<- TRUE
      
    })
    
    if(retry){
      
      k <- 1
      err <<- TRUE
      
      while(isTRUE(err)){
        
        tryCatch({
          
          k <- k + 1
          
          ### Resimulate.
          ergm_sims_het1[[i]] <<- simulate(true_ergm_results_het1, nsim = 1, set.seed = seed1)
          ### Define model.
          model <- ergm_sims_het1[[i]] ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex")
          ### Fit model to "i"th simulated network.
          ergm_results[[i]] <- ergm(formula = model, control = control.ergm(MCMC.interval = 1 * 1024, MCMC.burnin = (1 * 1024) * 16, MCMC.samplesize = 1 * 1024, 
                                                                            seed = seed1))
          ### Calculate model and auxillary statistics fit of model. 
          ergm_fits[[i]] <- gof(object = model, GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus, 
                                coef = ergm_results[[i]]$coef, control.gof.formula(seed = seed2, nsim = 1000))
          
          ### GOF. 
          
          ## Model.
          # t-statistic.
          meant <- apply(ergm_fits[[i]]$sim.model, 2, mean)
          sdt <- apply(ergm_fits[[i]]$sim.model, 2, sd)
          obst <- summary(model)
          t[[i]] <- (meant-obst) / sdt
          # Mahalanobis.
          range <- which(sdt != 0)
          covm <- cov(as.matrix(ergm_fits[[i]]$sim.model[, range]))
          rangem <- range[which(eigen(covm)$values > 0.001)]
          covmodel <- cov(as.matrix(ergm_fits[[i]]$sim.model[, rangem]))
          mahalmodel <- mahalanobis(meant[rangem], ergm_fits[[i]]$obs.model[rangem], covmodel, tol = 1e-20)
          dfm <- length(rangem)
          pm <- 1 - pchisq(mahalmodel, dfm)
          mahalm[[i]] <- c(mahalmodel, dfm, pm)
          
          ## In-degree.
          # t-statistic.
          meantideg <- apply(ergm_fits[[i]]$sim.ideg, 2, mean)
          sdtideg <- apply(ergm_fits[[i]]$sim.ideg, 2, sd)
          tideg[[i]] <- (meantideg - ergm_fits[[i]]$obs.ideg) / sdtideg
          # Mahalanobis.
          range <- which(sdtideg != 0)
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[,range])) 
          rangei <- range[which(eigen(covc)$values > 0.001)] 
          covideg <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[, rangei]))
          mahalideg <- mahalanobis(meantideg[rangei], ergm_fits[[i]]$obs.ideg[rangei], covideg, tol = 1e-20)
          dfi <- length(rangei)
          pi <- 1 - pchisq(mahalideg, dfi)
          ideg[[i]] <- c(mahalideg, dfi, pi)
          
          ## Out-degree.
          # t-statistic.
          meantodeg <- apply(ergm_fits[[i]]$sim.odeg,2,mean)
          sdtodeg <- apply(ergm_fits[[i]]$sim.odeg,2,sd)
          todeg[[i]] <- (meantodeg-ergm_fits[[i]]$obs.odeg)/sdtodeg
          # Mahalanobis.
          range <- which(sdtodeg != 0)
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, range]))
          rangeo <- range[which(eigen(covc)$values > 0.001)] 
          covodeg <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, rangeo]))
          mahalodeg <- mahalanobis(meantodeg[rangeo], ergm_fits[[i]]$obs.odeg[rangeo], covodeg, tol = 1e-20)
          dfo <- length(rangeo)
          po <- 1 - pchisq(mahalodeg, dfo)
          odeg[[i]] <- c(mahalodeg, dfo, po)
          
          ## Geodesic distance.
          # t-statistic.
          meantdist <- apply(ergm_fits[[i]]$sim.dist, 2, mean)
          sdtdist <- apply(ergm_fits[[i]]$sim.dist, 2, sd)
          tdist[[i]] <- (meantdist - ergm_fits[[i]]$obs.dist) / sdtdist
          # Mahalanobis.
          range <- which(sdtdist != 0)
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, range]))
          ranged <- range[which(eigen(covc)$values > 0.001)] 
          covdist <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, ranged]))
          mahaldist <- mahalanobis(meantdist[ranged], ergm_fits[[i]]$obs.dist[ranged], covdist, tol = 1e-20)
          dfdist <- length(ranged)
          pdist <- 1 - pchisq(mahaldist, dfdist)
          dist[[i]] <- c(mahaldist,dfdist,pdist)
          
          ## Edgewise shared partners.
          # t-statistic.
          meantesp <- apply(ergm_fits[[i]]$sim.esp, 2, mean)
          sdtesp <- apply(ergm_fits[[i]]$sim.esp, 2, sd)
          tesp[[i]] <- (meantesp - ergm_fits[[i]]$obs.esp) / sdtesp
          # Mahalanobis.
          range <- which(sdtesp != 0) 
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, range]))
          rangee <- range[which(eigen(covc)$values > 0.001)] 
          covesp <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, rangee]))
          mahalesp <- mahalanobis(meantesp[rangee], ergm_fits[[i]]$obs.esp[rangee], covesp, tol = 1e-20)
          dfesp <- length(rangee)
          pesp <- 1 - pchisq(mahalesp, dfesp)
          esp[[i]] <- c(mahalesp, dfesp, pesp)
          
          ## Dyadwise shared partners.
          # t-statistic.
          meantdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, mean)
          sdtdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, sd)
          tdsp[[i]] <- c(meantdsp - ergm_fits[[i]]$obs.dsp) / sdtdsp
          # Mahalanobis.
          range <- which(sdtdsp != 0)
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, range]))
          rangeds <- range[which(eigen(covc)$values > 0.001)] 
          covdsp <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, rangeds]))
          mahaldsp <- mahalanobis(meantdsp[rangeds], ergm_fits[[i]]$obs.dsp[rangeds], covdsp, tol = 1e-20)
          dfdsp <- length(rangeds)
          pdsp <- 1 - pchisq(mahaldsp, dfdsp)
          dsp[[i]] <- c(mahaldsp, dfdsp, pdsp)
          
          ## Triadcensus.
          # t-statistic.
          meanttriad <- apply(ergm_fits[[i]]$sim.triad, 2, mean)
          sdttriad <- apply(ergm_fits[[i]]$sim.triad, 2, sd)
          ttriad[[i]] <- (meanttriad - ergm_fits[[i]]$obs.triad) / sdttriad
          # Mahalanobis.
          range <- which(sdttriad != 0) 
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, range]))
          ranget <- range[which(eigen(covc)$values > 0.001)] 
          covtriad <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, ranget]))
          mahaltriad <- mahalanobis(meanttriad[ranget], ergm_fits[[i]]$obs.triad[ranget], covtriad, tol = 1e-20)
          dftriad <- length(ranget)
          ptriad <- 1 - pchisq(mahaltriad, dftriad)
          triad[[i]] <- c(mahaltriad, dftriad, ptriad)
          
          ## Bias and rmse.
          # Differences.
          if (i >= ((nsim * 0.80) + 1)){
            diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het1$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")]     
          } else{
            diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het0$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")]     
          }
          ## AIC / BIC.
          # AIC.
          aic[i] <- summary(ergm_results[[i]])[17]  
          # BIC.
          bic[i] <- summary(ergm_results[[i]])[18]
          
          err <<- FALSE
          
        }, error = function(e){
          
          cat(paste0("Condition: param 'cov' @ size ", nsim, " @ het_1. Error in iteration ", i, ".", k, ". Re-simulating and re-fitting.\n"))
          err <<- TRUE
          
        })
        
      }
      
    }
  }
  
  for (i in ((nsim * 0.80) + 1) : nsim){
    
    cat(paste0("\nCondition: param 'cov' @ size ", nsim, " @ het_1. Iteration ", i, ".\n")) 
    
    tryCatch({
      
      ### Define model.
      model <- ergm_sims_het1[[i]] ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex")
      ### Fit model to "i"th simulated network.
      ergm_results[[i]] <- ergm(formula = model, control = control.ergm(MCMC.interval = 1 * 1024, MCMC.burnin = (1 * 1024) * 16, MCMC.samplesize = 1 * 1024, 
                                                                        seed = seed1))
      ### Calculate model and auxillary statistics fit of model. 
      ergm_fits[[i]] <- gof(object = model, GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus, 
                            coef = ergm_results[[i]]$coef, control.gof.formula(seed = seed2, nsim = 1000))
      
      ### GOF. 
      
      ## Model.
      # t-statistic.
      meant <- apply(ergm_fits[[i]]$sim.model, 2, mean)
      sdt <- apply(ergm_fits[[i]]$sim.model, 2, sd)
      obst <- summary(model)
      t[[i]] <- (meant-obst) / sdt
      # Mahalanobis.
      range <- which(sdt != 0)
      covm <- cov(as.matrix(ergm_fits[[i]]$sim.model[, range]))
      rangem <- range[which(eigen(covm)$values > 0.001)]
      covmodel <- cov(as.matrix(ergm_fits[[i]]$sim.model[, rangem]))
      mahalmodel <- mahalanobis(meant[rangem], ergm_fits[[i]]$obs.model[rangem], covmodel, tol = 1e-20)
      dfm <- length(rangem)
      pm <- 1 - pchisq(mahalmodel, dfm)
      mahalm[[i]] <- c(mahalmodel, dfm, pm)
      
      ## In-degree.
      # t-statistic.
      meantideg <- apply(ergm_fits[[i]]$sim.ideg, 2, mean)
      sdtideg <- apply(ergm_fits[[i]]$sim.ideg, 2, sd)
      tideg[[i]] <- (meantideg - ergm_fits[[i]]$obs.ideg) / sdtideg
      # Mahalanobis.
      range <- which(sdtideg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[,range])) 
      rangei <- range[which(eigen(covc)$values > 0.001)] 
      covideg <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[, rangei]))
      mahalideg <- mahalanobis(meantideg[rangei], ergm_fits[[i]]$obs.ideg[rangei], covideg, tol = 1e-20)
      dfi <- length(rangei)
      pi <- 1 - pchisq(mahalideg, dfi)
      ideg[[i]] <- c(mahalideg, dfi, pi)
      
      ## Out-degree.
      # t-statistic.
      meantodeg <- apply(ergm_fits[[i]]$sim.odeg,2,mean)
      sdtodeg <- apply(ergm_fits[[i]]$sim.odeg,2,sd)
      todeg[[i]] <- (meantodeg-ergm_fits[[i]]$obs.odeg)/sdtodeg
      # Mahalanobis.
      range <- which(sdtodeg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, range]))
      rangeo <- range[which(eigen(covc)$values > 0.001)] 
      covodeg <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, rangeo]))
      mahalodeg <- mahalanobis(meantodeg[rangeo], ergm_fits[[i]]$obs.odeg[rangeo], covodeg, tol = 1e-20)
      dfo <- length(rangeo)
      po <- 1 - pchisq(mahalodeg, dfo)
      odeg[[i]] <- c(mahalodeg, dfo, po)
      
      ## Geodesic distance.
      # t-statistic.
      meantdist <- apply(ergm_fits[[i]]$sim.dist, 2, mean)
      sdtdist <- apply(ergm_fits[[i]]$sim.dist, 2, sd)
      tdist[[i]] <- (meantdist - ergm_fits[[i]]$obs.dist) / sdtdist
      # Mahalanobis.
      range <- which(sdtdist != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, range]))
      ranged <- range[which(eigen(covc)$values > 0.001)] 
      covdist <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, ranged]))
      mahaldist <- mahalanobis(meantdist[ranged], ergm_fits[[i]]$obs.dist[ranged], covdist, tol = 1e-20)
      dfdist <- length(ranged)
      pdist <- 1 - pchisq(mahaldist, dfdist)
      dist[[i]] <- c(mahaldist,dfdist,pdist)
      
      ## Edgewise shared partners.
      # t-statistic.
      meantesp <- apply(ergm_fits[[i]]$sim.esp, 2, mean)
      sdtesp <- apply(ergm_fits[[i]]$sim.esp, 2, sd)
      tesp[[i]] <- (meantesp - ergm_fits[[i]]$obs.esp) / sdtesp
      # Mahalanobis.
      range <- which(sdtesp != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, range]))
      rangee <- range[which(eigen(covc)$values > 0.001)] 
      covesp <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, rangee]))
      mahalesp <- mahalanobis(meantesp[rangee], ergm_fits[[i]]$obs.esp[rangee], covesp, tol = 1e-20)
      dfesp <- length(rangee)
      pesp <- 1 - pchisq(mahalesp, dfesp)
      esp[[i]] <- c(mahalesp, dfesp, pesp)
      
      ## Dyadwise shared partners.
      # t-statistic.
      meantdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, mean)
      sdtdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, sd)
      tdsp[[i]] <- c(meantdsp - ergm_fits[[i]]$obs.dsp) / sdtdsp
      # Mahalanobis.
      range <- which(sdtdsp != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, range]))
      rangeds <- range[which(eigen(covc)$values > 0.001)] 
      covdsp <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, rangeds]))
      mahaldsp <- mahalanobis(meantdsp[rangeds], ergm_fits[[i]]$obs.dsp[rangeds], covdsp, tol = 1e-20)
      dfdsp <- length(rangeds)
      pdsp <- 1 - pchisq(mahaldsp, dfdsp)
      dsp[[i]] <- c(mahaldsp, dfdsp, pdsp)
      
      ## Triadcensus.
      # t-statistic.
      meanttriad <- apply(ergm_fits[[i]]$sim.triad, 2, mean)
      sdttriad <- apply(ergm_fits[[i]]$sim.triad, 2, sd)
      ttriad[[i]] <- (meanttriad - ergm_fits[[i]]$obs.triad) / sdttriad
      # Mahalanobis.
      range <- which(sdttriad != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, range]))
      ranget <- range[which(eigen(covc)$values > 0.001)] 
      covtriad <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, ranget]))
      mahaltriad <- mahalanobis(meanttriad[ranget], ergm_fits[[i]]$obs.triad[ranget], covtriad, tol = 1e-20)
      dftriad <- length(ranget)
      ptriad <- 1 - pchisq(mahaltriad, dftriad)
      triad[[i]] <- c(mahaltriad, dftriad, ptriad)
      
      ## Bias and rmse.
      # Differences.
      if (i >= ((nsim * 0.80) + 1)){
        diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het1$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")]     
      } else{
        diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het0$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")]     
      }
      ## AIC / BIC.
      # AIC.
      aic[i] <- summary(ergm_results[[i]])[17]  
      # BIC.
      bic[i] <- summary(ergm_results[[i]])[18]
      
      retry <<- FALSE
      
    }, error = function(e){
      
      cat(paste0("Condition: param 'cov' @ size ", nsim, " @ het_1. Error in iteration ", i, ".1. Re-simulating and re-fitting.\n"))
      next_counter <<- next_counter + 1
      retry <<- TRUE
      
    })
    
    if(retry){
      
      k <- 1
      err <<- TRUE
      
      while(isTRUE(err)){
        
        tryCatch({
          
          k <- k + 1
          
          ### Resimulate.
          ergm_sims_het1[[i]] <<- simulate(true_ergm_results_het1, nsim = 1, set.seed = seed1)
          ### Define model.
          model <- ergm_sims_het1[[i]] ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex")
          ### Fit model to "i"th simulated network.
          ergm_results[[i]] <- ergm(formula = model, control = control.ergm(MCMC.interval = 1 * 1024, MCMC.burnin = (1 * 1024) * 16, MCMC.samplesize = 1 * 1024, 
                                                                            seed = seed1))
          ### Calculate model and auxillary statistics fit of model. 
          ergm_fits[[i]] <- gof(object = model, GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus, 
                                coef = ergm_results[[i]]$coef, control.gof.formula(seed = seed2, nsim = 1000))
          
          ### GOF. 
          
          ## Model.
          # t-statistic.
          meant <- apply(ergm_fits[[i]]$sim.model, 2, mean)
          sdt <- apply(ergm_fits[[i]]$sim.model, 2, sd)
          obst <- summary(model)
          t[[i]] <- (meant-obst) / sdt
          # Mahalanobis.
          range <- which(sdt != 0)
          covm <- cov(as.matrix(ergm_fits[[i]]$sim.model[, range]))
          rangem <- range[which(eigen(covm)$values > 0.001)]
          covmodel <- cov(as.matrix(ergm_fits[[i]]$sim.model[, rangem]))
          mahalmodel <- mahalanobis(meant[rangem], ergm_fits[[i]]$obs.model[rangem], covmodel, tol = 1e-20)
          dfm <- length(rangem)
          pm <- 1 - pchisq(mahalmodel, dfm)
          mahalm[[i]] <- c(mahalmodel, dfm, pm)
          
          ## In-degree.
          # t-statistic.
          meantideg <- apply(ergm_fits[[i]]$sim.ideg, 2, mean)
          sdtideg <- apply(ergm_fits[[i]]$sim.ideg, 2, sd)
          tideg[[i]] <- (meantideg - ergm_fits[[i]]$obs.ideg) / sdtideg
          # Mahalanobis.
          range <- which(sdtideg != 0)
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[,range])) 
          rangei <- range[which(eigen(covc)$values > 0.001)] 
          covideg <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[, rangei]))
          mahalideg <- mahalanobis(meantideg[rangei], ergm_fits[[i]]$obs.ideg[rangei], covideg, tol = 1e-20)
          dfi <- length(rangei)
          pi <- 1 - pchisq(mahalideg, dfi)
          ideg[[i]] <- c(mahalideg, dfi, pi)
          
          ## Out-degree.
          # t-statistic.
          meantodeg <- apply(ergm_fits[[i]]$sim.odeg,2,mean)
          sdtodeg <- apply(ergm_fits[[i]]$sim.odeg,2,sd)
          todeg[[i]] <- (meantodeg-ergm_fits[[i]]$obs.odeg)/sdtodeg
          # Mahalanobis.
          range <- which(sdtodeg != 0)
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, range]))
          rangeo <- range[which(eigen(covc)$values > 0.001)] 
          covodeg <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, rangeo]))
          mahalodeg <- mahalanobis(meantodeg[rangeo], ergm_fits[[i]]$obs.odeg[rangeo], covodeg, tol = 1e-20)
          dfo <- length(rangeo)
          po <- 1 - pchisq(mahalodeg, dfo)
          odeg[[i]] <- c(mahalodeg, dfo, po)
          
          ## Geodesic distance.
          # t-statistic.
          meantdist <- apply(ergm_fits[[i]]$sim.dist, 2, mean)
          sdtdist <- apply(ergm_fits[[i]]$sim.dist, 2, sd)
          tdist[[i]] <- (meantdist - ergm_fits[[i]]$obs.dist) / sdtdist
          # Mahalanobis.
          range <- which(sdtdist != 0)
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, range]))
          ranged <- range[which(eigen(covc)$values > 0.001)] 
          covdist <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, ranged]))
          mahaldist <- mahalanobis(meantdist[ranged], ergm_fits[[i]]$obs.dist[ranged], covdist, tol = 1e-20)
          dfdist <- length(ranged)
          pdist <- 1 - pchisq(mahaldist, dfdist)
          dist[[i]] <- c(mahaldist,dfdist,pdist)
          
          ## Edgewise shared partners.
          # t-statistic.
          meantesp <- apply(ergm_fits[[i]]$sim.esp, 2, mean)
          sdtesp <- apply(ergm_fits[[i]]$sim.esp, 2, sd)
          tesp[[i]] <- (meantesp - ergm_fits[[i]]$obs.esp) / sdtesp
          # Mahalanobis.
          range <- which(sdtesp != 0) 
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, range]))
          rangee <- range[which(eigen(covc)$values > 0.001)] 
          covesp <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, rangee]))
          mahalesp <- mahalanobis(meantesp[rangee], ergm_fits[[i]]$obs.esp[rangee], covesp, tol = 1e-20)
          dfesp <- length(rangee)
          pesp <- 1 - pchisq(mahalesp, dfesp)
          esp[[i]] <- c(mahalesp, dfesp, pesp)
          
          ## Dyadwise shared partners.
          # t-statistic.
          meantdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, mean)
          sdtdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, sd)
          tdsp[[i]] <- c(meantdsp - ergm_fits[[i]]$obs.dsp) / sdtdsp
          # Mahalanobis.
          range <- which(sdtdsp != 0)
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, range]))
          rangeds <- range[which(eigen(covc)$values > 0.001)] 
          covdsp <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, rangeds]))
          mahaldsp <- mahalanobis(meantdsp[rangeds], ergm_fits[[i]]$obs.dsp[rangeds], covdsp, tol = 1e-20)
          dfdsp <- length(rangeds)
          pdsp <- 1 - pchisq(mahaldsp, dfdsp)
          dsp[[i]] <- c(mahaldsp, dfdsp, pdsp)
          
          ## Triadcensus.
          # t-statistic.
          meanttriad <- apply(ergm_fits[[i]]$sim.triad, 2, mean)
          sdttriad <- apply(ergm_fits[[i]]$sim.triad, 2, sd)
          ttriad[[i]] <- (meanttriad - ergm_fits[[i]]$obs.triad) / sdttriad
          # Mahalanobis.
          range <- which(sdttriad != 0) 
          covc <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, range]))
          ranget <- range[which(eigen(covc)$values > 0.001)] 
          covtriad <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, ranget]))
          mahaltriad <- mahalanobis(meanttriad[ranget], ergm_fits[[i]]$obs.triad[ranget], covtriad, tol = 1e-20)
          dftriad <- length(ranget)
          ptriad <- 1 - pchisq(mahaltriad, dftriad)
          triad[[i]] <- c(mahaltriad, dftriad, ptriad)
          
          ## Bias and rmse.
          # Differences.
          if (i >= ((nsim * 0.80) + 1)){
            diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het1$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")]     
          } else{
            diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het0$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")]     
          }
          ## AIC / BIC.
          # AIC.
          aic[i] <- summary(ergm_results[[i]])[17]  
          # BIC.
          bic[i] <- summary(ergm_results[[i]])[18]
          
          err <<- FALSE 
          
        }, error = function(e){
          
          cat(paste0("Condition: param 'cov' @ size ", nsim, " @ het_1. Error in iteration ", i, ".", k, ". Re-simulating and re-fitting.\n"))
          err <<- TRUE
          
        })
        
      }
      
    }
  }
  
  ### Meta-analysis.
  ma_coeffs <- matrix(nrow = nsim, ncol = length(diffs[[1]]))
  ma_ses <- matrix(nrow = nsim, ncol = length(diffs[[1]]))
  ma_model <- list()
  for (i in 1 : nsim){
    for (j in 1 : length(diffs[[1]])){
      ma_coeffs[i, ] <- ergm_results[[i]]$coef[1 : length(diffs[[1]])]
      ma_ses[i, ] <- sqrt(diag(ergm_results[[i]]$covar))
      suppressWarnings(ma_model[[j]] <- rma(yi = ma_coeffs[, j], sei = ma_ses[, j], control = list(stepadj = 0.5, maxiter = 1000)))
    }
  }
  
  ### Bias and rmse. 
  ## Bias.
  bias_matrix <- matrix(data = NA, nrow = nsim, ncol = length(diffs[[1]]))
  for (i in 1 : nsim){
    bias_matrix[i, ] <- diffs[[i]]
  }
  bias <- colSums(bias_matrix) / nsim
  
  ## Rmse.
  rmse <- sqrt(colSums(bias_matrix ** 2) / nsim)
  
  ### Mean AIC and BIC.
  mean_aic <- mean(unlist(aic)) # Mean "meta-analysis" AIC.
  mean_bic <- mean(unlist(bic)) # Mean "meta-analysis" BIC.
  
  return(list(next_counter = next_counter,
              model_GOF = list(t = t, mahalm = mahalm), indegree_aux = list(tideg = tideg, ideg = ideg), 
              outdegree_aux = list(todeg = todeg, odeg = odeg), dist_aux = list(tdist = tdist, dist = dist), 
              esp_aux = list(tesp = tesp, esp = esp), dsp_aux = list(tdsp = tdsp, dsp = dsp), 
              triad_aux = list(ttriad = ttriad, triad = triad), aic = list(aic = unlist(aic), mean_aic = mean_aic), 
              bic = list(bic = unlist(bic), mean_bic = mean_bic), ma_model = ma_model, 
              bias = bias, rmse = rmse, seed1 = seed1, seed2 = seed2))
  
}

gwodegree_sim_het1 <- function(nsim){
  
  ### Define objects / lists for storing function output. 
  ## Flag.
  flag <- TRUE
  ## Next counter.
  next_counter <- 0
  ## ERGM results. 
  ergm_results <- vector(mode = "list", length = nsim)
  ## ERGM fits. 
  ergm_fits <- vector(mode = "list", length = nsim)
  ## List for storing model statistics. 
  t <- vector(mode = "list", length = nsim)
  mahalm <- vector(mode = "list", length = nsim)
  ## List for storing in-degree GOF. 
  tideg <- vector(mode = "list", length = nsim)
  ideg <- vector(mode = "list", length = nsim)
  ## List for storing out-degree GOF. 
  todeg <- vector(mode = "list", length = nsim)
  odeg <- vector(mode = "list", length = nsim)
  ## List for storing geodesic distance GOF. 
  tdist <- vector(mode = "list", length = nsim)
  dist <- vector(mode = "list", length = nsim)
  ## List for storing edgewise shared partners GOF. 
  tesp <- vector(mode = "list", length = nsim)
  esp <- vector(mode = "list", length = nsim)
  ## List for storing dyadwise shared partners GOF. 
  tdsp <- vector(mode = "list", length = nsim)
  dsp <- vector(mode = "list", length = nsim)
  ## List for storing triadic census GOF. 
  ttriad <- vector(mode = "list", length = nsim)
  triad <- vector(mode = "list", length = nsim)
  ## List for model parameter estimate differences relative to true model. 
  diffs <- vector(mode = "list", length = nsim)
  ## List for storing AIC and BIC.
  aic <- vector(mode = "list", length = nsim)
  bic <- vector(mode = "list", length = nsim)
  ## Seeds for reproducibility.
  seed1 <- as.integer(Sys.time()) %% 1000
  seed2 <- as.integer(Sys.time()) %% 10000 
  
  for (i in 1 : nsim){
    
    cat(paste0("\nCondition: param 'gwodeg' @ size ", nsim, " @ het_1. Iteration ", i, ".\n"))    
    
    tryCatch({
      
      ### Define model.
      model <- ergm_sims_het1[[i]] ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T)
      ### Fit model to "i"th simulated network.
      ergm_results[[i]] <- ergm(formula = model, control = control.ergm(MCMC.interval = 1 * 1024, MCMC.burnin = (1 * 1024) * 16, MCMC.samplesize = 1 * 1024, 
                                                                        seed = seed1))
      ### Calculate model and auxillary statistics fit of model. 
      ergm_fits[[i]] <- gof(object = model, GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus, 
                            coef = ergm_results[[i]]$coef, control.gof.formula(seed = seed2, nsim = 1000))
      
      ### GOF. 
      
      ## Model.
      # t-statistic.
      meant <- apply(ergm_fits[[i]]$sim.model, 2, mean)
      sdt <- apply(ergm_fits[[i]]$sim.model, 2, sd)
      obst <- summary(model)
      t[[i]] <- (meant-obst) / sdt
      # Mahalanobis.
      range <- which(sdt != 0)
      covm <- cov(as.matrix(ergm_fits[[i]]$sim.model[, range]))
      rangem <- range[which(eigen(covm)$values > 0.001)]
      covmodel <- cov(as.matrix(ergm_fits[[i]]$sim.model[, rangem]))
      mahalmodel <- mahalanobis(meant[rangem], ergm_fits[[i]]$obs.model[rangem], covmodel, tol = 1e-20)
      dfm <- length(rangem)
      pm <- 1 - pchisq(mahalmodel, dfm)
      mahalm[[i]] <- c(mahalmodel, dfm, pm)
      
      ## In-degree.
      # t-statistic.
      meantideg <- apply(ergm_fits[[i]]$sim.ideg, 2, mean)
      sdtideg <- apply(ergm_fits[[i]]$sim.ideg, 2, sd)
      tideg[[i]] <- (meantideg - ergm_fits[[i]]$obs.ideg) / sdtideg
      # Mahalanobis.
      range <- which(sdtideg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[,range])) 
      rangei <- range[which(eigen(covc)$values > 0.001)] 
      covideg <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[, rangei]))
      mahalideg <- mahalanobis(meantideg[rangei], ergm_fits[[i]]$obs.ideg[rangei], covideg, tol = 1e-20)
      dfi <- length(rangei)
      pi <- 1 - pchisq(mahalideg, dfi)
      ideg[[i]] <- c(mahalideg, dfi, pi)
      
      ## Out-degree.
      # t-statistic.
      meantodeg <- apply(ergm_fits[[i]]$sim.odeg,2,mean)
      sdtodeg <- apply(ergm_fits[[i]]$sim.odeg,2,sd)
      todeg[[i]] <- (meantodeg-ergm_fits[[i]]$obs.odeg)/sdtodeg
      # Mahalanobis.
      range <- which(sdtodeg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, range]))
      rangeo <- range[which(eigen(covc)$values > 0.001)] 
      covodeg <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, rangeo]))
      mahalodeg <- mahalanobis(meantodeg[rangeo], ergm_fits[[i]]$obs.odeg[rangeo], covodeg, tol = 1e-20)
      dfo <- length(rangeo)
      po <- 1 - pchisq(mahalodeg, dfo)
      odeg[[i]] <- c(mahalodeg, dfo, po)
      
      ## Geodesic distance.
      # t-statistic.
      meantdist <- apply(ergm_fits[[i]]$sim.dist, 2, mean)
      sdtdist <- apply(ergm_fits[[i]]$sim.dist, 2, sd)
      tdist[[i]] <- (meantdist - ergm_fits[[i]]$obs.dist) / sdtdist
      # Mahalanobis.
      range <- which(sdtdist != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, range]))
      ranged <- range[which(eigen(covc)$values > 0.001)] 
      covdist <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, ranged]))
      mahaldist <- mahalanobis(meantdist[ranged], ergm_fits[[i]]$obs.dist[ranged], covdist, tol = 1e-20)
      dfdist <- length(ranged)
      pdist <- 1 - pchisq(mahaldist, dfdist)
      dist[[i]] <- c(mahaldist,dfdist,pdist)
      
      ## Edgewise shared partners.
      # t-statistic.
      meantesp <- apply(ergm_fits[[i]]$sim.esp, 2, mean)
      sdtesp <- apply(ergm_fits[[i]]$sim.esp, 2, sd)
      tesp[[i]] <- (meantesp - ergm_fits[[i]]$obs.esp) / sdtesp
      # Mahalanobis.
      range <- which(sdtesp != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, range]))
      rangee <- range[which(eigen(covc)$values > 0.001)] 
      covesp <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, rangee]))
      mahalesp <- mahalanobis(meantesp[rangee], ergm_fits[[i]]$obs.esp[rangee], covesp, tol = 1e-20)
      dfesp <- length(rangee)
      pesp <- 1 - pchisq(mahalesp, dfesp)
      esp[[i]] <- c(mahalesp, dfesp, pesp)
      
      ## Dyadwise shared partners.
      # t-statistic.
      meantdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, mean)
      sdtdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, sd)
      tdsp[[i]] <- c(meantdsp - ergm_fits[[i]]$obs.dsp) / sdtdsp
      # Mahalanobis.
      range <- which(sdtdsp != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, range]))
      rangeds <- range[which(eigen(covc)$values > 0.001)] 
      covdsp <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, rangeds]))
      mahaldsp <- mahalanobis(meantdsp[rangeds], ergm_fits[[i]]$obs.dsp[rangeds], covdsp, tol = 1e-20)
      dfdsp <- length(rangeds)
      pdsp <- 1 - pchisq(mahaldsp, dfdsp)
      dsp[[i]] <- c(mahaldsp, dfdsp, pdsp)
      
      ## Triadcensus.
      # t-statistic.
      meanttriad <- apply(ergm_fits[[i]]$sim.triad, 2, mean)
      sdttriad <- apply(ergm_fits[[i]]$sim.triad, 2, sd)
      ttriad[[i]] <- (meanttriad - ergm_fits[[i]]$obs.triad) / sdttriad
      # Mahalanobis.
      range <- which(sdttriad != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, range]))
      ranget <- range[which(eigen(covc)$values > 0.001)] 
      covtriad <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, ranget]))
      mahaltriad <- mahalanobis(meanttriad[ranget], ergm_fits[[i]]$obs.triad[ranget], covtriad, tol = 1e-20)
      dftriad <- length(ranget)
      ptriad <- 1 - pchisq(mahaltriad, dftriad)
      triad[[i]] <- c(mahaltriad, dftriad, ptriad)
      
      ## Bias and rmse.
      if (i >= ((nsim * 0.80) + 1)){
        # Differences.
        diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1")] - true_ergm_results_het1$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1")]     
      } else{
        diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1")] - true_ergm_results_het0$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1")]     
      }
      ## AIC / BIC.
      # AIC.
      aic[i] <- summary(ergm_results[[i]])[17]  
      # BIC.
      bic[i] <- summary(ergm_results[[i]])[18]
      
    }, error = function(e){
      
      cat(paste0("Condition: param 'gwodeg' @ size ", nsim, " @ het_1. Error in iteration ", i, ". Skipping to next iteration.\n"))
      next_counter <<- next_counter + 1
      flag <<- FALSE
      
      ergm_results[i] <<- list(NULL) 
      ergm_fits[i] <<- list(NULL) 
      t[i] <<- list(NULL) 
      mahalm[i] <<- list(NULL) 
      tideg[i] <<- list(NULL) 
      ideg[i] <<- list(NULL) 
      todeg[i] <<- list(NULL) 
      odeg[i] <<- list(NULL) 
      tdist[i] <<- list(NULL) 
      dist[i] <<- list(NULL) 
      tesp[i] <<- list(NULL) 
      esp[i] <<- list(NULL) 
      tdsp[i] <<- list(NULL) 
      dsp[i] <<- list(NULL) 
      ttriad[i] <<- list(NULL) 
      triad[i] <<- list(NULL) 
      diffs[i] <<- list(NULL) 
      aic[i] <<- NA
      bic[i] <<- NA  
      
    })
    
    if (!flag) next
    
  }
  
  print(ergm_results)
  print(diffs)
  ergm_results <- ergm_results[!sapply(ergm_results, is.null)]
  ergm_fits <- ergm_fits[!sapply(ergm_fits, is.null)]
  t <- t[!sapply(t, is.null)]
  mahalm <- mahalm[!sapply(mahalm, is.null)]
  tideg <- tideg[!sapply(tideg, is.null)]
  ideg <- ideg[!sapply(ideg, is.null)]
  todeg <- todeg[!sapply(todeg, is.null)]
  odeg <- odeg[!sapply(odeg, is.null)]
  tdist <- tdist[!sapply(tdist, is.null)]
  dist <- dist[!sapply(dist, is.null)]
  tesp <- tesp[!sapply(tesp, is.null)]
  esp <- esp[!sapply(esp, is.null)]
  tdsp <- tdsp[!sapply(tdsp, is.null)]
  dsp <- dsp[!sapply(dsp, is.null)]
  ttriad <- ttriad[!sapply(ttriad, is.null)]
  triad <- triad[!sapply(triad, is.null)]
  diffs <- diffs[!sapply(diffs, is.null)]
  aic <- aic[!is.na(aic)]
  bic <- bic[!is.na(bic)]
  print(ergm_results)
  print(diffs)
  print(length(diffs))
  
  ### Meta-analysis.
  ma_coeffs <- matrix(nrow = length(diffs), ncol = length(diffs[[1]]))
  ma_ses <- matrix(nrow = length(diffs), ncol = length(diffs[[1]]))
  ma_model <- list()
  for (i in 1 : length(diffs)){
    for (j in 1 : length(diffs[[1]])){
      ma_coeffs[i, ] <- ergm_results[[i]]$coef[1 : length(diffs[[1]])]
      ma_ses[i, ] <- sqrt(diag(ergm_results[[i]]$covar[1 : length(diffs[[1]]), 1 : length(diffs[[1]])]))
      suppressWarnings(ma_model[[j]] <- rma(yi = ma_coeffs[, j], sei = ma_ses[, j], control = list(stepadj = 0.5, maxiter = 1000)))
    }
  }
  
  ### Bias and rmse. 
  ## Bias.
  bias_matrix <- matrix(data = NA, nrow = length(diffs), ncol = length(diffs[[1]]))
  for (i in 1 : length(diffs)){
    bias_matrix[i, ] <- diffs[[i]]
  }
  bias <- colSums(bias_matrix) / length(diffs)
  
  ## rmse.
  rmse <- sqrt(colSums(bias_matrix ** 2) / length(diffs))
  
  ### Mean AIC and BIC.
  mean_aic <- mean(unlist(aic)) # Mean "meta-analysis" AIC.
  mean_bic <- mean(unlist(bic)) # Mean "meta-analysis" BIC.
  
  return(list(next_counter = next_counter,
              model_GOF = list(t = t, mahalm = mahalm), indegree_aux = list(tideg = tideg, ideg = ideg), 
              outdegree_aux = list(todeg = todeg, odeg = odeg), dist_aux = list(tdist = tdist, dist = dist), 
              esp_aux = list(tesp = tesp, esp = esp), dsp_aux = list(tdsp = tdsp, dsp = dsp), 
              triad_aux = list(ttriad = ttriad, triad = triad), aic = list(aic = unlist(aic), mean_aic = mean_aic), 
              bic = list(bic = unlist(bic), mean_bic = mean_bic), ma_model = ma_model, 
              bias = bias, rmse = rmse, seed1 = seed1, seed2 = seed2))
  
}

dtriad_sim_het1 <- function(nsim){
  
  ### Define objects / lists for storing function output. 
  ## Flag.
  flag <- TRUE
  ## Next counter.
  next_counter <- 0
  ## ERGM results. 
  ergm_results <- vector(mode = "list", length = nsim)
  ## ERGM fits. 
  ergm_fits <- vector(mode = "list", length = nsim)
  ## List for storing model statistics. 
  t <- vector(mode = "list", length = nsim)
  mahalm <- vector(mode = "list", length = nsim)
  ## List for storing in-degree GOF. 
  tideg <- vector(mode = "list", length = nsim)
  ideg <- vector(mode = "list", length = nsim)
  ## List for storing out-degree GOF. 
  todeg <- vector(mode = "list", length = nsim)
  odeg <- vector(mode = "list", length = nsim)
  ## List for storing geodesic distance GOF. 
  tdist <- vector(mode = "list", length = nsim)
  dist <- vector(mode = "list", length = nsim)
  ## List for storing edgewise shared partners GOF. 
  tesp <- vector(mode = "list", length = nsim)
  esp <- vector(mode = "list", length = nsim)
  ## List for storing dyadwise shared partners GOF. 
  tdsp <- vector(mode = "list", length = nsim)
  dsp <- vector(mode = "list", length = nsim)
  ## List for storing triadic census GOF. 
  ttriad <- vector(mode = "list", length = nsim)
  triad <- vector(mode = "list", length = nsim)
  ## List for model parameter estimate differences relative to true model. 
  diffs <- vector(mode = "list", length = nsim)
  ## List for storing AIC and BIC.
  aic <- vector(mode = "list", length = nsim)
  bic <- vector(mode = "list", length = nsim)
  ## Seeds for reproducibility.
  seed1 <- as.integer(Sys.time()) %% 1000
  seed2 <- as.integer(Sys.time()) %% 10000 
  
  for (i in 1 : nsim){
    
    cat(paste0("\nCondition: param 'dtriad' @ size ", nsim, " @ het_1. Iteration ", i, ".\n"))   
    
    tryCatch({
      
      ### Define model.
      model <- ergm_sims_het1[[i]] ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex") + triadcensus(14)
      ### Fit model to "i"th simulated network.
      ergm_results[[i]] <- ergm(formula = model, control = control.ergm(MCMC.interval = 1 * 1024, MCMC.burnin = (1 * 1024) * 16, MCMC.samplesize = 1 * 1024, 
                                                                        seed = seed1))
      ### Calculate model and auxillary statistics fit of model. 
      ergm_fits[[i]] <- gof(object = model, GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus, 
                            coef = ergm_results[[i]]$coef, control.gof.formula(seed = seed2, nsim = 1000))
      
      if (any(is.infinite(ergm_results[[i]]$coef))){
        
        ergm_results[i] <- list(NULL) 
        ergm_fits[i] <- list(NULL) 
        stop()
        
      } 
      
      ### GOF. 
      
      ## Model.
      # t-statistic.
      meant <- apply(ergm_fits[[i]]$sim.model, 2, mean)
      sdt <- apply(ergm_fits[[i]]$sim.model, 2, sd)
      obst <- summary(model)
      t[[i]] <- (meant-obst) / sdt
      # Mahalanobis.
      range <- which(sdt != 0)
      covm <- cov(as.matrix(ergm_fits[[i]]$sim.model[, range]))
      rangem <- range[which(eigen(covm)$values > 0.001)]
      covmodel <- cov(as.matrix(ergm_fits[[i]]$sim.model[, rangem]))
      mahalmodel <- mahalanobis(meant[rangem], ergm_fits[[i]]$obs.model[rangem], covmodel, tol = 1e-20)
      dfm <- length(rangem)
      pm <- 1 - pchisq(mahalmodel, dfm)
      mahalm[[i]] <- c(mahalmodel, dfm, pm)
      
      ## In-degree.
      # t-statistic.
      meantideg <- apply(ergm_fits[[i]]$sim.ideg, 2, mean)
      sdtideg <- apply(ergm_fits[[i]]$sim.ideg, 2, sd)
      tideg[[i]] <- (meantideg - ergm_fits[[i]]$obs.ideg) / sdtideg
      # Mahalanobis.
      range <- which(sdtideg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[,range])) 
      rangei <- range[which(eigen(covc)$values > 0.001)] 
      covideg <- cov(as.matrix(ergm_fits[[i]]$sim.ideg[, rangei]))
      mahalideg <- mahalanobis(meantideg[rangei], ergm_fits[[i]]$obs.ideg[rangei], covideg, tol = 1e-20)
      dfi <- length(rangei)
      pi <- 1 - pchisq(mahalideg, dfi)
      ideg[[i]] <- c(mahalideg, dfi, pi)
      
      ## Out-degree.
      # t-statistic.
      meantodeg <- apply(ergm_fits[[i]]$sim.odeg,2,mean)
      sdtodeg <- apply(ergm_fits[[i]]$sim.odeg,2,sd)
      todeg[[i]] <- (meantodeg-ergm_fits[[i]]$obs.odeg)/sdtodeg
      # Mahalanobis.
      range <- which(sdtodeg != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, range]))
      rangeo <- range[which(eigen(covc)$values > 0.001)] 
      covodeg <- cov(as.matrix(ergm_fits[[i]]$sim.odeg[, rangeo]))
      mahalodeg <- mahalanobis(meantodeg[rangeo], ergm_fits[[i]]$obs.odeg[rangeo], covodeg, tol = 1e-20)
      dfo <- length(rangeo)
      po <- 1 - pchisq(mahalodeg, dfo)
      odeg[[i]] <- c(mahalodeg, dfo, po)
      
      ## Geodesic distance.
      # t-statistic.
      meantdist <- apply(ergm_fits[[i]]$sim.dist, 2, mean)
      sdtdist <- apply(ergm_fits[[i]]$sim.dist, 2, sd)
      tdist[[i]] <- (meantdist - ergm_fits[[i]]$obs.dist) / sdtdist
      # Mahalanobis.
      range <- which(sdtdist != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, range]))
      ranged <- range[which(eigen(covc)$values > 0.001)] 
      covdist <- cov(as.matrix(ergm_fits[[i]]$sim.dist[, ranged]))
      mahaldist <- mahalanobis(meantdist[ranged], ergm_fits[[i]]$obs.dist[ranged], covdist, tol = 1e-20)
      dfdist <- length(ranged)
      pdist <- 1 - pchisq(mahaldist, dfdist)
      dist[[i]] <- c(mahaldist,dfdist,pdist)
      
      ## Edgewise shared partners.
      # t-statistic.
      meantesp <- apply(ergm_fits[[i]]$sim.esp, 2, mean)
      sdtesp <- apply(ergm_fits[[i]]$sim.esp, 2, sd)
      tesp[[i]] <- (meantesp - ergm_fits[[i]]$obs.esp) / sdtesp
      # Mahalanobis.
      range <- which(sdtesp != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, range]))
      rangee <- range[which(eigen(covc)$values > 0.001)] 
      covesp <- cov(as.matrix(ergm_fits[[i]]$sim.esp[, rangee]))
      mahalesp <- mahalanobis(meantesp[rangee], ergm_fits[[i]]$obs.esp[rangee], covesp, tol = 1e-20)
      dfesp <- length(rangee)
      pesp <- 1 - pchisq(mahalesp, dfesp)
      esp[[i]] <- c(mahalesp, dfesp, pesp)
      
      ## Dyadwise shared partners.
      # t-statistic.
      meantdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, mean)
      sdtdsp <- apply(ergm_fits[[i]]$sim.dsp, 2, sd)
      tdsp[[i]] <- c(meantdsp - ergm_fits[[i]]$obs.dsp) / sdtdsp
      # Mahalanobis.
      range <- which(sdtdsp != 0)
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, range]))
      rangeds <- range[which(eigen(covc)$values > 0.001)] 
      covdsp <- cov(as.matrix(ergm_fits[[i]]$sim.dsp[, rangeds]))
      mahaldsp <- mahalanobis(meantdsp[rangeds], ergm_fits[[i]]$obs.dsp[rangeds], covdsp, tol = 1e-20)
      dfdsp <- length(rangeds)
      pdsp <- 1 - pchisq(mahaldsp, dfdsp)
      dsp[[i]] <- c(mahaldsp, dfdsp, pdsp)
      
      ## Triadcensus.
      # t-statistic.
      meanttriad <- apply(ergm_fits[[i]]$sim.triad, 2, mean)
      sdttriad <- apply(ergm_fits[[i]]$sim.triad, 2, sd)
      ttriad[[i]] <- (meanttriad - ergm_fits[[i]]$obs.triad) / sdttriad
      # Mahalanobis.
      range <- which(sdttriad != 0) 
      covc <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, range]))
      ranget <- range[which(eigen(covc)$values > 0.001)] 
      covtriad <- cov(as.matrix(ergm_fits[[i]]$sim.triad[, ranget]))
      mahaltriad <- mahalanobis(meanttriad[ranget], ergm_fits[[i]]$obs.triad[ranget], covtriad, tol = 1e-20)
      dftriad <- length(ranget)
      ptriad <- 1 - pchisq(mahaltriad, dftriad)
      triad[[i]] <- c(mahaltriad, dftriad, ptriad)
      
      ## Bias and rmse.
      if (i >= ((nsim * 0.80) + 1)){
        # Differences.
        diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het1$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")]     
      } else{
        diffs[[i]] <- ergm_results[[i]]$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")] - true_ergm_results_het0$coef[c("edges", "mutual", "gwesp.fixed.0.25", "gwodeg.fixed.0.1", "nodeocov.identclass", "nodematch.sex")]     
      }
      ## AIC / BIC.
      # AIC.
      aic[i] <- summary(ergm_results[[i]])[17]  
      # BIC.
      bic[i] <- summary(ergm_results[[i]])[18]
      
    }, error = function(e){
      
      if (e$message == ""){
        
        cat(paste0("Condition: param 'dtriad' @ size ", nsim, " @ het_1. Infinite coefficient in iteration ", i,". Skipping to next iteration.\n")) 
        next_counter <<- next_counter + 1
        flag <<- FALSE
        
      } else{
        
        cat(paste0("Condition: param 'dtriad' @ size ", nsim, " @ het_1. Error in iteration ", i, ". Skipping to next iteration.\n"))
        next_counter <<- next_counter + 1
        flag <<- FALSE
        
        ergm_results[i] <<- list(NULL) 
        ergm_fits[i] <<- list(NULL) 
        t[i] <<- list(NULL) 
        mahalm[i] <<- list(NULL) 
        tideg[i] <<- list(NULL) 
        ideg[i] <<- list(NULL) 
        todeg[i] <<- list(NULL) 
        odeg[i] <<- list(NULL) 
        tdist[i] <<- list(NULL) 
        dist[i] <<- list(NULL) 
        tesp[i] <<- list(NULL) 
        esp[i] <<- list(NULL) 
        tdsp[i] <<- list(NULL) 
        dsp[i] <<- list(NULL) 
        ttriad[i] <<- list(NULL) 
        triad[i] <<- list(NULL) 
        diffs[i] <<- list(NULL) 
        aic[i] <<- NA
        bic[i] <<- NA  
        
      }
      
    })
    
    if (!flag) next
    
  }
  
  print(ergm_results)
  print(diffs)
  ergm_results <- ergm_results[!sapply(ergm_results, is.null)]
  ergm_fits <- ergm_fits[!sapply(ergm_fits, is.null)]
  t <- t[!sapply(t, is.null)]
  mahalm <- mahalm[!sapply(mahalm, is.null)]
  tideg <- tideg[!sapply(tideg, is.null)]
  ideg <- ideg[!sapply(ideg, is.null)]
  todeg <- todeg[!sapply(todeg, is.null)]
  odeg <- odeg[!sapply(odeg, is.null)]
  tdist <- tdist[!sapply(tdist, is.null)]
  dist <- dist[!sapply(dist, is.null)]
  tesp <- tesp[!sapply(tesp, is.null)]
  esp <- esp[!sapply(esp, is.null)]
  tdsp <- tdsp[!sapply(tdsp, is.null)]
  dsp <- dsp[!sapply(dsp, is.null)]
  ttriad <- ttriad[!sapply(ttriad, is.null)]
  triad <- triad[!sapply(triad, is.null)]
  diffs <- diffs[!sapply(diffs, is.null)]
  aic <- aic[!is.na(aic)]
  bic <- bic[!is.na(bic)]
  print(ergm_results)
  print(diffs)
  print(length(diffs))
  
  ### Meta-analysis.
  ma_coeffs <- matrix(nrow = length(diffs), ncol = length(diffs[[1]]))
  ma_ses <- matrix(nrow = length(diffs), ncol = length(diffs[[1]]))
  ma_model <- list()
  for (i in 1 : length(diffs)){
    for (j in 1 : length(diffs[[1]])){
      ma_coeffs[i, ] <- ergm_results[[i]]$coef[1 : length(diffs[[1]])]
      ma_ses[i, ] <- sqrt(diag(ergm_results[[i]]$covar[1 : length(diffs[[1]]), 1 : length(diffs[[1]])]))
      suppressWarnings(ma_model[[j]] <- rma(yi = ma_coeffs[, j], sei = ma_ses[, j], control = list(stepadj = 0.5, maxiter = 1000)))
    }
  }
  
  ### Bias and rmse. 
  ## Bias.
  bias_matrix <- matrix(data = NA, nrow = length(diffs), ncol = length(diffs[[1]]))
  for (i in 1 : length(diffs)){
    bias_matrix[i, ] <- diffs[[i]]
  }
  bias <- colSums(bias_matrix) / length(diffs)
  
  ## rmse.
  rmse <- sqrt(colSums(bias_matrix ** 2) / length(diffs))
  
  ### Mean AIC and BIC.
  mean_aic <- mean(unlist(aic)) # Mean "meta-analysis" AIC.
  mean_bic <- mean(unlist(bic)) # Mean "meta-analysis" BIC.
  
  return(list(next_counter = next_counter,
              model_GOF = list(t = t, mahalm = mahalm), indegree_aux = list(tideg = tideg, ideg = ideg), 
              outdegree_aux = list(todeg = todeg, odeg = odeg), dist_aux = list(tdist = tdist, dist = dist), 
              esp_aux = list(tesp = tesp, esp = esp), dsp_aux = list(tdsp = tdsp, dsp = dsp), 
              triad_aux = list(ttriad = ttriad, triad = triad), aic = list(aic = unlist(aic), mean_aic = mean_aic), 
              bic = list(bic = unlist(bic), mean_bic = mean_bic), ma_model = ma_model, 
              bias = bias, rmse = rmse, seed1 = seed1, seed2 = seed2))
  
}

############################
##### Simulation study #####
############################

simulation_study_fun <- function(nsim){
  
  if (nsim %% 5 != 0){
    cat(paste0("Illegal 'nsim' input (not modular divisible).\n"))
    stop()
  }
  
  cat(paste0("\nInitializing condition het_0.\n"))
  
  ergm_sims_het0 <<- simulate(true_ergm_results_het0, nsim = nsim, set.seed = seed1)
  
  cov_results_het0 <- cov_sim_het0(nsim = nsim)
  
  gwodegree_results_het0 <- gwodegree_sim_het0(nsim = nsim)
  
  dtriad_results_het0 <- dtriad_sim_het0(nsim = nsim)
  
  cat(paste0("\nFinished with condition het_0. Initializing condition het_1.\n"))
  
  ergm_sims_het1 <<- vector(mode = "list", length = nsim)
  for (i in 1 : (0.8 * nsim)){
    ergm_sims_het1[[i]] <<- simulate(true_ergm_results_het0, nsim = 1, set.seed = seed1)  
  }
  for (i in (0.8 * nsim) : nsim){
    ergm_sims_het1[[i]] <<- simulate(true_ergm_results_het1, nsim = 1, set.seed = seed1)  
  }
  
  cov_results_het1 <- cov_sim_het1(nsim = nsim)
  
  gwodegree_results_het1 <- gwodegree_sim_het1(nsim = nsim)
  
  dtriad_results_het1 <- dtriad_sim_het1(nsim = nsim)
  
  cat(paste0("\nFinished with condition het_1. Returning output and terminating.\n"))
  
  return(list(cov_results_het0 = cov_results_het0, gwodegree_results_het0 = gwodegree_results_het0, 
              dtriad_results_het0 = dtriad_results_het0, cov_results_het1 = cov_results_het1,
              gwodegree_results_het1 = gwodegree_results_het1, dtriad_results_het1 = dtriad_results_het1))
  
}

cell_25  <- vector(mode = "list", length = 5)
cell_75  <- vector(mode = "list", length = 5)
save(cell_25, cell_75, file = "simdata5.RData")

timeR(
  for (i in 1 : 5){
    
    cat(paste0("\nThis is iteration ", i, ".\n"))
    
    cell_25[[i]]  <- simulation_study_fun(nsim = 25)
    cell_75[[i]]  <- simulation_study_fun(nsim = 75)
    
    resaveR(cell_25, cell_75, file = "simdata5.RData")
    
    cat(paste0("\nFinished with iteration ", i, ". Returning output.\n"))
    
  }
)


