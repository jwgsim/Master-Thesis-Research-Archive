#########################################################
########## Goodness-of-fit (GOF) demonstration ##########
#########################################################

#############################
##### Required Packages #####
#############################

### If necessary, install and require:  
## "network" package;
if(!require(network)) install.packages("network")
library(network) # Network data storage.
## "sna" package;
if(!require(sna)) install.packages("sna")
library(sna) # Network analysis routines.
## "ergm" package;
if(!require(ergm)) install.packages("ergm")
library(ergm) # Fitting & evaluating ERGMs.
## "coda" package;
if(!require(coda)) install.packages("coda")
library(coda) # MCMC diagnostics for goodness of fit.

##################################
##### Importing Network data #####
##################################

### Set working directory. 
setwd("Your_Working_Directory_Here")

### Read in the network data on the 19th network from the Vermeij data. Note that this is an empirical and not a simulated 
### network. The primary analysis strategy of the thesis project is to simulate a network which is similar in structure 
### to an empirical one (like the one presented here). This provides the advantage of knowing the "true" underlying parameter
### values, and will enable calculating quantities such as the bias when applying an ERGM to re-estimate these parameters. This 
### is then extended to a network sample context (but note that here we are exclusively operating on the single network level). 
### For the sake of convenience, an empirical network is here used, because the simulation setup is still a work in progress. 

## Import network file. 
not <- as.matrix(read.table(file = paste("C:\\Academia\\MSBBSS\\Jaar 2\\Thesis\\Thesis Data\\not", 19, ".net", sep = "")))
size <- dim(not)[1] 
not <- not[1 : size, 1 : size]
## Import covariate file.
cov <- read.table(file = paste("C:\\Academia\\MSBBSS\\Jaar 2\\Thesis\\Thesis Data\\cov", 19, ".dat", sep = ""),
                  col.names = c("girl", "identcGMC", "identcCMC", "identcCM", "percgirl")) 

### Transform the matrix into a network object and attach two of the actor variables to the network object:
notnet <- network(not, directed = T) # Specify as directed network.
notnet %v% "girl" <- cov$girl # Add sex actor variable. 
notnet %v% "identclass" <- cov$identcGMC # Add class membership actor variable.

###################################################################################
##### Exponential random graph model (ERGM) analysis of the classroom network #####
###################################################################################

### Specify the configuration of the ERGM. Note that both covariate and structural effects are included. 
model <- notnet ~ edges + mutual + nodeicov("identclass") + nodeocov("identclass") + nodematch("girl") + gwesp(.25,fixed = T) + m2star

### ERGM Estimation.
results <- ergm(model, control= control.ergm(MCMC.samplesize = 10000, seed = 3791), verbose = T) # Takes a short while. Note that 
# these parameter estimates and their standard errors will ultimately serve as inputs for the meta-analysed ERGM over a simulated
# network sample. Put differently, one ERGM with the same specification is fitted to each network in a sample (is done here for
# one such a network), and a meta-analysis is then used to obtain a single ERGM specification over that network sample.

### Summary of results.
summary(results)

##############################################################################
##### Calculating goodness-of-fit on the level of the individual network #####
##############################################################################

### Information criteria.
# Note that it is now possible to obtain the AIC and SIC for this ERGM over this network. 
summary(results)["aic"] # AIC.
summary(results)["bic"] # SIC.
# Note however that in their current form these values are not useful for diagnosing the GOF of this single network because they are 
# relative measures. This means that the same ERGM has to be applied to the remaining networks in the sample for the AIC and SIC to 
# become meaningful. In that case, the AIC and SIC can be used to evaluate the relative fit of the ERGM to each networks in the sample. 
# Note that different ERGMs can also be applied to either the same network or different networks in the sample, after which the AIC and
# SIC can be used to identify the superior model. For the sake of keeping the scope manageable, this possibility is excluded on the 
# single network level. 

# From these AIC and SIC values the AIC and SIC on the level of the network sample can subsequently also be obtained. The exact procedure
# is presented in the manuscript. Note however that the same problem that was defined earlier is here re-introduced, i.e., we need at least 
# two AIC and SIC at the level of the network sample to meaningfully interpret these values. As such, a differently specificied ERGM 
# should be applied to the sample to obtain sample level AIC and SIC values which can be compared to the previously obtained 
# AIC and SIC values. This comparison can then be used to identify the superior ERGM on the level of the network sample.

### Auxilary statistics GOF. 
## Start by quantifying whether the ERGM converged. The GOF function draws a sample of graphs for the parameters in the 
## specified ERGM, as described in the manuscript.
fit <- gof(model, GOF =~ model, coef = results$coef, control.gof.ergm(seed = 3792, MCMC.interval = 1000))
## It is then possible to compute t-values for each convergence statistic, by placing the observed value in the sample of graphs. 
meant <- apply(fit$sim.model, 2, mean) # Mean of each convergence statistic over the sample of graphs.
sdt <- apply(fit$sim.model, 2, sd) # Standard deviation of each convergence statistic over the sample of graphs.
obst <- summary(model) # Observed value for each convergence statistic.
t <- (meant - obst) / sdt # t-value for each convergence statistic.
## Compute mahalanobis distance over the fitted ERGM.
range <- which(sdt != 0) # To prevent singularity.
covm <- cov(as.matrix(fit$sim.model[, range]))
rangem <- range[which(eigen(covm)$values > 0.001)] # To prevent singularity.
covmodel <- cov(as.matrix(fit$sim.model[, rangem]))
mahalmodel <- mahalanobis(meant[rangem], fit$obs.model[rangem], covmodel) # Mahalanobis distance over the converged ERGM.
dfm <- length(rangem) # Degrees of freedom.
pm <- 1 - pchisq(mahalmodel, dfm) # p-value associated with Mahalanobis distance.
mahalm <- c(mahalmodel, dfm, pm) # Summary of three previous values.
names(mahalm)<-c("mahalmodel", "dfm", "pm") # Assign names.
# The t-values and the Mahalanobis distance are all close to zero. This is an indication that the ERGM converged. 

### Check model fit on auxiliary statistics. Here: indegree, outdegree, distance, edgewise-shared partners, 
### dyadwise shared partners, and triad census. To reiterate, this is the auxilary statistics GOF on the level of the 
### network as described in the manuscript. 
fitmodel <- gof(model, GOF =~ idegree + odegree + distance + espartners + dspartners + triadcensus,
                coef = results$coef, control.gof.formula(seed = 3792, nsim = 1000))

### These can be plotted to inspect visually whether the ERGM fits to the network. 
plot(fitmodel) 
# It is clear that the indegree is not fit well by the model. Note however that it is not very clear whether the 
# remaining paramters do (mostly because of a poorly positioned y-axis). 

### As was done for the convergence statistics, it is now similarly possible to calculate t-values and Mahalanobis distances 
### for each auxilary statistic over the sample of graphs. This provides a more formal way (relative to visualizing) of 
### quantifying the GOF of the ERGM to the network.

## Indegree.
# t-value.
meantideg <- apply(fitmodel$sim.ideg, 2, mean)
sdtideg <- apply(fitmodel$sim.ideg, 2, sd)
tideg <- (meantideg - fitmodel$obs.ideg) / sdtideg
# Mahalanobis-distance.
range <- which(sdtideg != 0) 
covc <- cov(as.matrix(fitmodel$sim.ideg[, range])) 
rangei <- range[which(eigen(covc)$values > 0.001)] 
covideg <- cov(as.matrix(fitmodel$sim.ideg[, rangei]))
mahalideg <- mahalanobis(meantideg[rangei], fitmodel$obs.ideg[rangei], covideg)
dfi <- length(rangei)
pi <- 1 - pchisq(mahalideg,dfi)
ideg <- c(mahalideg ,dfi, pi)
names(ideg) <- c("mahalideg", "dfi", "pi")

## Outdegree.
# t-value.
meantodeg <- apply(fitmodel$sim.odeg, 2, mean)
sdtodeg <- apply(fitmodel$sim.odeg, 2, sd)
todeg <- (meantodeg - fitmodel$obs.odeg) / sdtodeg
# Mahalanobis-distance.
range <- which(sdtodeg != 0) 
covc <- cov(as.matrix(fitmodel$sim.odeg[, range]))
rangeo <- range[which(eigen(covc)$values > 0.001)]
covodeg <- cov(as.matrix(fitmodel$sim.odeg[, rangeo]))
mahalodeg <- mahalanobis(meantodeg[rangeo], fitmodel$obs.odeg[rangeo], covodeg)
dfo <- length(rangeo)
po <- 1 - pchisq(mahalodeg, dfo)
odeg <- c(mahalodeg, dfo, po)
names(odeg) <- c("mahalodeg", "dfo", "po")

## Geodesic distance.
# t-value.
meantdist <- apply(fitmodel$sim.dist, 2, mean)
sdtdist <- apply(fitmodel$sim.dist, 2, sd)
tdist <- (meantdist - fitmodel$obs.dist) / sdtdist
# Mahalanobis-distance.
range <- which(sdtdist != 0) 
covc <- cov(as.matrix(fitmodel$sim.dist[, range]))
ranged <- range[which(eigen(covc)$values > 0.001)] 
covdist <- cov(as.matrix(fitmodel$sim.dist[, ranged]))
mahaldist <- mahalanobis(meantdist[ranged], fitmodel$obs.dist[ranged], covdist)
dfdist <- length(ranged)
pdist <- 1 - pchisq(mahaldist, dfdist)
dist <- c(mahaldist, dfdist, pdist)
names(dist) <- c("mahaldist", "dfdist", "pdist")

## Edgewise shared partners.
# t-value.
meantesp <- apply(fitmodel$sim.esp, 2, mean)
sdtesp <- apply(fitmodel$sim.esp, 2, sd)
tesp <- (meantesp - fitmodel$obs.esp) / sdtesp
# Mahalanobis-distance.
range <- which(sdtesp != 0) 
covc <- cov(as.matrix(fitmodel$sim.esp[, range]))
rangee <- range[which(eigen(covc)$values > 0.001)]
covesp <- cov(as.matrix(fitmodel$sim.esp[, rangee]))
mahalesp <- mahalanobis(meantesp[rangee], fitmodel$obs.esp[rangee], covesp)
dfesp <- length(rangee)
pesp <- 1 - pchisq(mahalesp, dfesp)
esp <- c(mahalesp, dfesp, pesp)

## Dyadwise shared partners
# t-value.
meantdsp <- apply(fitmodel$sim.dsp, 2, mean)
sdtdsp <- apply(fitmodel$sim.dsp, 2, sd)
tdsp <- c(meantdsp - fitmodel$obs.dsp) / sdtdsp
# Mahalanobis-distance.
range <- which(sdtdsp != 0) 
covc <- cov(as.matrix(fitmodel$sim.dsp[, range]))
rangeds <- range[which(eigen(covc)$values > 0.001)] 
covdsp <- cov(as.matrix(fitmodel$sim.dsp[, rangeds]))
mahaldsp <- mahalanobis(meantdsp[rangeds], fitmodel$obs.dsp[rangeds], covdsp)
dfdsp <- length(rangeds)
pdsp <- 1 - pchisq(mahaldsp, dfdsp)
dsp <- c(mahaldsp, dfdsp, pdsp)
names(dsp) <- c("mahaldsp", "dfdsp", "pdsp")

## Triadcensus
# t-value.
meanttriad <- apply(fitmodel$sim.triad, 2, mean)
sdttriad <- apply(fitmodel$sim.triad, 2, sd)
ttriad <- (meanttriad - fitmodel$obs.triad) / sdttriad
# Mahalanobis-distance.
range <- which(sdttriad != 0) 
covc <- cov(as.matrix(fitmodel$sim.triad[, range]))
ranget <- range[which(eigen(covc)$values > 0.001)] 
covtriad <- cov(as.matrix(fitmodel$sim.triad[, ranget]))
mahaltriad <- mahalanobis(meanttriad[ranget], fitmodel$obs.triad[ranget], covtriad)
dftriad <- length(ranget)
ptriad <- 1 - pchisq(mahaltriad, dftriad)
triad <- c(mahaltriad, dftriad, ptriad)

######################
##### Conclusion #####
######################

### The previous can be summarized by way of the Mahalanobis p-values: 
summ <- round(rbind(pi, po, pdist, pesp, pdsp, ptriad), 2)
rownames(summ) <- c("p-indegree", "p-outdegree", "p-distance", "p-edgewise", "p-dyadwise", "p-triad-census")
summ

# It becomes clear from the Mahalanobis p-values that the indegree, dyadwise shared partners, and triad census statistics are not 
# well-modeled by the ERGM. The outdegree, distance, and edgewise shared partners are. 

# In summary, we have seen how to determine the goodness-of-fit (GOF) of an exponential random graph model (ERGM) to a single network 
# with information criteria as an index (with some additional textual explanation), or auxilary statistics as an index. 

# Note that GOF on the single network is readily extended to the meta-analyzed ERGM at the network sample level. It was already discussed 
# how the AIC and SIC can be extended to the network level. With respect to the auxilary statistics approach, a similar extension from the 
# individual network level to the sample network level can be made. After having obtained the meta-analyzed ERGM from the parameter estimates
# and standard errors from the set of individual networks in the sample, a similar procedure as outlined earlier can applied. More specifically,
# t-values and the Mahalanobis distance can be applied at the level of the meta-analyzed ERGM to inspect how well it fits to the sample as a
# whole. Doing so was unfortunately not possible here, because the implementation of the meta-analyzed ERGM and the associated calculations 
# for obtaining auxilary goodness-of-fit are still a work-in-progress. 
