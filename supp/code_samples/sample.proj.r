# Forecast model
# Since the forecast model can be implemented identically to a purely retrospective model
# as far as JAGS code goes, we are including R code as well to show the differences in implementation.
library(rjags)

# Settings
nRun <- 100000 # Number of times to run the MCMC, after burnin
nThin <- 25 # Thinning rate (only store 1 iteration out of every nThin)
nChain <- 3 # Number of separate chains to run MCMC (for testing convergence)
nBurn <- 10000 # Number of times to run the MCMC before getting output
nAdapt <- 1000 # Number of times to run the MCMC for adaptation at beginning
year1 <- 1966 # First year of estimation data
year2 <- 2010 # Last year of estimation data/first year of projection
year3 <- 2100 # Last year of projection
thresholds <- c(1,10,50) # Quasi-extirpation levels to monitor

# Read in input files
ydata <- as.matrix(read.csv('oven3a.csv', row.names=1))
observers <- as.matrix(read.csv('observers.csv', row.names=1))
wind.factor <- as.matrix(read.csv('wind.factor.csv', row.names=1))
# Fix formatting
wind.factor[wind.factor==" 0"] <- "0"
wind.factor[wind.factor==" 1"] <- "1"
wind.factor[wind.factor==" 2"] <- "2"

# More settings
species = 'oven'
aou <- 6740
nSites <- nrow(ydata)
nYears <- year3 - year1 + 1
nPast <- year2 - year1 + 1
nFuture <- year3 - year2

# File names to save
mod.name <- paste(species, 'sample.proj.JAGS', sep='.')
file.name <- paste(mod.name, 'RData', sep='.')
file.name2 <- paste(mod.name, 'gzip', sep='.')
file.name2

# For projection, fill in the future counts with all NA
ydata <- cbind(ydata, matrix(NA, nrow=nSites, ncol=nFuture))

fakeN <- ydata + 10 # Try to set up plausible initial values for N
fillNAs <- function(N, default) { # Need to replace NAs in initial N with values
  tMax <- ncol(N)
  defaults <- which(is.na(N[,1]) & is.na(N[,2]))
  N[defaults,1] <- default # Assign missing values from the 1st year which also have missing values for the 2nd year with default supplied
  subsequent <- which(is.na(N[,1]))
  N[subsequent,1] <- N[subsequent,2] # Assign missing values from the 1st year which don't have missing values for the 2nd year with the 2nd year's value
  for (t in 2:(tMax-1)) {
    # Assign missing values from this year which also have missing values for the next year with the values from the previous year
    previous <- which(is.na(N[,t]) & is.na(N[,t+1]))
    N[previous,t] <- N[previous,t-1]
    # Assign missing values from this year which don't have missing values for the next year with the mean of the previous and next year
    means <- which(is.na(N[,t]))
    N[means,t] <- ceiling((N[means, t-1] + N[means, t+1])/2)
  }
  # Assign missing values from the last year with the values from the previous year
  previous <- which(is.na(N[,tMax]))
  N[previous, tMax] <- N[previous, tMax - 1]
  return(N)
}
fakeN <- fillNAs(fakeN, 24)
head(fakeN)

# Replace missing observer numbers with a new value, calculate number of observers
obs.nr5 <- cbind(observers, matrix(NA, nrow=nSites, ncol=nFuture))
nObs <- max(as.vector(obs.nr5), na.rm=T) + 1
obs.nr5[is.na(obs.nr5)] = nObs
# Set wind speed category values, counting routes not run as wind speed 0
wind1 <- cbind(matrix(as.integer(!is.na(wind.factor) & wind.factor=="1"), nSites, nPast),
  matrix(0, nSites, nFuture))
wind2 <- cbind(matrix(as.integer(!is.na(wind.factor) & wind.factor=="2"), nSites, nPast),
  matrix(0, nSites, nFuture))
wind3 <- cbind(matrix(as.integer(!is.na(wind.factor) & wind.factor=="3+"), nSites, nPast),
  matrix(0, nSites, nFuture))

# Ricker + immigration model, NB initial abundance
# Logit-normal random observer effect
# Wind (factor) effects on p
# Projection model included
sink(file="dm.proj.oven1.txt")
cat("model {
lambda ~ dunif(0, 200) # Average initial abundance
alpha ~ dunif(0, 200) # Initial abundance overdispersion
# JAGS doesn't allow lambda based parameterization of NB, so reparameterize
P <- alpha/(alpha + lambda)
# Carrying capacity (K)
K ~ dunif(0, 200)
r ~ dunif(0, 2) # Instantaneous growth rate (return to equil. abundance)
iota ~ dunif(0, 10) # Immigration rate
# Detection probability parameters
p0 ~ dnorm(0, 0.01) # Intercept
sigma.p ~ dunif(0,10) # SD due to random observer effects
tau.p <- 1 / (sigma.p * sigma.p) # Convert SD to precision
p1 ~ dnorm(0, 0.01) # Effect of wind speed 1 on detection probability
p2 ~ dnorm(0, 0.01) # Effect of wind speed 2 on detection probability
p3 ~ dnorm(0, 0.01) # Effect of wind speed 3 or higher on detection probability

# Random observer effects, by observer
for (o in 1:nObs) {
  logitpObs[o] ~ dnorm(p0,tau.p)
}

for(i in 1:nSites) {
  N[i,1] ~ dnegbin(P, alpha) #Initial abundance
  # Detection probability
  logit(p[i,1]) <- logitpObs[obsMat[i,1]] + p1 * wind1[i,1] + p2 * wind2[i,1] +
    p3 * wind3[i,1]
  # Observed counts, based on true abundance and detection probabilities
  y[i,1] ~ dbin(p[i,1], N[i,1])
  for(t in 2:nYears) {
    # Expected abundance this time step (Ricker + immigration model)
    muN[i,t-1] <- N[i,t-1] * exp(r * (1 - N[i,t - 1] / K)) + iota
    # Due to demographic stochasticity, actual abundance differs from expected
    N[i,t] ~ dpois(muN[i,t-1])
    # Detection probability
    logit(p[i,t]) <- logitpObs[obsMat[i,t]] + p1 * wind1[i,t] +
      p2 * wind2[i,t] + p3 * wind3[i,t]
    # Observed counts, based on true abundance and detection probabilities
    y[i,t] ~ dbin(p[i,t], N[i,t])
  }
}

# If don't have enough memory to just trace N, instead trace derived variables
# of interest
for(t in 1:nYears) {
  meanN[t] <- mean(N[,t]) # Average by year
  for(j in 1:nThresholds) {
  # Proportion of sites (quasi-)extirpated by year (not cumulative over time)
    PQE[j,t] <- sum(N[,t]<thresholds[j])/nSites
  }
}

# If you're also interested in spatial patterns,
# snapshots for three years, but all sites
N1 <- N[,1]
N2 <- N[,nPast]
N3 <- N[,nYears]

}", fill=TRUE)
sink()

# Data to feed into JAGS
dat.proj.oven1 <- list(nSites=nSites, nYears=nYears, nPast=nPast,
  nObs=nObs, y=ydata, thresholds=thresholds, nThresholds=length(thresholds),
  obsMat = obs.nr5, wind1=wind1, wind2=wind2, wind3=wind3)
# Initial values for parameters
init.proj.oven1 <- function(chain) list(lambda=runif(1, 20, 75), alpha=runif(1, 0, 20),
  K=runif(1, 10, 50), iota=runif(1), r=runif(1, 0, 0.1),
  p0=rnorm(1), p1=rnorm(1), p2=rnorm(1), p3=rnorm(1), sigma.p=runif(1),
  N=fakeN, .RNG.name="base::Wichmann-Hill", .RNG.seed=year1 + aou + chain)
# Three options on what parameters to trace
pars.proj.oven1 <- c("lambda", "alpha", "K", "r", "iota", "p0", "p1", "p2",
  "p3", "sigma.p")
pars.proj.oven1.1 <- c("lambda", "alpha", "K", "r", "iota", "p0", "p1", "p2",
  "p3", "sigma.p", "meanN", "PQE", "N1", "N2", "N3")
pars.proj.oven1.2 <- c("lambda", "alpha", "K", "r", "iota", "p0", "p1", "p2",
  "p3", "sigma.p", "N")

set.seed(year1 + aou)
jm.proj.oven1 <- jags.model("dm.proj.oven1.txt", dat.proj.oven1,
  init.proj.oven1, n.chains=nChain, n.adapt=nAdapt)
update(jm.proj.oven1, n.iter=nBurn)
jc.proj.oven1 <- coda.samples(jm.proj.oven1, pars.proj.oven1.1, n.iter=nRun,
  thin=nThin)
summary(jc.proj.oven1)
save.image(file.name)
png(paste(mod.name, 'converge.\%i.png', sep='.'))
plot(jc.proj.oven1)
dev.off()
gelman.diag(jc.proj.oven1)
save(jc.proj.oven1, file=file.name2)