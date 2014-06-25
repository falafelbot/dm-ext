library(unmarked)

# Read in input files
ydata <- as.matrix(read.csv('oven3a.csv', row.names=1))
observers <- as.matrix(read.csv('observers.csv', row.names=1))
first.run <- as.matrix(read.csv('first.run.csv', row.names=1))
wind.factor <- as.matrix(read.csv('wind.factor.csv', row.names=1))
# Fix formatting
wind.factor[wind.factor==" 0"] <- "0"
wind.factor[wind.factor==" 1"] <- "1"
wind.factor[wind.factor==" 2"] <- "2"
# Combine obsCovs together into list
obs4 <- list(observers = observers, first.run = first.run, wind.factor = 
  wind.factor)
  
# Basic settings and labels
year1 <- 66
aou <- 'oven'
gam.model <- list(~1)[[1]] 
om.model <- list(~1)[[1]] 
lam.model <- list(~1)[[1]] 
K <- 600
iota.model <- list(~1) [[1]]

# Input data for pcountOpen
bird.frame <- unmarkedFramePCO(ydata, obsCovs=obs4, numPrimary=45)

# Run mixture models for initial abundance  
mods1 <- rep('', 3)
outputs1 = list()
for (i in 1:3) {
  mixture <- c('P','NB','ZIP')[i]
  dynamics <- c("constant", "autoreg", "trend", "ricker", "gompertz")[3]
  p.model <- list(~1, ~first.run, ~wind.factor, ~wind.factor + first.run)[[1]]
  immigration <- FALSE
  mod.name <- paste(gsub("\\*", ".", gsub(".stand", "", gsub("[ 1]", "", 
    paste(paste(aou, year1, sep=''), mixture, dynamics, 'lam', 
    substring(deparse(lam.model), 2), 'gam', substring(deparse(gam.model), 2), 
    'om', substring(deparse(om.model), 2), 'p', substring(deparse(p.model), 2), 
    ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
    'K', sep='.')))), K, sep="")
  mods1[i] <- mod.name
  file.name <- paste(mod.name, 'gzip', sep='.')
  print(file.name)
  fm <- pcountOpen(data=bird.frame, mixture=mixture, gammaformula=gam.model, 
    lambdaformula=lam.model, omegaformula=om.model, pformula=p.model, K=K, 
    dynamics=dynamics, immigration=immigration, iotaformula=iota.model)
  outputs1[[i]] <- fm
  print(summary(fm))
  save(fm, file=file.name)
}
names(outputs1) <- mods1

AIC.table <- function(mods, sort=T) {
  tab <- data.frame(model = names(mods), npar = NA, AIC = NA, delta = NA, 
    weight = NA)
  for (i in 1:length(mods)) {
    tab$npar[i] <- length(coef(mods[[i]]))
    tab$AIC[i] <- mods[[i]]@AIC
  }
  tab$delta <- tab$AIC - min(tab$AIC)
  tab$weight <- exp(-tab$delta/2)/sum(exp(-tab$delta/2))
  if (sort)
    tab <- tab[order(tab$AIC), ]
  return(tab)
}

# Compare these models
mix.aic <- AIC.table(outputs1)
mix.aic
mix.num <- as.integer(rownames(mix.aic)[1])
mix.num

# Run models for detection probability  
mods2 <- c(mods1[mix.num], rep('', 3))
outputs2 <- list(outputs1[[mix.num]])
for (i in 2:4) {
  mixture <- c('P','NB','ZIP')[mix.num]
  dynamics <- c("constant", "autoreg", "trend", "ricker", "gompertz")[3]
  p.model <- list(~1, ~first.run, ~wind.factor, ~wind.factor + first.run)[[i]]
  immigration <- FALSE
  mod.name <- paste(gsub("\\*", ".", gsub(".stand", "", gsub("[ 1]", "", 
    paste(paste(aou, year1, sep=''), mixture, dynamics, 'lam', 
    substring(deparse(lam.model), 2), 'gam', substring(deparse(gam.model), 2), 
    'om', substring(deparse(om.model), 2), 'p', substring(deparse(p.model), 2), 
    ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
    'K', sep='.')))), K, sep="")
  mods2[i] <- mod.name
  file.name <- paste(mod.name, 'gzip', sep='.')
  print(file.name)
  fm <- pcountOpen(data=bird.frame, mixture=mixture, gammaformula=gam.model, 
    lambdaformula=lam.model, omegaformula=om.model, pformula=p.model, K=K, 
    dynamics=dynamics, immigration=immigration, iotaformula=iota.model)
  outputs2[[i]] <- fm
  print(summary(fm))
  save(fm, file=file.name)
}
names(outputs2) <- mods2

# Compare these models
p.aic <- AIC.table(outputs2)
p.aic
p.num <- as.integer(rownames(p.aic)[1])
p.num


# Run models for population dynamics  
mods3 <- c(mods2[p.num], rep('', 8))
outputs3 <- list(outputs2[[p.num]])
# Constant dynamics
mixture <- c('P','NB','ZIP')[mix.num]
dynamics <- c("constant", "autoreg", "trend", "ricker", "gompertz")[1]
p.model <- list(~1, ~first.run, ~wind.factor, ~wind.factor + first.run)[[p.num]]
immigration <- FALSE
mod.name <- paste(gsub("\\*", ".", gsub(".stand", "", gsub("[ 1]", "", 
  paste(paste(aou, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'gam', substring(deparse(gam.model), 2), 
  'om', substring(deparse(om.model), 2), 'p', substring(deparse(p.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'K', sep='.')))), K, sep="")
mods3[2] <- mod.name
file.name <- paste(mod.name, 'gzip', sep='.')
print(file.name)
fm <- pcountOpen(data=bird.frame, mixture=mixture, gammaformula=gam.model, 
  lambdaformula=lam.model, omegaformula=om.model, pformula=p.model, K=K, 
  dynamics=dynamics, immigration=immigration, iotaformula=iota.model)
outputs3[[2]] <- fm
print(summary(fm))
save(fm, file=file.name)

# Autoregressive dynamics
mixture <- c('P','NB','ZIP')[mix.num]
dynamics <- c("constant", "autoreg", "trend", "ricker", "gompertz")[2]
p.model <- list(~1, ~first.run, ~wind.factor, ~wind.factor + first.run)[[p.num]]
immigration <- FALSE
mod.name <- paste(gsub("\\*", ".", gsub(".stand", "", gsub("[ 1]", "", 
  paste(paste(aou, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'gam', substring(deparse(gam.model), 2), 
  'om', substring(deparse(om.model), 2), 'p', substring(deparse(p.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'K', sep='.')))), K, sep="")
mods3[3] <- mod.name
file.name <- paste(mod.name, 'gzip', sep='.')
print(file.name)
fm <- pcountOpen(data=bird.frame, mixture=mixture, gammaformula=gam.model, 
  lambdaformula=lam.model, omegaformula=om.model, pformula=p.model, K=K, 
  dynamics=dynamics, immigration=immigration, iotaformula=iota.model)
outputs3[[3]] <- fm
print(summary(fm))
save(fm, file=file.name)

# Ricker dynamics
mixture <- c('P','NB','ZIP')[mix.num]
dynamics <- c("constant", "autoreg", "trend", "ricker", "gompertz")[4]
p.model <- list(~1, ~first.run, ~wind.factor, ~wind.factor + first.run)[[p.num]]
immigration <- FALSE
mod.name <- paste(gsub("\\*", ".", gsub(".stand", "", gsub("[ 1]", "", 
  paste(paste(aou, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'gam', substring(deparse(gam.model), 2), 
  'om', substring(deparse(om.model), 2), 'p', substring(deparse(p.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'K', sep='.')))), K, sep="")
mods3[4] <- mod.name
file.name <- paste(mod.name, 'gzip', sep='.')
print(file.name)
# Ricker and Gompertz generally seem to fail with default starting values, 
# so we set them
starts <- coef(outputs2[[p.num]])
starts <- c(starts[1], -2.35, starts[1]+0.3, starts[3:length(starts)])
fm <- pcountOpen(data=bird.frame, mixture=mixture, gammaformula=gam.model, 
  lambdaformula=lam.model, omegaformula=om.model, pformula=p.model, K=K, 
  dynamics=dynamics, immigration=immigration, iotaformula=iota.model, 
  starts=starts)
outputs3[[4]] <- fm
print(summary(fm))
save(fm, file=file.name)

# Gompertz dynamics
mixture <- c('P','NB','ZIP')[mix.num]
dynamics <- c("constant", "autoreg", "trend", "ricker", "gompertz")[5]
p.model <- list(~1, ~first.run, ~wind.factor, ~wind.factor + first.run)[[p.num]]
immigration <- FALSE
mod.name <- paste(gsub("\\*", ".", gsub(".stand", "", gsub("[ 1]", "", 
  paste(paste(aou, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'gam', substring(deparse(gam.model), 2), 
  'om', substring(deparse(om.model), 2), 'p', substring(deparse(p.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'K', sep='.')))), K, sep="")
mods3[5] <- mod.name
file.name <- paste(mod.name, 'gzip', sep='.')
print(file.name)
# Ricker and Gompertz generally seem to fail with default starting values, so we set them
starts <- coef(outputs2[[p.num]])
starts <- c(starts[1], -2.35, starts[1]+0.3, starts[3:length(starts)])
fm <- pcountOpen(data=bird.frame, gammaformula=gam.model, 
  lambdaformula=lam.model, omegaformula=om.model, pformula=p.model, K=K, 
  mixture=mixture, dynamics=dynamics, immigration=immigration, 
  iotaformula=iota.model, starts=starts)
outputs3[[5]] <- fm
print(summary(fm))
save(fm, file=file.name)

# Trend + immigration dynamics
mixture <- c('P','NB','ZIP')[mix.num]
dynamics <- c("constant", "autoreg", "trend", "ricker", "gompertz")[3]
p.model <- list(~1, ~first.run, ~wind.factor, ~wind.factor + first.run)[[p.num]]
immigration <- TRUE
mod.name <- paste(gsub("\\*", ".", gsub(".stand", "", gsub("[ 1]", "", 
  paste(paste(aou, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'gam', substring(deparse(gam.model), 2), 
  'om', substring(deparse(om.model), 2), 'p', substring(deparse(p.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'K', sep='.')))), K, sep="")
mods3[6] <- mod.name
file.name <- paste(mod.name, 'gzip', sep='.')
print(file.name)
fm <- pcountOpen(data=bird.frame, mixture=mixture, gammaformula=gam.model, 
  lambdaformula=lam.model, omegaformula=om.model, pformula=p.model, K=K, 
  dynamics=dynamics, immigration=immigration, iotaformula=iota.model)
outputs3[[6]] <- fm
print(summary(fm))
save(fm, file=file.name)

# Autoregressive + immigration dynamics
mixture <- c('P','NB','ZIP')[mix.num]
dynamics <- c("constant", "autoreg", "trend", "ricker", "gompertz")[2]
p.model <- list(~1, ~first.run, ~wind.factor, ~wind.factor + first.run)[[p.num]]
immigration <- TRUE
mod.name <- paste(gsub("\\*", ".", gsub(".stand", "", gsub("[ 1]", "", 
  paste(paste(aou, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'gam', substring(deparse(gam.model), 2), 
  'om', substring(deparse(om.model), 2), 'p', substring(deparse(p.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'K', sep='.')))), K, sep="")
mods3[7] <- mod.name
file.name <- paste(mod.name, 'gzip', sep='.')
print(file.name)
fm <- pcountOpen(data=bird.frame, mixture=mixture, gammaformula=gam.model, 
  lambdaformula=lam.model, omegaformula=om.model, pformula=p.model, K=K, 
  dynamics=dynamics, immigration=immigration, iotaformula=iota.model)
outputs3[[7]] <- fm
print(summary(fm))
save(fm, file=file.name)

# Ricker + immigration dynamics
mixture <- c('P','NB','ZIP')[mix.num]
dynamics <- c("constant", "autoreg", "trend", "ricker", "gompertz")[4]
p.model <- list(~1, ~first.run, ~wind.factor, ~wind.factor + first.run)[[p.num]]
immigration <- TRUE
mod.name <- paste(gsub("\\*", ".", gsub(".stand", "", gsub("[ 1]", "", 
  paste(paste(aou, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'gam', substring(deparse(gam.model), 2), 
  'om', substring(deparse(om.model), 2), 'p', substring(deparse(p.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'K', sep='.')))), K, sep="")
mods3[8] <- mod.name
file.name <- paste(mod.name, 'gzip', sep='.')
print(file.name)
# Ricker and Gompertz generally seem to fail with default starting values, 
# so we set them
starts <- coef(outputs2[[p.num]])
starts <- c(starts[1], -2.35, starts[1]+0.3, starts[3:(length(starts)-1)], 
  -1.5, starts[length(starts)])
fm <- pcountOpen(data=bird.frame, gammaformula=gam.model, 
  lambdaformula=lam.model, omegaformula=om.model, pformula=p.model, K=K, 
  mixture=mixture, dynamics=dynamics, immigration=immigration, 
  iotaformula=iota.model, starts=starts)
outputs3[[8]] <- fm
print(summary(fm))
save(fm, file=file.name)

# Gompertz + immigration dynamics
mixture <- c('P','NB','ZIP')[mix.num]
dynamics <- c("constant", "autoreg", "trend", "ricker", "gompertz")[5]
p.model <- list(~1, ~first.run, ~wind.factor, ~wind.factor + first.run)[[p.num]]
immigration <- TRUE
mod.name <- paste(gsub("\\*", ".", gsub(".stand", "", gsub("[ 1]", "", 
  paste(paste(aou, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'gam', substring(deparse(gam.model), 2), 
  'om', substring(deparse(om.model), 2), 'p', substring(deparse(p.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'K', sep='.')))), K, sep="")
mods3[9] <- mod.name
file.name <- paste(mod.name, 'gzip', sep='.')
print(file.name)
# Ricker and Gompertz generally seem to fail with default starting values, 
# so we set them
starts <- coef(outputs2[[p.num]])
starts <- c(starts[1], -2.35, starts[1]+0.3, starts[3:(length(starts)-1)], 
  -1.5, starts[length(starts)])
fm <- pcountOpen(data=bird.frame, gammaformula=gam.model, 
  lambdaformula=lam.model, omegaformula=om.model, pformula=p.model, K=K, 
  mixture=mixture, dynamics=dynamics, immigration=immigration, 
  iotaformula=iota.model, starts=starts)
outputs3[[9]] <- fm
print(summary(fm))
save(fm, file=file.name)

names(outputs3) <- mods3

# Compare these models
dyn.aic <- AIC.table(outputs3)
dyn.aic
dyn.num <- as.integer(rownames(dyn.aic)[1])
dyn.num
