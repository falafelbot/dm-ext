# BBS analysis code using JAGS
# Read in input files
ydata <- as.matrix(read.csv('oven3a.csv', row.names=1))
obs5 <- as.matrix(read.csv('observers.csv', row.names=1))
first <- as.matrix(read.csv('first.run.csv', row.names=1))
wind.factor <- as.matrix(read.csv('wind.factor.csv', row.names=1))
# Fix formatting
wind.factor[wind.factor==" 0"] <- "0"
wind.factor[wind.factor==" 1"] <- "1"
wind.factor[wind.factor==" 2"] <- "2"

library(rjags)
library(runjags)

# Settings
year1 <- 66
species <- 'oven'
nSites <- nrow(ydata)
nYears <- ncol(ydata)


fillN3 <- function(N, p, lambda, alpha) { # Need to replace NAs in initial N with values, adjust others by detection probability
  for (i in 1:nrow(N)) {
    notNA <- which(!is.na(N[i, ]))
    size <- N[i, notNA] 
    size[size==0] <- 0.5
    addTo <- rnbinom(length(notNA), size=size, prob=p)
    N[i, notNA] <- N[i, notNA] + addTo
    for (t in 1:ncol(N))
      if (is.na(N[i,t])) {
        if (t == ncol(N) || is.na(N[i, t+1])) {
          if (t==1) {
            N[i, t] <- rnbinom(1, size=alpha, mu=lambda) #If no previous or subsequent value, use supplied default
          }
          else {
            N[i, t] <- N[i, t-1] #If no subsequent value, use previous
          }
        }
        else if (t==1)
          N[i, t] <- N[i, t+1] #If no previous value, use subsequent
        else
          N[i, t] <- ceiling((N[i, t-1] + N[i, t+1])/2) #Otherwise, use mean
      }
  }
  return(N)
}

N <- fillN3(ydata, 0.18, 33.4, 0.65)
head(N)
summary(N)

# Other input variables
nObs <- max(as.vector(obs5), na.rm=T) + 1
obs5[is.na(obs5)] <- nObs
wind1 <- matrix(as.integer(!is.na(wind.factor) & wind.factor=="1"), nSites, nYears)
summary(wind1)
wind2 <- matrix(as.integer(!is.na(wind.factor) & wind.factor=="2"), nSites, nYears)
wind3 <- matrix(as.integer(!is.na(wind.factor) & wind.factor=="3+"), nSites, nYears)
first[is.na(first)] <- 1

# Repeat top model from frequentist analysis
# Pick filename based on model
mixture <- c('P','NB','ZIP')[2]
dynamics <- c("constant", "autoreg", "trend", "ricker", "gompertz")[4]
p.model <- list(~1, ~first.run, ~wind.factor, ~wind.factor + first.run)[[4]] 
r.model <- list(~1)[[1]] 
K.model <- list(~1)[[1]] 
lam.model <- list(~1)[[1]] 
immigration <- TRUE
iota.model <- list(~1) [[1]]
p.rand <- c("None","Normal","Beta")[1]
mod.name <- paste(paste(species, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'r', substring(deparse(r.model), 2), 
  'K', substring(deparse(K.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'p', substring(deparse(p.model), 2), 'rand', p.rand, 'JAGS', sep='.') 
file.name <- paste(mod.name, 'RData', sep='.')
file.name

sink(file="dm.ricki1.txt")
cat("
model {
lambda ~ dunif(0, 200)
alpha ~ dunif(0, 200)
p.nb <- alpha/(alpha + lambda)
r ~ dunif(0, 5)
K ~ dunif(0, 200)
p0 ~ dnorm(0, 0.01)
p1 ~ dnorm(0, 0.01)
p2 ~ dnorm(0, 0.01)
p3 ~ dnorm(0, 0.01)
p1st ~ dnorm(0, 0.01)
iota ~ dunif(0, 15)
for(i in 1:nSites) {
  N[i,1] ~ dnegbin(p.nb, alpha)
  for(t in 2:nYears) {
      muN[i,t-1] <- N[i,t-1]*exp(r*(1-N[i,t-1]/K))+iota
      N[i,t] ~ dpois(muN[i,t-1])
    }
  for(t in 1:nYears) {
    logit(p[i,t]) <- p0 + p1*wind1[i,t] + p2*wind2[i,t] + p3*wind3[i,t] + p1st*first[i,t]
    y[i,t] ~ dbin(p[i,t], N[i,t])
    }
  }
}
", fill=TRUE)
sink()

dat.ricki1 <- list(nSites=nSites, nYears=nYears, y=ydata, wind1=wind1, 
  wind2=wind2, wind3=wind3, first=first)
init.ricki1 <- function() list(lambda=runif(1, 10, 60), alpha=runif(1, 0, 60), 
  r=runif(1), K=runif(1, 30, 120), iota=runif(1),
  p0=rnorm(1, -1.5), p1=rnorm(1), p2=rnorm(1), p3=rnorm(1), p1st=rnorm(1), N=N)
pars.ricki1 <- c("lambda", "alpha", "r", "K", "p0", "p1", "p2", "p3", "p1st", "iota")

jm.ricki1 <- jags.model("dm.ricki1.txt", dat.ricki1, init.ricki1, n.chains=5, n.adapt=5000)
jc.ricki1 <- coda.samples(jm.ricki1, pars.ricki1, n.iter=100000)

summary(jc.ricki1)
windows(record=T)
plot(jc.ricki1)
warnings()
gelman.diag(jc.ricki1)

jc.ricki1a <- window(jc.ricki1, 51001)
summary(jc.ricki1a)
gelman.diag(jc.ricki1a)

save.image(file.name)

# Add random observer effects
rm(jc.ricki1)
rm(jc.ricki1a)
p.rand <- c("None","Normal","Beta")[2]
mod.name <- paste(paste(species, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'r', substring(deparse(r.model), 2), 
  'K', substring(deparse(K.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'p', substring(deparse(p.model), 2), 'rand', p.rand, 'JAGS', sep='.') 
file.name = paste(mod.name, 'RData', sep='.')
file.name

sink(file="dm.ricki2.txt")
cat("
model {
lambda ~ dunif(0, 200)
alpha ~ dunif(0, 200)
p.nb <- alpha/(alpha + lambda)
r ~ dunif(0, 5)
K ~ dunif(0, 200)
p0 ~ dnorm(0, 0.01)
p1 ~ dnorm(0, 0.01)
p2 ~ dnorm(0, 0.01)
p3 ~ dnorm(0, 0.01)
p1st ~ dnorm(0, 0.01)
iota ~ dunif(0, 15)
sigma.p ~ dunif(0,15)
tau.p <- 1 / (sigma.p * sigma.p)

for (o in 1:nObs) {
  eta[o] ~ dnorm(0,tau.p)
}

for(i in 1:nSites) {
  N[i,1] ~ dnegbin(p.nb, alpha)
  for(t in 2:nYears) {
      muN[i,t-1] <- N[i,t-1]*exp(r*(1-N[i,t-1]/K))+iota
      N[i,t] ~ dpois(muN[i,t-1])
    }
  for(t in 1:nYears) {
    logit(p[i,t]) <- p0 + p1*wind1[i,t] + p2*wind2[i,t] + p3*wind3[i,t] + p1st*first[i,t] + eta[obsID[i,t]]
    y[i,t] ~ dbin(p[i,t], N[i,t])
    }
  }
}
", fill=TRUE)
sink()

dat.ricki2 <- list(nSites=nSites, nYears=nYears, y=ydata, wind1=wind1, 
  wind2=wind2, wind3=wind3, first=first, nObs=nObs, obsID=obs5)
init.ricki2 <- function() list(lambda=runif(1, 10, 60), alpha=runif(1, 0, 60), 
  r=runif(1, 0, 0.1), K=runif(1, 40, 70),
  p0=rnorm(1, -1.5), p1=rnorm(1), p2=rnorm(1), p3=rnorm(1), p1st=rnorm(1), 
  iota=runif(1), sigma.p=runif(1, 0, 2), N=N)
pars.ricki2 <- c("lambda", "alpha", "r", "K", "iota", "p0", "p1", "p2", "p3", 
  "p1st", "sigma.p")

jm.ricki2 <- jags.model("dm.ricki2.txt", dat.ricki2, init.ricki2, n.chains=5, n.adapt=3000)
jc.ricki2 <- coda.samples(jm.ricki2, pars.ricki2, n.iter=40000)

summary(jc.ricki2)
windows(record=T)
plot(jc.ricki2)
warnings()

save.image(file.name)

# Regional environmental stochasticity
rm(jc.ricki2)
p.rand <- c("None","Normal","Beta")[1]
mod.name <- paste(paste(species, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'r', substring(deparse(r.model), 2), 
  'K', substring(deparse(K.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'p', substring(deparse(p.model), 2), 'rand', p.rand, 'JAGS.ES', sep='.') 
file.name = paste(mod.name, 'RData', sep='.')
file.name

sink(file="dm.ricki3.txt")
cat("
model {
lambda ~ dunif(0, 400)
alpha ~ dunif(0, 200)
p.nb <- alpha/(alpha + lambda)
r ~ dunif(0, 3)
sigma.nu ~ dunif(0, 10)
tau.nu <- 1 / (sigma.nu * sigma.nu)
K ~ dunif(0, 400)
p0 ~ dnorm(0, 0.01)
p1 ~ dnorm(0, 0.01)
p2 ~ dnorm(0, 0.01)
p3 ~ dnorm(0, 0.01)
p1st ~ dnorm(0, 0.01)
iota ~ dunif(0, 15)

for(t in 2:nYears) {
  nu[t-1] ~ dnorm(0, tau.nu)
}


for(i in 1:nSites) {
  N[i,1] ~ dnegbin(p.nb, alpha)
  for(t in 2:nYears) {
      muN[i,t-1] <- N[i,t-1]*exp(nu[t-1]+r*(1-N[i,t-1]/K))+iota
      N[i,t] ~ dpois(muN[i,t-1])
    }
  for(t in 1:nYears) {
    logit(p[i,t]) <- p0 + p1*wind1[i,t] + p2*wind2[i,t] + p3*wind3[i,t] + p1st*first[i,t]
    y[i,t] ~ dbin(p[i,t], N[i,t])
    }
  }
}
", fill=TRUE)
sink()

dat.ricki3 <- list(nSites=nSites, nYears=nYears, y=ydata, wind1=wind1, 
  wind2=wind2, wind3=wind3, first=first) 
init.ricki3 <- function() list(lambda=runif(1, 10, 60), alpha=runif(1, 0, 60), 
  r=runif(1, 0, 0.2), sigma.nu=runif(1), K=runif(1, 30, 120), iota=runif(1),
  p0=rnorm(1, -1.5), p1=rnorm(1), p2=rnorm(1), p3=rnorm(1), p1st=rnorm(1), N=N)
pars.ricki3 <- c("lambda", "alpha", "r", "sigma.nu", "K", "iota", 
  "p0", "p1", "p2", "p3", "p1st")

jm.ricki3 <- jags.model("dm.ricki3.txt", dat.ricki3, init.ricki3, n.chains=5, n.adapt=1000)
update(jm.ricki3, 50000)
jc.ricki3 <- coda.samples(jm.ricki3, pars.ricki3, n.iter=130000)

summary(jc.ricki3)
save.image(file.name)
HPDinterval(jc.ricki3)
gelman.diag(jc.ricki3)
png(paste(mod.name, '%i.png', sep='.'))
plot(jc.ricki3)
dev.off()

# Regional environmental stochasticity AND random observer effects
summary(jc.ricki3)
p.rand <- c("None","Normal","Beta")[2]
mod.name <- paste(paste(species, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'r', substring(deparse(r.model), 2), 
  'K', substring(deparse(K.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'p', substring(deparse(p.model), 2), 'rand', p.rand, 'JAGS.ES', sep='.') 
file.name <- paste(mod.name, 'RData', sep='.')
file.name

sink(file="dm.ricki4.txt")
cat("
model {
lambda ~ dunif(0, 400)
alpha ~ dunif(0, 200)
p.nb <- alpha/(alpha + lambda)
r ~ dunif(0, 3)
sigma.nu ~ dunif(0, 10)
tau.nu <- 1 / (sigma.nu * sigma.nu)
K ~ dunif(0, 1000)
p0 ~ dnorm(0, 0.01)
p1 ~ dnorm(0, 0.01)
p2 ~ dnorm(0, 0.01)
p3 ~ dnorm(0, 0.01)
p1st ~ dnorm(0, 0.01)
iota ~ dunif(0, 15)
sigma.p ~ dunif(0,15)
tau.p <- 1 / (sigma.p * sigma.p)

for (o in 1:nObs) {
  eta[o] ~ dnorm(0,tau.p)
}

for(t in 2:nYears) {
  nu[t-1] ~ dnorm(0, tau.nu)
}


for(i in 1:nSites) {
  N[i,1] ~ dnegbin(p.nb, alpha)
  for(t in 2:nYears) {
      muN[i,t-1] <- N[i,t-1]*exp(nu[t-1]+r*(1-N[i,t-1]/K))+iota
      N[i,t] ~ dpois(muN[i,t-1])
    }
  for(t in 1:nYears) {
    logit(p[i,t]) <- p0 + p1*wind1[i,t] + p2*wind2[i,t] + p3*wind3[i,t] + p1st*first[i,t] + eta[obsID[i,t]]
    y[i,t] ~ dbin(p[i,t], N[i,t])
    }
  }
}
", fill=TRUE)
sink()

dat.ricki4 <- list(nSites=nSites, nYears=nYears, y=ydata, wind1=wind1, 
  wind2=wind2, wind3=wind3, nObs=nObs, obsID=obs5, first=first)
init.ricki4 <- function() list(lambda=runif(1, 10, 60), alpha=runif(1, 0, 60), 
  r=runif(1, 0, 0.2), sigma.nu=runif(1), K=runif(1, 30, 120), iota=runif(1),
  p0=rnorm(1, -1.5), p1=rnorm(1), p2=rnorm(1), p3=rnorm(1), p1st=rnorm(1), sigma.p=runif(1), N=N)
pars.ricki4 <- c("lambda", "alpha", "r", "sigma.nu", "K", "iota", 
  "p0", "p1", "p2", "p3", "p1st", "sigma.p")

jm.ricki4 <- jags.model("dm.ricki4.txt", dat.ricki4, init.ricki4, n.chains=5, n.adapt=1000)
update(jm.ricki4, 50000)
jc.ricki4 <- coda.samples(jm.ricki4, pars.ricki4, n.iter=130000)

summary(jc.ricki4)
save.image(file.name)
HPDinterval(jc.ricki4)
gelman.diag(jc.ricki4)

jpeg(paste(mod.name, '%i.jpg', sep='.'))
plot(jc.ricki4)
dev.off()

# Track population sizes (includes random observer effects but not environmental stochasticity)
rm(jc.ricki4)
p.rand <- c("None","Normal","Beta")[2]
mod.name <- paste(paste(species, year1, sep=''), mixture, dynamics, 'lam', 
  substring(deparse(lam.model), 2), 'r', substring(deparse(r.model), 2), 
  'K', substring(deparse(K.model), 2), 
  ifelse(immigration, "iotaT", "iotaF"), substring(deparse(iota.model), 2),
  'p', substring(deparse(p.model), 2), 'rand', p.rand, 'JAGS.N', sep='.') 
file.name <- paste(mod.name, 'RData', sep='.')
file.name

sink(file="dm.ricki5.txt")
cat("
var Ndecade[nSites,5]
model {
lambda ~ dunif(0, 200)
alpha ~ dunif(0, 200)
p.nb <- alpha/(alpha + lambda)
r ~ dunif(0, 5)
K ~ dunif(0, 200)
p0 ~ dnorm(0, 0.01)
p1 ~ dnorm(0, 0.01)
p2 ~ dnorm(0, 0.01)
p3 ~ dnorm(0, 0.01)
p1st ~ dnorm(0, 0.01)
iota ~ dunif(0, 15)
sigma.p ~ dunif(0,15)
tau.p <- 1 / (sigma.p * sigma.p)

for (o in 1:nObs) {
  eta[o] ~ dnorm(0,tau.p)
}

for(i in 1:nSites) {
  N[i,1] ~ dnegbin(p.nb, alpha)
  for(t in 2:nYears) {
      muN[i,t-1] <- N[i,t-1]*exp(r*(1-N[i,t-1]/K))+iota
      N[i,t] ~ dpois(muN[i,t-1])
    }
  for(t in 1:nYears) {
    logit(p[i,t]) <- p0 + p1*wind1[i,t] + p2*wind2[i,t] + p3*wind3[i,t] + p1st*first[i,t] + eta[obsID[i,t]]
    y[i,t] ~ dbin(p[i,t], N[i,t])
    }
  }

for(t in 1:nYears) {
  Nmean[t] <- mean(N[,t])
}

for (t2 in 1:length(decades)) {
  for(i in 1:nSites) {
    Ndecade[i,t2] <- N[i,decades[t2]]
  }
}
}
", fill=TRUE)
sink()

dat.ricki5 <- list(nSites=nSites, nYears=nYears, y=ydata, wind1=wind1, 
  wind2=wind2, wind3=wind3, nObs=nObs, obsID=obs5, first=first,
  decades=as.integer(seq(5, 45, 10)))
init.ricki5 <- function(chain) list(lambda=runif(1, 10, 60), alpha=runif(1, 0, 60), 
  r=runif(1), K=runif(1, 30, 120), iota=runif(1),
  p0=rnorm(1, -1.5), p1=rnorm(1), p2=rnorm(1), p3=rnorm(1), p1st=rnorm(1), sigma.p=runif(1), 
  N=N, .RNG.seed=chain, .RNG.name = c("base::Wichmann-Hill","base::Marsaglia-Multicarry",
  "base::Super-Duper","base::Mersenne-Twister")[chain])
pars.ricki5 <- c("lambda", "alpha", "r", "K", "p0", "p1", "p2", "p3", "p1st", 
  "sigma.p", "iota", "Nmean", "Ndecade")

jr.ricki5 <- run.jags("dm.ricki5.txt", data=dat.ricki5, inits=init.ricki5, 
  n.chains=4, adapt=1000, burnin=99000, sample=40000, monitor=pars.ricki5, 
  thin=2, method="parallel")

print(jr.ricki5)
save.image(file.name)

jc.ricki5 <- as.mcmc(jr.ricki5)
mat.ricki5 <- as.matrix(jc.ricki5)
NmeanCols<- which(colnames(mat.ricki5)=="Nmean[1]") + 1:nYears - 1
Nmean.frame<-data.frame(year=1966:2010, Mean=colMeans(mat.ricki5[,NmeanCols]), HPDinterval(jc.ricki5[,NmeanCols]))
Nmean.frame
NdecCols <- which(colnames(mat.ricki5)=="Ndecade[1,1]") + 1:(5*nSites) - 1
mean.N1 <- colMeans(mat.ricki5[,NdecCols])
mean.N1
mean.N2 <- matrix(mean.N1, nrow=122, ncol=5)
mean.N2
colnames(mean.N2) <- paste("N", seq(1970,2010,10), sep="")
max(mean.N2)
max.plot <- round(max(mean.N2)/5, -1)*5
Ndec.frame <- with(rt2a, data.frame(statenum, Route, Lati, Longi, mean.N2))
head(Ndec.frame)


# Measures of variability
sd.N1 <- apply(mat.ricki5[,NdecCols], 2, sd)
sd.N1
mean(sd.N1)
mean(mean.N1)
mean(mat.ricki5[,"lambda"])
sd(mat.ricki5[,"lambda"])
sd(mat.ricki5[,"lambda"])/mean(mat.ricki5[,"lambda"])
sd(mat.ricki5[,"alpha"])/mean(mat.ricki5[,"alpha"])
sd(mat.ricki5[,"iota"])/mean(mat.ricki5[,"iota"])
sd(mat.ricki5[,"r"])/mean(mat.ricki5[,"r"])
summary(sd.N1/mean.N1)
sd.N2 <- apply(mean.N2, 2, sd) #Variation between sites in a year
sd.N2

# Make plot
library(ggplot2)
library(gridExtra)
library(maps)

all.states <- map_data("state")
states <- subset(all.states, region %in% c('maryland', 'virginia')) 

max.radius <- 4
leg.sizes <- c(1, 10, 100, 250)
gN.sp <- vector("list", 5)
gN.sp[[1]] <- ggplot(Ndec.frame, aes(x=Longi,y=Lati)) + theme_bw() + labs(x="", y="")+
  geom_polygon(data=states, aes(x=long, y=lat, group=group), color="black", fill="white") +
  geom_point(aes(size=N1970), alpha=0.5, color="black")+scale_size_area(name="Abundance", max_size=max.radius, breaks=leg.sizes)+
  geom_text(x=-83.2, y=39.75, label="1970")
gN.sp[[2]] <- ggplot(Ndec.frame, aes(x=Longi,y=Lati)) + theme_bw() + labs(x="", y="")+
  geom_polygon(data=states, aes(x=long, y=lat, group=group), color="black", fill="white") +
  geom_point(aes(size=N1980), alpha=0.5, color="black")+scale_size_area(name="Abundance", max_size=max.radius, breaks=leg.sizes)+
  geom_text(x=-83.2, y=39.75, label="1980")
gN.sp[[3]] <- ggplot(Ndec.frame, aes(x=Longi,y=Lati)) + theme_bw() + labs(x="", y="")+
  geom_polygon(data=states, aes(x=long, y=lat, group=group), color="black", fill="white") +
  geom_point(aes(size=N1990), alpha=0.5, color="black")+scale_size_area(name="Abundance", max_size=max.radius, breaks=leg.sizes)+
  geom_text(x=-83.2, y=39.75, label="1990")
gN.sp[[4]] <- ggplot(Ndec.frame, aes(x=Longi,y=Lati)) + theme_bw() + labs(x="", y="")+
  geom_polygon(data=states, aes(x=long, y=lat, group=group), color="black", fill="white") +
  geom_point(aes(size=N2000), alpha=0.5, color="black")+scale_size_area(name="Abundance", max_size=max.radius, breaks=leg.sizes)+
  geom_text(x=-83.2, y=39.75, label="2000")
gN.sp[[5]] <- ggplot(Ndec.frame, aes(x=Longi,y=Lati)) + theme_bw() + labs(x="", y="")+
  geom_polygon(data=states, aes(x=long, y=lat, group=group), color="black", fill="white") +
  geom_point(aes(size=N2010), alpha=0.5, color="black")+scale_size_area(name="Abundance", max_size=max.radius, breaks=leg.sizes)+
  geom_text(x=-83.2, y=39.75, label="2010")

gN.mean <- ggplot(Nmean.frame, aes(x=year,y=Mean)) + theme_bw() + labs(x="Year",
  y="Mean Abundance")+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="black",
  alpha=0.5,colour=NA) + geom_line(colour="black", size=1) 

png(filename = "OVEN_N_by_route_year8.png", height = 6.5, width = 6.5, units = "in", res=600)
print(arrangeGrob(gN.sp[[1]] + theme(legend.position=c(0.145, 0.56), 
  legend.background=element_blank(), axis.text=element_blank(), 
  axis.line=element_blank(), axis.ticks=element_blank(), panel.border=element_blank(),
  panel.grid=element_blank(), plot.margin=unit(rep(0.1, 4), "lines")),
  gN.sp[[2]] + theme(legend.position="none", axis.text=element_blank(), 
  axis.line=element_blank(), axis.ticks=element_blank(), panel.border=element_blank(),
  panel.grid=element_blank(), plot.margin=unit(rep(0.1, 4), "lines")),
  gN.sp[[3]] + theme(legend.position="none", axis.text=element_blank(), 
  axis.line=element_blank(), axis.ticks=element_blank(), panel.border=element_blank(),
  panel.grid=element_blank(), plot.margin=unit(rep(0.1, 4), "lines")),
  gN.sp[[4]] + theme(legend.position="none", axis.text=element_blank(), 
  axis.line=element_blank(), axis.ticks=element_blank(), panel.border=element_blank(),
  panel.grid=element_blank(), plot.margin=unit(rep(0.1, 4), "lines")),
  gN.sp[[5]] + theme(legend.position="none", axis.text=element_blank(), 
  axis.line=element_blank(), axis.ticks=element_blank(), panel.border=element_blank(),
  panel.grid=element_blank(), plot.margin=unit(rep(0.1, 4), "lines")),
  gN.mean, main=" ", nrow=3, ncol=2))
dev.off()

save.image(file.name)
