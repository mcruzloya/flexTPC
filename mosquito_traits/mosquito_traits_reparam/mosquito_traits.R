library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')

set.seed(42)
setwd('/Users/cruzloya/git/flexTPC/mosquito_traits_reparam')

lit.col = "purple"
flex.col = "darkgreen"

# FlexTPC model for thermal responses.
flexTPC <- function(T, Tmin, Tmax, rmax, alpha, beta) {
  result <- rep(0, length(T))
  Tidx = (T > Tmin) & (T < Tmax)
  s = alpha * (1 - alpha) / beta^2
  result[Tidx] <- rmax * exp(s * (alpha * log( (T[Tidx] - Tmin) / alpha) 
                                  +  (1 - alpha) * log( (Tmax - T[Tidx]) / (1 - alpha))
                                  - log(Tmax - Tmin)) ) 
  return(result)
}

# Briere1 TPC model.
briere <- function(T, Tmin, Tmax, c) {
  result <- c * T * (T - Tmin) * sqrt((Tmax - T) * (T > Tmin) * (T < Tmax))
  return(result)
}

# Quadratic TPC model.
quad <- function(T, Tmin, Tmax, c) {
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(result)
}

# Quadratic TPC model with upper limit at 1 for proportions.
quad_lim <- function(T, Tmin, Tmax, c) {
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(pmin(result, 1))
}

# Linear TPC model with upper limit.
linear_lim <- function(T, m, Tmax, Tlim) {
  result <- m * (Tmax - T)
  result[T < Tlim] <- m * (Tmax - Tlim) # Set maximum value at Tlim
  result[T > Tmax] <- 0 # Set maximum value to zero
  return(result)
}

# Function to extract log-likelihood array from JAGS output.
get_loglik <- function(jags.out, n.iter=31250, n.chains=4, N=7) {
  llarray <- array(0, dim=c(n.iter, n.chains, N))
  for(chain in 1:n.chains) {
   llarray[,chain,] <- MCMCchains(jags.out, params="log_lik", chain_num=chain) 
  }
  return(llarray)
}

##### Set MCMC Settings common to all traits
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 4 # number of chains

## Egg viability
#data.EV <- read.csv("./TraitData_EV.csv")
#data.EV.Cpip <- subset(data.EV, data.EV$host.code == "Cpip")
#data.EV.Cqui <- subset(data.EV, data.EV$host.code == "Cqui")
data.EV.Cpip <- read.csv("EV_Cpip.csv")
data.EV.Cqui <- read.csv("EV_Cqui.csv")

sink("quad_EV_binom.txt")
cat("
    model{
    
    ## Priors
    Tmin ~ dnorm(10, 1/5^2) 
    Tmax ~ dnorm(35, 1/5^2)
    rmax ~ dunif(0, 1)

    ## Derived quantities
    c <- 4 * rmax / (Tmax - Tmin)^2
    Topt <- (Tmin + Tmax) / 2
    
    ## Likelihood
    for(i in 1:N.obs){
    p[i] <- min(c * (temp[i] - Tmin) * (Tmax - temp[i]) * (temp[i] < Tmax) * (Tmin < temp[i]),
    1)
    n[i] ~ dbin(max(p[i], 10^-10), N[i])
    
    }
    
    } # close model
    ",fill=TRUE)
sink()

sink("briere_EV_binom.txt")
cat("
    model{
    
    ## Priors
    Tmin ~ dnorm(10, 1/5^2) 
    Tmax ~ dnorm(35, 1/5^2)
    rmax ~ dunif(0, 1)

    ## Derived quantities
    Topt <- (4*Tmax + 3*Tmin)/10 + sqrt(4 * Tmax^2 + (9/4) * Tmin^2 - 4*Tmax*Tmin) / 5
    c <- rmax / (Topt * (Topt - Tmin) * sqrt(max(Tmax - Topt, 10^-20)))
    
    ## Likelihood
    for(i in 1:N.obs){
    p[i] <- min(c * temp[i] * (temp[i] - Tmin) * 
                     sqrt((Tmax - temp[i]) * (Tmax > temp[i])) * (Tmin < temp[i]), 1)
    n[i] ~ dbin(max(p[i], 10^-10), N[i])
    
    }
    
    } # close model
    ",fill=TRUE)
sink()

sink("flex_EV_binom.txt")
cat("
    model{
    
    ## Priors
    Tmin ~ dnorm(10, 1/5^2) 
    Tmax ~ dnorm(35, 1/5^2)
    rmax ~ dunif(0, 1)
    alpha ~ dunif(0, 1)
    beta ~ dgamma(0.2^2 / 0.4^2, 0.2 / 0.4^2)
    
    
    # Derived quantities
    s <- alpha * (1 - alpha) / beta^2
    Topt <- alpha * Tmax + beta * Tmin

    ## Likelihood
    for(i in 1:N.obs){
    p[i] <- (Tmax > temp[i]) * (Tmin < temp[i]) * rmax * exp(s * (alpha * log( max(temp[i] - Tmin, 10^-20) / alpha) 
                          + (1 - alpha) * log( max(Tmax - temp[i], 10^-20) / (1 - alpha))
                          - log(Tmax - Tmin)) ) 
    
    n[i] ~ dbin(max(p[i], 10^-10), N[i])
    
    }
    
    } # close model
    ",fill=TRUE)
sink()

##### Organize Data for JAGS
N.Cpip <- data.EV.Cpip$N
n.Cpip <- data.EV.Cpip$hatched
temp.Cpip <- data.EV.Cpip$temperature
N.obs.Cpip <- length(N.Cpip)

N.Cqui <- data.EV.Cqui$N
n.Cqui <- data.EV.Cqui$hatched
temp.Cqui <- data.EV.Cqui$temperature
N.obs.Cqui <- length(N.Cqui)

###### EV - Culex pipiens
##### Bundle Data
jag.data<-list(N=N.Cpip, n=n.Cpip, temp = temp.Cpip, N.obs=N.obs.Cpip)
jag.data

#### Quadratic model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=30, max=40),
  rmax = runif(1, min=0, max=1))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "Topt", "c")

##### Run JAGS

quad.EV.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="quad_EV_binom.txt", n.thin=nt, n.chains=nc, 
                    n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
quad.EV.Cpip.out
mcmcplot(quad.EV.Cpip.out)

#### flexTPC model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=1),
  alpha = runif(1, min=0.2, max=0.8),
  beta = runif(1, min=0.05, max=0.3))}


inits
##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "s", "Topt")



flex.EV.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="flex_EV_binom.txt", n.thin=nt, n.chains=nc, 
                         n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
flex.EV.Cpip.out
mcmcplot(flex.EV.Cpip.out)


###### EV - Culex quinquefasciatus
##### Bundle Data
jag.data<-list(N=N.Cqui, n=n.Cqui, temp = temp.Cqui, N.obs=N.obs.Cqui)
jag.data

#### Quadratic model
inits<-function(){list(
  Tmin = runif(1, min=5, max=15),
  Tmax = runif(1, min=30, max=40),
  rmax = runif(1, min=0, max=1))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "Topt", "c")

##### Run JAGS

briere.EV.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="briere_EV_binom.txt", n.thin=nt, n.chains=nc, 
                         n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
briere.EV.Cqui.out
mcmcplot(briere.EV.Cqui.out)

#### flexTPC model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=1),
  alpha = runif(1, min=0.2, max=0.8),
  beta = runif(1, min=0.1, max=0.5))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "s", "Topt")

flex.EV.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="flex_EV_binom.txt", n.thin=nt, n.chains=nc, 
                         n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
flex.EV.Cqui.out
mcmcplot(flex.EV.Cqui.out)

## Extract log-likelihood
#briere.EV.Cqui.ll <- get_loglik(briere.EV.Cqui.out, N=19)
#flex.EV.Cqui.ll <- get_loglik(flex.EV.Cqui.out, N=19)

#briere.EV.Cqui.loo <- loo(briere.EV.Cqui.ll, cores=4, r_eff=relative_eff(exp(briere.EV.Cqui.ll)))
#flex.EV.Cqui.loo <- loo(flex.EV.Cqui.ll, cores=4, r_eff=relative_eff(exp(flex.EV.Cqui.ll)))

#flex.EV.Cqui.loo

#loo_compare(list(briere=briere.EV.Cqui.loo, flexTPC=flex.EV.Cqui.loo))


par(mfcol=c(2, 4))
data.EV.Cpip$sderr <- sqrt(data.EV.Cpip$p.hatched * (1 - data.EV.Cpip$p.hatched) / data.EV.Cpip$N)

plot(data.EV.Cpip$temperature, data.EV.Cpip$p.hatched, pch=20, xlim=c(0, 45),
     ylim=c(0,1), xlab="Temperature [°C]", ylab="Proportion hatched", main="EV (Culex pipiens)")
arrows(data.EV.Cpip$temperature, data.EV.Cpip$p.hatched, y1=data.EV.Cpip$p.hatched+data.EV.Cpip$sderr,
       angle=90, length=0.06, lwd=2)
arrows(data.EV.Cpip$temperature, data.EV.Cpip$p.hatched, y1=data.EV.Cpip$p.hatched-data.EV.Cpip$sderr,
       angle=90, length=0.06, lwd=2)

temps <- seq(0, 45, 0.05) 

# Plot quadratic model curves
chains.EV.Cpip.quad <- MCMCchains(quad.EV.Cpip.out, params=c("Tmin", "Tmax", "c"))

curves <- apply(chains.EV.Cpip.quad, 1, 
                function(x) quad_lim(temps, x[1], x[2], x[3]))
meancurve.EV.Cpip.quad <- apply(curves, 1, mean)
CI.EV.Cpip.quad  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.EV.Cpip.quad, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.EV.Cpip.quad[1,], rev(CI.EV.Cpip.quad[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.EV.Cpip.flex <- MCMCchains(flex.EV.Cpip.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.EV.Cpip.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.EV.Cpip.flex <- apply(curves, 1, mean)
CI.EV.Cpip.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.EV.Cpip.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.EV.Cpip.flex[1,], rev(CI.EV.Cpip.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)

#### EV - Culex quinquefasciatus

data.EV.Cqui$sderr <- sqrt(data.EV.Cqui$p.hatched * (1 - data.EV.Cqui$p.hatched) / data.EV.Cqui$N)
plot(data.EV.Cqui$temperature, data.EV.Cqui$p.hatched, pch=20, xlim=c(0, 45),
     ylim=c(0,1), xlab="Temperature [°C]", ylab="Proportion hatched", main="EV (Culex quinquefasciatus)")
arrows(data.EV.Cqui$temperature, data.EV.Cqui$p.hatched, y1=data.EV.Cqui$p.hatched+data.EV.Cqui$sderr,
       angle=90, length=0.06, lwd=2)
arrows(data.EV.Cqui$temperature, data.EV.Cqui$p.hatched, y1=data.EV.Cqui$p.hatched-data.EV.Cqui$sderr,
       angle=90, length=0.06, lwd=2)

# Plot Briere model curves
chains.EV.Cqui.briere <- MCMCchains(briere.EV.Cqui.out, params=c("Tmin", "Tmax", "c"))

curves <- apply(chains.EV.Cqui.briere, 1, 
                function(x) briere(temps, x[1], x[2], x[3]))
meancurve.EV.Cqui.briere <- apply(curves, 1, mean)
CI.EV.Cqui.briere  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.EV.Cqui.briere, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.EV.Cqui.briere[1,], rev(CI.EV.Cqui.briere[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.EV.Cqui.flex <- MCMCchains(flex.EV.Cqui.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.EV.Cqui.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.EV.Cqui.flex <- apply(curves, 1, mean)
CI.EV.Cqui.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.EV.Cqui.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.EV.Cqui.flex[1,], rev(CI.EV.Cqui.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)




###### Probability of larval survival.
data.pLA <- read.csv("./TraitData_pLA.csv")
data.pLA.Cpip <- subset(data.pLA, data.pLA$host.code == "Cpip")
data.pLA.Cqui <- subset(data.pLA, data.pLA$host.code == "Cqui")
data.pLA.Cpip
data.pLA.Cqui

sink("quad_pLA_norm.txt")
cat("
    model{
    
    ## Priors
    Tmin ~ dnorm(5, 1/2.5^2) 
    Tmax ~ dnorm(35, 1/5^2)
    rmax ~ dunif(0, 1)
    sigma ~ dunif(0, 1)

    ## Derived quantities
    c <- 4 * rmax / (Tmax - Tmin)^2
    Topt <- (Tmin + Tmax) / 2
    
    ## Likelihood
    for(i in 1:N.obs){
    mu[i] <- min(c * (temp[i] - Tmin) * (Tmax - temp[i]) * (temp[i] < Tmax) * (Tmin < temp[i]),
    1)
    y[i] ~ dnorm(mu[i], 1 / sigma^2)
    }
    
    } # close model
    ",fill=TRUE)
sink()

sink("flex_pLA_norm.txt")
cat("
    model{
    
    ## Priors
    Tmin ~ dnorm(5, 1/2.5^2) 
    Tmax ~ dnorm(35, 1/5^2)
    rmax ~ dunif(0, 1)
    alpha ~ dunif(0, 1)
    beta ~ dgamma(0.2^2 / 0.4^2, 0.2 / 0.4^2)
    sigma ~ dunif(0, 1)
    
    # Derived quantities
    s <- alpha * (1 - alpha) / beta^2
    Topt <- alpha * Tmax + beta * Tmin

    ## Likelihood
    for(i in 1:N.obs){
    mu[i] <- (Tmax > temp[i]) * (Tmin < temp[i]) * rmax * exp(s * (alpha * log( max(temp[i] - Tmin, 10^-20) / alpha) 
                          + (1 - alpha) * log( max(Tmax - temp[i], 10^-20) / (1 - alpha))
                          - log(Tmax - Tmin)) ) 
    
    y[i] ~ dnorm(mu[i], 1 / sigma^2)
    }
    
    } # close model
    ",fill=TRUE)
sink()

##### Organize Data for JAGS
y.Cpip <- data.pLA.Cpip$trait
temp.Cpip <- data.pLA.Cpip$T
N.obs.Cpip <- length(y.Cpip)

y.Cqui <- data.pLA.Cqui$trait
temp.Cqui <- data.pLA.Cqui$T
N.obs.Cqui <- length(y.Cqui)

###### pLA - Culex pipiens
##### Bundle Data
jag.data<-list(y=y.Cpip, temp = temp.Cpip, N.obs=N.obs.Cpip)
jag.data

#### Quadratic model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=30, max=40),
  rmax = runif(1, min=0, max=1))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "c", "sigma")

##### Run JAGS
quad.pLA.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                        model.file="quad_pLA_norm.txt", n.thin=nt, n.chains=nc, 
                        n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
quad.pLA.Cpip.out
#mcmcplot(quad.pLA.Cpip.out)

#### flexTPC model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=1),
  alpha = runif(1, min=0.2, max=0.8),
  beta = runif(1, min=0.1, max=0.5))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "s", "Topt", "sigma")

jag.data<-list(y=y.Cpip, temp = temp.Cpip, N.obs=N.obs.Cpip)
flex.pLA.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="flex_pLA_norm.txt", n.thin=nt, n.chains=nc, 
                         n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
flex.pLA.Cpip.out
#mcmcplot(flex.pLA.Cpip.out)


##### Bundle Data
jag.data<-list(y=y.Cqui, temp = temp.Cqui, N.obs=N.obs.Cqui)
jag.data

#### Quadratic model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=30, max=40),
  rmax = runif(1, min=0, max=1))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "c", "sigma")

##### Run JAGS
quad.pLA.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                        model.file="quad_pLA_norm.txt", n.thin=nt, n.chains=nc, 
                        n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
quad.pLA.Cqui.out
#mcmcplot(quad.pLA.Cqui.out)

#### flexTPC model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=1),
  alpha = runif(1, min=0.2, max=0.8),
  beta = runif(1, min=0.1, max=0.5))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "s", "Topt", "sigma")

jag.data<-list(y=y.Cqui, temp = temp.Cqui, N.obs=N.obs.Cqui)

flex.pLA.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="flex_pLA_norm.txt", n.thin=nt, n.chains=nc, 
                         n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
flex.pLA.Cqui.out
mcmcplot(flex.pLA.Cqui.out)

### pLA - Cx pipiens plots

plot(data.pLA.Cpip$T, data.pLA.Cpip$trait, pch=20, xlim=c(0, 45),
     ylim=c(0,1), xlab="Temperature [°C]",
     ylab="Larval survival proportion",
     main="pLA (Culex pipiens)")

# Plot linear model curves
chains.pLA.Cpip.quad <- MCMCchains(quad.pLA.Cpip.out, params=c("Tmin", "Tmax", "c"))

curves <- apply(chains.pLA.Cpip.quad, 1, 
                function(x) quad_lim(temps, x[1], x[2], x[3]))
meancurve.pLA.Cpip.quad <- apply(curves, 1, mean)
CI.pLA.Cpip.quad  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.pLA.Cpip.quad, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.pLA.Cpip.quad[1,], rev(CI.pLA.Cpip.quad[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.pLA.Cpip.flex <- MCMCchains(flex.pLA.Cpip.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.pLA.Cpip.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.pLA.Cpip.flex <- apply(curves, 1, mean)
CI.pLA.Cpip.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.pLA.Cpip.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.pLA.Cpip.flex[1,], rev(CI.pLA.Cpip.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)

#### pLA - Culex quinquefasciatus

plot(data.pLA.Cqui$T, data.pLA.Cqui$trait, pch=20, xlim=c(0, 45),
     ylim=c(0, 1), xlab="Temperature [°C]", ylab="Larval survival proportion", 
     main="pLA (Culex quinquefasciatus)")


# Plot quadratic model curves
chains.pLA.Cqui.quad <- MCMCchains(quad.pLA.Cqui.out, params=c("Tmin", "Tmax", "c"))

curves <- apply(chains.pLA.Cqui.quad, 1, 
                function(x) quad_lim(temps, x[1], x[2], x[3]))
meancurve.pLA.Cqui.quad <- apply(curves, 1, mean)
CI.pLA.Cqui.quad  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.pLA.Cqui.quad, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.pLA.Cqui.quad[1,], rev(CI.pLA.Cqui.quad[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.pLA.Cqui.flex <- MCMCchains(flex.pLA.Cqui.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.pLA.Cqui.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.pLA.Cqui.flex <- apply(curves, 1, mean)
CI.pLA.Cqui.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.pLA.Cqui.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.pLA.Cqui.flex[1,], rev(CI.pLA.Cqui.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)


####### Mosquito Development Rate
data.MDR <- read.csv("./TraitData_MDR.csv")
data.MDR.Cpip <- subset(data.MDR, data.MDR$host.code == "Cpip")
data.MDR.Cqui <- subset(data.MDR, data.MDR$host.code == "Cqui")
data.MDR.Cpip
data.MDR.Cqui

sink("briere_MDR_norm.txt")
cat("
    model{
    
    ## Priors
    Tmin ~ dnorm(5, 1/2.5^2) 
    Tmax ~ dnorm(35, 1/5^2)
    rmax ~ dunif(0, 1)
    sigma ~ dunif(0, 1)

    ## Derived quantities
    Topt <- (4*Tmax + 3*Tmin)/10 + sqrt(4 * Tmax^2 + (9/4) * Tmin^2 - 4*Tmax*Tmin) / 5
    c <- rmax / (Topt * (Topt - Tmin) * sqrt(max(Tmax - Topt, 10^-20)))
    
    ## Likelihood
    for(i in 1:N.obs){
    mu[i] <- c * temp[i] * (temp[i] - Tmin) * sqrt((Tmax - temp[i]) * (Tmax > temp[i])) * (Tmin < temp[i])
    y[i] ~ dnorm(mu[i], 1 / sigma^2)
    }
    
    } # close model
    ",fill=TRUE)
sink()

sink("flex_MDR_norm.txt")
cat("
    model{
    
    ## Priors
    Tmin ~ dnorm(5, 1/2.5^2) 
    Tmax ~ dnorm(35, 1/5^2)
    rmax ~ dunif(0, 1)
    alpha ~ dunif(0, 1)
    
    beta ~ dgamma(0.2^2 / 0.4^2, 0.2 / 0.4^2)
    sigma ~ dunif(0, 1)
    
    # Derived quantities
    s <- alpha * (1 - alpha) / beta^2
    Topt <- alpha * Tmax + beta * Tmin

    ## Likelihood
    for(i in 1:N.obs){
    mu[i] <- (Tmax > temp[i]) * (Tmin < temp[i]) * rmax * exp(s * (alpha * log( max(temp[i] - Tmin, 10^-20) / alpha) 
                          + (1 - alpha) * log( max(Tmax - temp[i], 10^-20) / (1 - alpha))
                          - log(Tmax - Tmin)) ) 
    
    y[i] ~ dnorm(mu[i], 1 / sigma^2)
    }
    
    } # close model
    ",fill=TRUE)
sink()

##### Organize Data for JAGS
y.Cpip <- 1 / data.MDR.Cpip$trait
temp.Cpip <- data.MDR.Cpip$T
N.obs.Cpip <- length(y.Cpip)

y.Cqui <- 1 / data.MDR.Cqui$trait
temp.Cqui <- data.MDR.Cqui$T
N.obs.Cqui <- length(y.Cqui)

###### pLA - Culex pipiens
##### Bundle Data
jag.data<-list(y=y.Cpip, temp = temp.Cpip, N.obs=N.obs.Cpip)
jag.data

#### Quadratic model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=30, max=40),
  rmax = runif(1, min=0, max=1))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "c", "sigma")

##### Run JAGS
briere.MDR.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                          model.file="briere_MDR_norm.txt", n.thin=nt, n.chains=nc, 
                          n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
briere.MDR.Cpip.out
#mcmcplot(briere.MDR.Cpip.out)

#### flexTPC model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=1),
  alpha = runif(1, min=0.2, max=0.8),
  beta = runif(1, min=0.1, max=0.5))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "s", "Topt", "sigma")

jag.data<-list(y=y.Cpip, temp = temp.Cpip, N.obs=N.obs.Cpip)
flex.MDR.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                          model.file="flex_MDR_norm.txt", n.thin=nt, n.chains=nc, 
                          n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
flex.MDR.Cpip.out
#mcmcplot(flex.MDR.Cpip.out)


##### Bundle Data
jag.data<-list(y=y.Cqui, temp = temp.Cqui, N.obs=N.obs.Cqui)
jag.data

#### Quadratic model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=30, max=40),
  rmax = runif(1, min=0, max=1))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "c", "sigma")

##### Run JAGS
briere.MDR.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                          model.file="briere_MDR_norm.txt", n.thin=nt, n.chains=nc, 
                          n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
briere.MDR.Cqui.out
#mcmcplot(briere.MDR.Cqui.out)

#### flexTPC model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=1),
  alpha = runif(1, min=0.2, max=0.8),
  beta = runif(1, min=0.1, max=0.5))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "s", "Topt", "sigma")

jag.data<-list(y=y.Cqui, temp = temp.Cqui, N.obs=N.obs.Cqui)

flex.MDR.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                          model.file="flex_MDR_norm.txt", n.thin=nt, n.chains=nc, 
                          n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
flex.MDR.Cqui.out
mcmcplot(flex.MDR.Cqui.out)

##### lf - Plot curves

plot(data.MDR.Cpip$T, 1 / data.MDR.Cpip$trait, pch=20, xlim=c(0, 45),
     ylim=c(0, 0.3), xlab="Temperature [°C]", ylab="Mosquito development rate [1 / day]", 
     main="MDR (Culex pipiens)")

temps <- seq(0, 45, 0.05)

## Cx pipiens

# Plot Briere1 model curves
chains.MDR.Cpip.briere <- MCMCchains(briere.MDR.Cpip.out, params=c("Tmin", "Tmax", "c"))

curves <- apply(chains.MDR.Cpip.briere, 1, 
                function(x) briere(temps, x[1], x[2], x[3]))
meancurve.MDR.Cpip.briere <- apply(curves, 1, mean)
CI.MDR.Cpip.briere  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.MDR.Cpip.briere, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.MDR.Cpip.briere[1,],
                                rev(CI.MDR.Cpip.briere[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.MDR.Cpip.flex <- MCMCchains(flex.MDR.Cpip.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.MDR.Cpip.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.MDR.Cpip.flex <- apply(curves, 1, mean)
CI.MDR.Cpip.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.MDR.Cpip.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.MDR.Cpip.flex[1,], 
                                rev(CI.MDR.Cpip.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)

#### MDR - Culex quinquefasciatus plots

plot(data.MDR.Cqui$T, 1 / data.MDR.Cqui$trait, pch=20, xlim=c(0, 45),
     ylim=c(0, 0.3), xlab="Temperature [°C]", ylab="Mosquito development rate [1 / day]", 
     main="MDR (Culex quinquefasciatus)")


# Plot linear model curves
chains.MDR.Cqui.briere <- MCMCchains(briere.MDR.Cqui.out, params=c("Tmin", "Tmax", "c"))

curves <- apply(chains.MDR.Cqui.briere, 1, 
                function(x) briere(temps, x[1], x[2], x[3]))
meancurve.MDR.Cqui.briere <- apply(curves, 1, mean)
CI.MDR.Cqui.briere  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.MDR.Cqui.briere, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.MDR.Cqui.briere[1,], rev(CI.MDR.Cqui.briere[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.MDR.Cqui.flex <- MCMCchains(flex.MDR.Cqui.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.MDR.Cqui.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.MDR.Cqui.flex <- apply(curves, 1, mean)
CI.MDR.Cqui.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.MDR.Cqui.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.MDR.Cqui.flex[1,], rev(CI.MDR.Cqui.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)

####### Adult lifespan
data.lf <- read.csv("./TraitData_lf.csv")

data.lf.Cpip.fem <- subset(data.lf, host.code == "Cpip" & Trait != "Male") # For some reason this is also excluding NAs, so work around below
data.lf.Cpip <- subset(data.lf, host.code == "Cpip") # For some reason this is also excluding NAs, so work around below

which(data.lf.Cpip$Trait == "Male")
nrow(data.lf.Cpip)
data.lf.Cpip.fem <- rbind(data.lf.Cpip[1:15,], data.lf.Cpip[21:37,])
data.lf.Cqui.fem <- subset(data.lf, host.code == "Cqui" & Trait != "Male")

sink("linear_lf_norm.txt")
cat("
    model{
    
    ## Priors
    Tmax ~ dnorm(35, 1/5^2)
    rmax ~ dunif(0, 150)
    sigma ~ dunif(0, 100)

    ## Derived quantities
    m <- rmax / (Tmax - Tlim)
    
    ## Likelihood
    for(i in 1:N.obs){
    mu[i] <- m * (Tmax - temp[i]) 
    y[i] ~ dnorm(mu[i], 1 / sigma^2)
    }
    
    } # close model
    ",fill=TRUE)
sink()

sink("flex_lf_norm.txt")
cat("
    model{
    
    ## Priors
    Tmin ~ dnorm(5, 1/2.5^2) 
    Tmax ~ dnorm(35, 1/5^2)
    rmax ~ dunif(0, 150)
    alpha ~ dunif(0, 1)
    beta ~ dgamma(0.2^2 / 0.4^2, 0.2 / 0.4^2)
    sigma ~ dunif(0, 100)
    
    # Derived quantities
    s <- alpha * (1 - alpha) / beta^2
    Topt <- alpha * Tmax + beta * Tmin

    ## Likelihood
    for(i in 1:N.obs){
    mu[i] <- (Tmax > temp[i]) * (Tmin < temp[i]) * rmax * exp(s * (alpha * log( max(temp[i] - Tmin, 10^-20) / alpha) 
                          + (1 - alpha) * log( max(Tmax - temp[i], 10^-20) / (1 - alpha))
                          - log(Tmax - Tmin)) ) 
    
    y[i] ~ dnorm(mu[i], 1 / sigma^2)
    }
    
    } # close model
    ",fill=TRUE)
sink()

##### Organize Data for JAGS
y.Cpip <- data.lf.Cpip.fem$trait
temp.Cpip <- data.lf.Cpip.fem$T
N.obs.Cpip <- length(y.Cpip)

y.Cqui <- data.lf.Cqui.fem$trait
temp.Cqui <- data.lf.Cqui.fem$T
N.obs.Cqui <- length(y.Cqui)

###### Lf - Culex pipiens
##### Bundle Data
jag.data<-list(y=y.Cpip, temp = temp.Cpip, N.obs=N.obs.Cpip, Tlim=15)
jag.data

#### Linear model
inits<-function(){list(
  Tmax = runif(1, min=30, max=40),
  rmax = runif(1, min=0, max=150))}

##### Parameters to Estimate
parameters <- c("Tmax", "rmax", "m", "sigma")

##### Run JAGS
lin.lf.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="linear_lf_norm.txt", n.thin=nt, n.chains=nc, 
                         n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
lin.lf.Cpip.out
#mcmcplot(lin.lf.Cpip.out)

#### flexTPC model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=150),
  alpha = runif(1, min=0.2, max=0.8),
  beta = runif(1, min=0.1, max=0.5))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "s", "Topt", "sigma")

jag.data<-list(y=y.Cpip, temp = temp.Cpip, N.obs=N.obs.Cpip)
flex.lf.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="flex_lf_norm.txt", n.thin=nt, n.chains=nc, 
                         n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
flex.lf.Cpip.out
#mcmcplot(flex.lf.Cpip.out)

### Cqui lifespan

##### Bundle Data
jag.data<-list(y=y.Cqui, temp = temp.Cqui, N.obs=N.obs.Cqui, Tlim=16)
jag.data

#### Linear model
inits<-function(){list(
  Tmax = runif(1, min=30, max=40),
  rmax = runif(1, min=0, max=150))}

##### Parameters to Estimate
parameters <- c("Tmax", "rmax", "m", "sigma")

##### Run JAGS
lin.lf.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                        model.file="linear_lf_norm.txt", n.thin=nt, n.chains=nc, 
                        n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
lin.lf.Cqui.out
#mcmcplot(lin.lf.Cqui.out)

#### flexTPC model
inits<-function(){list(
  Tmin = runif(1, min=2.5, max=7.5),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=150),
  alpha = runif(1, min=0.2, max=0.8),
  beta = runif(1, min=0.1, max=0.5))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "s", "Topt", "sigma")

jag.data<-list(y=y.Cqui, temp = temp.Cqui, N.obs=N.obs.Cqui)
flex.lf.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="flex_lf_norm.txt", n.thin=nt, n.chains=nc, 
                         n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
flex.lf.Cqui.out
#mcmcplot(flex.lf.Cqui.out)


##### lf - Plot curves

plot(data.lf.Cpip.fem$T, data.lf.Cpip.fem$trait, pch=20, xlim=c(0, 45),
     ylim=c(0, 150), xlab="Temperature [°C]", ylab="Female adult lifespan [days]", 
     main="lf (Culex pipiens)")

temps <- seq(0, 45, 0.05)

## Cx pipiens

# Plot linear model curves
chains.lf.Cpip.lin <- MCMCchains(lin.lf.Cpip.out, params=c("m", "Tmax"))

curves <- apply(chains.lf.Cpip.lin, 1, 
                function(x) linear_lim(temps, x[1], x[2], 15))
meancurve.lf.Cpip.lin <- apply(curves, 1, mean)
CI.lf.Cpip.lin  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.lf.Cpip.lin, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.lf.Cpip.lin[1,], rev(CI.lf.Cpip.lin[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.lf.Cpip.flex <- MCMCchains(flex.lf.Cpip.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.lf.Cpip.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.lf.Cpip.flex <- apply(curves, 1, mean)
CI.lf.Cpip.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.lf.Cpip.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.lf.Cpip.flex[1,], rev(CI.lf.Cpip.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)

#### lf - Culex quinquefasciatus

plot(data.lf.Cqui.fem$T, data.lf.Cqui.fem$trait, pch=20, xlim=c(0, 45),
     ylim=c(0, 150), xlab="Temperature [°C]", ylab="Female adult lifespan [days]", 
     main="lf (Culex quinquefasciatus)")


# Plot linear model curves
chains.lf.Cqui.lin <- MCMCchains(lin.lf.Cqui.out, params=c("m", "Tmax"))

curves <- apply(chains.lf.Cqui.lin, 1, 
                function(x) linear_lim(temps, x[1], x[2], 15))
meancurve.lf.Cqui.lin <- apply(curves, 1, mean)
CI.lf.Cqui.lin  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.lf.Cqui.lin, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.lf.Cqui.lin[1,], rev(CI.lf.Cqui.lin[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.lf.Cqui.flex <- MCMCchains(flex.lf.Cqui.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.lf.Cqui.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.lf.Cqui.flex <- apply(curves, 1, mean)
CI.lf.Cqui.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.lf.Cqui.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.lf.Cqui.flex[1,], rev(CI.lf.Cqui.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)




############## Final figure for paper

par(mfcol=c(2, 5))

### Common plotting parameters
cex.main = 2.0
cex.lab = 1.8
cex.axis = 1.5

# Left margin offset
lof <- 0.5


par(mar=(c(2, 4, 4, 2) + 0.1))
plot("", axes=FALSE, frame.plot=FALSE)
text(0.5, 0.57, "Culex", cex=1.8)
text(0.5, 0.43, "pipiens", cex=1.8)

par(mar=(c(5, 4, 2, 2) + 0.1))
plot("", axes=FALSE, frame.plot=FALSE)
text(0.5, 0.57, "Culex", cex=1.8)
text(0.5, 0.43, "quinquefasciatus", cex=1.8)


############ EV 
data.EV.Cpip$sderr <- sqrt(data.EV.Cpip$p.hatched * (1 - data.EV.Cpip$p.hatched) / data.EV.Cpip$N)

par(mar=(c(2, 4 + lof, 4, 2) + 0.1))

plot(data.EV.Cpip$temperature, data.EV.Cpip$p.hatched, pch=20, xlim=c(0, 45),
     ylim=c(0,1), xlab="", ylab="Proportion hatched", main="EV",
     axes=FALSE,
     frame.plot=TRUE,
     cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)
Axis(side=1, labels=FALSE, cex=cex.axis)
Axis(side=2, cex.axis=cex.axis)
arrows(data.EV.Cpip$temperature, data.EV.Cpip$p.hatched, y1=data.EV.Cpip$p.hatched+data.EV.Cpip$sderr,
       angle=90, length=0.06, lwd=2)
arrows(data.EV.Cpip$temperature, data.EV.Cpip$p.hatched, y1=data.EV.Cpip$p.hatched-data.EV.Cpip$sderr,
       angle=90, length=0.06, lwd=2)

temps <- seq(0, 45, 0.05) 

# Plot quadratic model curves
chains.EV.Cpip.quad <- MCMCchains(quad.EV.Cpip.out, params=c("Tmin", "Tmax", "c"))

curves <- apply(chains.EV.Cpip.quad, 1, 
                function(x) quad_lim(temps, x[1], x[2], x[3]))
meancurve.EV.Cpip.quad <- apply(curves, 1, mean)
CI.EV.Cpip.quad  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.EV.Cpip.quad, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.EV.Cpip.quad[1,], rev(CI.EV.Cpip.quad[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.EV.Cpip.flex <- MCMCchains(flex.EV.Cpip.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.EV.Cpip.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.EV.Cpip.flex <- apply(curves, 1, mean)
CI.EV.Cpip.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.EV.Cpip.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.EV.Cpip.flex[1,], rev(CI.EV.Cpip.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)

#### EV - Culex quinquefasciatus

data.EV.Cqui$sderr <- sqrt(data.EV.Cqui$p.hatched * (1 - data.EV.Cqui$p.hatched) / data.EV.Cqui$N)

par(mar=(c(5, 4 + lof, 2, 2) + 0.1))

plot(data.EV.Cqui$temperature, data.EV.Cqui$p.hatched, pch=20, xlim=c(0, 45),
     ylim=c(0,1), xlab="Temperature [°C]", ylab="Proportion hatched", main="",
     cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)
arrows(data.EV.Cqui$temperature, data.EV.Cqui$p.hatched, y1=data.EV.Cqui$p.hatched+data.EV.Cqui$sderr,
       angle=90, length=0.06, lwd=2)
arrows(data.EV.Cqui$temperature, data.EV.Cqui$p.hatched, y1=data.EV.Cqui$p.hatched-data.EV.Cqui$sderr,
       angle=90, length=0.06, lwd=2)

# Plot Briere model curves
chains.EV.Cqui.briere <- MCMCchains(briere.EV.Cqui.out, params=c("Tmin", "Tmax", "c"))

curves <- apply(chains.EV.Cqui.briere, 1, 
                function(x) briere(temps, x[1], x[2], x[3]))
meancurve.EV.Cqui.briere <- apply(curves, 1, mean)
CI.EV.Cqui.briere  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.EV.Cqui.briere, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.EV.Cqui.briere[1,], rev(CI.EV.Cqui.briere[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.EV.Cqui.flex <- MCMCchains(flex.EV.Cqui.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.EV.Cqui.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.EV.Cqui.flex <- apply(curves, 1, mean)
CI.EV.Cqui.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.EV.Cqui.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.EV.Cqui.flex[1,], rev(CI.EV.Cqui.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)



### pLA - Cx pipiens plots

par(mar=(c(2, 4 + lof, 4, 2) + 0.1))
plot(data.pLA.Cpip$T, data.pLA.Cpip$trait, pch=20, xlim=c(0, 45),
     ylim=c(0,1), xlab="",
     ylab="Larval survival proportion",
     main="pLA",
     axes=FALSE,
     frame.plot=TRUE,
     cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)
Axis(side=1, labels=FALSE, cex=cex.axis)
Axis(side=2, cex.axis=cex.axis)

# Plot linear model curves
chains.pLA.Cpip.quad <- MCMCchains(quad.pLA.Cpip.out, params=c("Tmin", "Tmax", "c"))

curves <- apply(chains.pLA.Cpip.quad, 1, 
                function(x) quad_lim(temps, x[1], x[2], x[3]))
meancurve.pLA.Cpip.quad <- apply(curves, 1, mean)
CI.pLA.Cpip.quad  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.pLA.Cpip.quad, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.pLA.Cpip.quad[1,], rev(CI.pLA.Cpip.quad[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.pLA.Cpip.flex <- MCMCchains(flex.pLA.Cpip.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.pLA.Cpip.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.pLA.Cpip.flex <- apply(curves, 1, mean)
CI.pLA.Cpip.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.pLA.Cpip.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.pLA.Cpip.flex[1,], rev(CI.pLA.Cpip.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)

#### pLA - Culex quinquefasciatus
par(mar=(c(5, 4 + lof, 2, 2) + 0.1))
plot(data.pLA.Cqui$T, data.pLA.Cqui$trait, pch=20, xlim=c(0, 45),
     ylim=c(0, 1), xlab="Temperature [°C]", ylab="Larval survival proportion", 
     main="", cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)

# Plot quadratic model curves
chains.pLA.Cqui.quad <- MCMCchains(quad.pLA.Cqui.out, params=c("Tmin", "Tmax", "c"))

curves <- apply(chains.pLA.Cqui.quad, 1, 
                function(x) quad_lim(temps, x[1], x[2], x[3]))
meancurve.pLA.Cqui.quad <- apply(curves, 1, mean)
CI.pLA.Cqui.quad  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.pLA.Cqui.quad, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.pLA.Cqui.quad[1,], rev(CI.pLA.Cqui.quad[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.pLA.Cqui.flex <- MCMCchains(flex.pLA.Cqui.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.pLA.Cqui.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.pLA.Cqui.flex <- apply(curves, 1, mean)
CI.pLA.Cqui.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.pLA.Cqui.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.pLA.Cqui.flex[1,], rev(CI.pLA.Cqui.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)


##### MDR - Plot curves

par(mar=(c(2, 4 + lof, 4, 2) + 0.1))
plot(data.MDR.Cpip$T, 1 / data.MDR.Cpip$trait, pch=20, xlim=c(0, 45),
     ylim=c(0, 0.2), xlab="Temperature", ylab="Development rate [1 / day]", 
     main="MDR",
     axes=FALSE,
     frame.plot=TRUE,
     cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)
Axis(side=1, labels=FALSE, cex=cex.axis)
Axis(side=2, cex.axis=cex.axis)

temps <- seq(0, 45, 0.05)

## Cx pipiens

# Plot Briere1 model curves
chains.MDR.Cpip.briere <- MCMCchains(briere.MDR.Cpip.out, params=c("Tmin", "Tmax", "c"))

curves <- apply(chains.MDR.Cpip.briere, 1, 
                function(x) briere(temps, x[1], x[2], x[3]))
meancurve.MDR.Cpip.briere <- apply(curves, 1, mean)
CI.MDR.Cpip.briere  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.MDR.Cpip.briere, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.MDR.Cpip.briere[1,],
                                rev(CI.MDR.Cpip.briere[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.MDR.Cpip.flex <- MCMCchains(flex.MDR.Cpip.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.MDR.Cpip.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.MDR.Cpip.flex <- apply(curves, 1, mean)
CI.MDR.Cpip.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.MDR.Cpip.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.MDR.Cpip.flex[1,], 
                                rev(CI.MDR.Cpip.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)

#### MDR - Culex quinquefasciatus plots
par(mar=(c(5, 4 + lof, 2, 2) + 0.1))
plot(data.MDR.Cqui$T, 1 / data.MDR.Cqui$trait, pch=20, xlim=c(0, 45),
     ylim=c(0, 0.2), xlab="Temperature [°C]", ylab="Development rate [1 / day]", 
     main="", cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)


# Plot linear model curves
chains.MDR.Cqui.briere <- MCMCchains(briere.MDR.Cqui.out, params=c("Tmin", "Tmax", "c"))

curves <- apply(chains.MDR.Cqui.briere, 1, 
                function(x) briere(temps, x[1], x[2], x[3]))
meancurve.MDR.Cqui.briere <- apply(curves, 1, mean)
CI.MDR.Cqui.briere  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.MDR.Cqui.briere, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.MDR.Cqui.briere[1,], rev(CI.MDR.Cqui.briere[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.MDR.Cqui.flex <- MCMCchains(flex.MDR.Cqui.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.MDR.Cqui.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.MDR.Cqui.flex <- apply(curves, 1, mean)
CI.MDR.Cqui.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.MDR.Cqui.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.MDR.Cqui.flex[1,], rev(CI.MDR.Cqui.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)


##### lf - Plot curves
par(mar=(c(2, 4 + lof, 4, 2) + 0.1))
plot(data.lf.Cpip.fem$T, data.lf.Cpip.fem$trait, pch=20, xlim=c(0, 45),
     ylim=c(0, 150), xlab="", ylab="Adult lifespan [days]", 
     main="lf",
     axes=FALSE,
     frame.plot=TRUE,
     cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)
Axis(side=1, labels=FALSE, cex=cex.axis)
Axis(side=2, cex.axis=cex.axis)

temps <- seq(0, 45, 0.05)

## Cx pipiens

# Plot linear model curves
chains.lf.Cpip.lin <- MCMCchains(lin.lf.Cpip.out, params=c("m", "Tmax"))

curves <- apply(chains.lf.Cpip.lin, 1, 
                function(x) linear_lim(temps, x[1], x[2], 15))
meancurve.lf.Cpip.lin <- apply(curves, 1, mean)
CI.lf.Cpip.lin  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.lf.Cpip.lin, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.lf.Cpip.lin[1,], rev(CI.lf.Cpip.lin[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.lf.Cpip.flex <- MCMCchains(flex.lf.Cpip.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.lf.Cpip.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.lf.Cpip.flex <- apply(curves, 1, mean)
CI.lf.Cpip.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.lf.Cpip.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.lf.Cpip.flex[1,], rev(CI.lf.Cpip.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)

#### lf - Culex quinquefasciatus
par(mar=(c(5, 4 + lof, 2, 2) + 0.1))
plot(data.lf.Cqui.fem$T, data.lf.Cqui.fem$trait, pch=20, xlim=c(0, 45),
     ylim=c(0, 150), xlab="Temperature [°C]", ylab="Adult lifespan [days]", 
     main="", cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)


# Plot linear model curves
chains.lf.Cqui.lin <- MCMCchains(lin.lf.Cqui.out, params=c("m", "Tmax"))

curves <- apply(chains.lf.Cqui.lin, 1, 
                function(x) linear_lim(temps, x[1], x[2], 15))
meancurve.lf.Cqui.lin <- apply(curves, 1, mean)
CI.lf.Cqui.lin  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.lf.Cqui.lin, col=lit.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.lf.Cqui.lin[1,], rev(CI.lf.Cqui.lin[2,])), 
        col=alpha(lit.col, 0.2), lty=0)

# Plot flexTPC model curves
chains.lf.Cqui.flex <- MCMCchains(flex.lf.Cqui.out, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
curves <- apply(chains.lf.Cqui.flex, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
meancurve.lf.Cqui.flex <- apply(curves, 1, mean)
CI.lf.Cqui.flex  <- apply(curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.lf.Cqui.flex, col=flex.col, lwd=1.5)
polygon(c(temps, rev(temps)), c(CI.lf.Cqui.flex[1,], rev(CI.lf.Cqui.flex[2,])), 
        col=alpha(flex.col, 0.2), lty=0)

