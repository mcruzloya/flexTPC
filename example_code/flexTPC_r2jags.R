## Code by: Mauricio Cruz-Loya
## Stanford University

## Fits flexTPC model for thermal performance curves with a Bayesian approach
## using JAGS (through the R2jags package).

library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')

set.seed(42)
setwd('/Users/cruzloya/git/flexTPC/example_code/')

# FlexTPC model for thermal performance curves.
flexTPC <- function(T, Tmin, Tmax, rmax, alpha, s) {
  beta <- 1 - alpha
  result <- rep(0, length(T))
  Tidx = (T > Tmin) & (T < Tmax)
  result[Tidx] <- rmax * exp(s * (alpha * log( (T[Tidx] - Tmin) / alpha) 
                                  + beta * log( (Tmax - T[Tidx]) / beta)
                                  - log(Tmax - Tmin)) ) 
  return(result)
}

## Defines JAGS model.
sink("flexTPC_norm.txt")
cat("
    model{
    
    ## Priors (please make sure to change these priors so they make sense for your system!)
    Tmin ~ dnorm(10, 1/5^2) # Minimum temperature (Celsius). Approx. prior 95% CI: [0, 20] 
    Tmax ~ dnorm(35, 1/5^2) # Maximum temperature (Celsius). Approx. prior 95% CI: [25, 45]
    rmax ~ dunif(0, 1)  # Maximum value of trait (i.e. trait value at Topt). Current prior assumes upper limit
                        # for maximum hatching rate of 1 /day since eggs will not hatch faster than in a day.
                        # Can be changed to e.g. a gamma prior if there is previous information for the trait
                        # from related species.
    alpha ~ dunif(0, 1) # Relative position of Topt relative to Tmin and Tmax.
                        # Uniform prior gives equal prior weight to curves of any skewness.
                        # Can also use e.g. beta prior to give more prior probability to left-skewed, symmetric
                        # or right-skewed curves depending on parameters of prior if it makes sense for the trait.
    logs ~ dnorm(0, 1 / 0.5^2) # Prior on base 10 logarithm of shape parameter. 
                               # Approx. prior logs 95% CI: [-1, 1] -> Approx s 95% CI: [0.1, 10] 
                               # This should work for most 'normal looking' TPCs, but can be relaxed
                               # to e.g. dnorm(0, 1) to allow TPCs that are very skinny or that are
                               # very flat in a wide temperature range.
    sigma ~ dunif(0, 1) # Standard deviation for the data points around the fitted TPC. 
                        # Upper limit of 1 for this problem, but please change this
                        # to make sense for the data you're fitting.
    
    # Derived quantities
    beta <- 1 - alpha
    Topt <- alpha * Tmax + beta * Tmin
    s <- 10^logs
    
    ## Likelihood
    for(i in 1:N.obs){
      # flexTPC model.
      mu[i] <- (Tmax > temp[i]) * (Tmin < temp[i]) * rmax * exp(s * (alpha * log( max(temp[i] - Tmin, 10^-20) / alpha) 
                            + beta * log( max(Tmax - temp[i], 10^-20) / beta)
                            - log(Tmax - Tmin)) ) 
      # Normal likelihood with constant variance.
      y[i] ~ dnorm(mu[i], 1 / sigma^2)
    }
    
    } # close model
    ",fill=TRUE)
sink()

##### Set MCMC Settings common to all traits
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 4 # number of chains

# Reads L. botrana dataset.
botrana <- read.csv("briere_data_L_botrana.csv")

##### Organize Data for JAGS
y <- 1 / botrana$eggs
temp <- botrana$T
N.obs <- length(y)

##### Bundle Data
jag.data<-list(y=y, temp = temp, N.obs=N.obs)

### Initial values for MCMC chains (chosen randomly from specified 
### distributions). Please remember to change these to something suitable
### for the data you're fitting.

inits<-function(){list(
  Tmin = runif(1, min=0, max=10),
  Tmax = runif(1, min=30, max=40),
  rmax = runif(1, min=0, max=0.5),
  alpha = runif(1, min=0, max=1),
  logs = rnorm(1, 0, 0.5),
  sigma = runif(1, min=0, max=0.2))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "logs", "s", "Topt", "sigma")

jags.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                          model.file="flexTPC_norm.txt", n.thin=nt, n.chains=nc, 
                          n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())

# Show summary, including posterior mean and 95% credible intervals for model parameters.
jags.out ## Check that n.eff > 10000 and Rhat < 1.01 for all parameters.
         ## Otherwise it may mean there's an issue with model (priors or likelihood),
         ## or that we need to increase the number of iterations.

# MCMC diagnostic plots (check to assess if it looks like the chains converged).
mcmcplot(jags.out)

# Plot data.
plot(temp, y, pch=20, xlab="Temperature [ºC]", ylab='Developmental rate [1 /day]', xlim=c(5, 35), ylim=c(0, 0.3),
     main="eggs")

# Extract MCMC chains for flexTPC parameters.
chains <- MCMCchains(jags.out, params=c("Tmin", "Tmax", "rmax", "alpha", "s"))

# Calculate TPCs for each posterior sample.
temps <- seq(0, 40, 0.1)
curves <- apply(chains, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(curves, 1, mean)
CI  <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="steelblue", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("steelblue", 0.2), lty=0)