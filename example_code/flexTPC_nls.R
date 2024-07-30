## Code by: Mauricio Cruz-Loya
## Stanford University

## Fits flexTPC model for thermal performance curves using
## nonlinear least squares. This example uses hatching
## rate data for the grapevine moth L. botrana.
library("nlstools")
library("nlraa")
library("nlme")

set.seed(42)

# FlexTPC model for thermal performance curves.
# (parametrized with Topt instead of alpha in order to easily get
# confidence intervals for Topt from nls function)
flexTPC_Topt <- function(T, Tmin, Tmax, rmax, Topt, B) {
  alpha <- (Topt - Tmin) / (Tmax - Tmin)
  beta <- B / (Tmax - Tmin)
  s <- alpha * (1 - alpha) / beta^2
  result <- rep(0, length(T))
  
  ## Return zeroes when Topt is outside the interval [Tmin, Tmax].
  #if((alpha < 0) | (alpha > 1)) {
  #  return(result)
  #}
  
  Tidx = (T > Tmin) & (T < Tmax)
  result[Tidx] <- rmax * exp(s * (alpha * log( (T[Tidx] - Tmin) / alpha) 
                                  + (1 - alpha) * log( (Tmax - T[Tidx]) / (1 - alpha))
                                  - log(Tmax - Tmin)) ) 
  return(result)
}


# Reads L. botrana dataset.
botrana <- read.csv("briere_data_L_botrana.csv")

##### Organize Data for JAGS
y <- 1 / botrana$eggs
temp <- botrana$T

fit.df <- data.frame(temp=temp, y=y)

## Note: flexTPC parametrized with Topt instead of alpha to get 95% confidence
## interval for Topt directly.

## It is recommended to try to pick good starting parameters for
## your data by using the preview function.
preview(y ~ flexTPC_Topt(temp, Tmin, Tmax, rmax, Topt, B), 
        start = list(Tmin=8, Tmax=32, Topt=30, rmax=0.2, B=5), data=fit.df)


eggs.flexTPC <- nls(y ~ flexTPC_Topt(temp, Tmin, Tmax, rmax, Topt, B), 
                    start = list(Tmin=8, Tmax=32, Topt=30, rmax=0.2, B=5),
                    lower= c(-10, 30, 5, 0, 0),
                    upper= c(15, 100, 40, 1, 20),
                    algorithm="port",
                    data=fit.df)

eggs.flexTPC
plotfit(eggs.flexTPC, smooth=TRUE)

## Two methods for calculating confidence intervals for parameters.

# Calculate asymptotic 95% confidence intervals for model parameters based on
# standard errors.
round(confint2(eggs.flexTPC), 3)

# Bootstrap confidence intervals (preferred method, but slower).
bts <- boot_nls(eggs.flexTPC, R=10000)
ci <- apply(bts$t, 2, quantile, c(0.025, 0.975), na.rm=TRUE)
colnames(ci) <- c("Tmin", "Tmax", "Topt", "rmax", "B")
t(round(ci, 3))