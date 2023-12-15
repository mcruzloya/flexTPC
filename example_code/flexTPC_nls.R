## Code by: Mauricio Cruz-Loya
## Stanford University

## Fits flexTPC model for thermal performance curves using
## nonlinear least squares. This example uses hatching
## rate data for the grapevine moth L. botrana.
library("nlstools")

set.seed(42)
setwd('/Users/cruzloya/git/flexTPC/example_code/')

# FlexTPC model for thermal performance curves.
# (parametrized with Topt instead of alpha in order to easily get
# confidence intervals for Topt from nls function)
flexTPC_Topt <- function(T, Tmin, Tmax, rmax, Topt, s) {
  alpha <- (Topt - Tmin) / (Tmax - Tmin)
  beta <- 1 - alpha
  result <- rep(0, length(T))
  
  ## Return zeroes when Topt is outside the interval [Tmin, Tmax].
  if((alpha < 0) | (alpha > 1)) {
    return(result)
  }
  
  Tidx = (T > Tmin) & (T < Tmax)
  result[Tidx] <- rmax * exp(s * (alpha * log( (T[Tidx] - Tmin) / alpha) 
                                  + beta * log( (Tmax - T[Tidx]) / beta)
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
## interval from nls.

## It is recommended to try to pick good starting parameters for
## your data by using the preview function.
preview(y ~ flexTPC_Topt(temp, Tmin, Tmax, rmax, Topt, s), 
        start = list(Tmin=8, Tmax=32, Topt=30, rmax=0.2, s=2), data=fit.df)


eggs.flexTPC <- nls(y ~ flexTPC_Topt(temp, Tmin, Tmax, rmax, Topt, s), 
                    start = list(Tmin=8, Tmax=32, Topt=30, rmax=0.2, s=2),
                    lower= c(-10, 30, 5, 0, 0),
                    upper= c(15, 100, 40, 1, 10),
                    algorithm="port",
                    data=fit.df)

eggs.flexTPC
plotfit(eggs.flexTPC, smooth=TRUE)

# Calculate approximate 95% confidence intervals for model parameters.
round(confint2(eggs.flexTPC), 2)
