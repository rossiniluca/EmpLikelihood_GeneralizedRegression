############################################################################################################################################################################################
# Main file for the simulation experiments for the generalized logistic regression and for standard one
# by using the empirical likelihood based equation.
# PLEASE CITE THE PAPER:
# Dalla Valle, L., Leisen, F., Rossini, L. and Zhu, W. (2019) - 
# "Bayesian Analysis of immigration in Europe with generalized logistic regression" - Journal of Applied Statistics
# https://www.tandfonline.com/doi/full/10.1080/02664763.2019.1642310
############################################################################################################################################################################################

rm(list=ls(all=TRUE))

# Load the necessary library
library(emplik)
library(numDeriv)
library(lhs)

# Set the directory where save the results of the Gibbs sampler
# setwd("")

##############################################################################################
##                              LIKELIHOOD OF GENERALIZED LOGISTIC
# likelihood of generalized logistic regression (p especially)
# with x y pre-specified
##############################################################################################
loglikelihood_p_xy <- function(param, y, X){
  
  # last one in param is: p
  n.beta <- length(param) - 1
  p      <- param[n.beta+1]
  
  x.beta <- X%*%param[1:n.beta]
  x.beta <- ifelse(x.beta>10, 10, x.beta)
  prob.p <- pbeta(expit(x.beta), p, p)
  L      <- sum(dbinom(y,1,prob.p,log=TRUE))
  -L
}
expit<-function(x){exp(x)/(1+exp(x))}

##############################################################################################
##                          NUMERICAL DERIVATIVE OF LIKELIHOOD FOR GENERALIZED
# numerically evaluate the derivative of likelihood function of generalized logit regression
##############################################################################################
loglike_p_numDev <- function(param, y, X){
  GradVal <- matrix(NA, length(y), length(param))
  for(i in 1:length(y)){
  GradVal[i,] <- grad(loglikelihood_p_xy , param, y=y[i], X=X[i,])
  }
  GradVal
}


##############################################################################################
##                          GIBBS SAMPLING FOR GENERALIZED LOGIT REGRESSION
##############################################################################################
Bayes.logistic.el<-function(y, X, beta.mn=rep(0,ncol(X)), beta.sd=1, hyper.a, hyper.b, n.samples){
  # EL sampler code for the model:
  # y[i]        ~ Bern(p[i])
  # Generalized-Logit(p[i]) = X%*%beta
  # beta[j]     ~ N(prior.mn,prior.sd)
  # prior for p ~ Gamma(hyper.a, hyper.b)
  nc <- ncol(X)
  n  <- length(y)
  K  = n.samples; # number of posterior samples
  L  = n.samples; # number of candidates for one posterior sample
  # latin hypercube sampling
  randomCDF = randomLHS(L, nc+1)
  beta.prop = matrix(NA, L, nc)
  for(j in 1:nc){
    beta.prop[ , j] = qnorm(randomCDF[ , j], mean = beta.mn[j], sd = beta.sd)
  }
  p.prop = qgamma(randomCDF[ , nc+1], hyper.a, hyper.b)
  
  prob = numeric(L)
  for (i in 1:L){
    el = try(el.test(loglike_p_numDev(c(beta.prop[i, ], p.prop[i]), y, X), rep(0, nc+1)), silent = TRUE)
    if ('try-error' %in% class(el)) {
      prob[i] = 0
    }else{
      prob[i] = exp(-el$'-2LLR'/2)
    }
  }
  pos             = sample(L, K, replace = TRUE, prob)
  beta.pos.sample = beta.prop[pos, ]
  p.pos.sample    = p.prop[pos]
  
  # Return the posterior sample beta and p, estimated parameters from bootstrap samples
  list(beta = beta.pos.sample, p=p.pos.sample, prop = prob, beta.prop = beta.prop, p.prop = p.prop)
}


##############################################################################################
##                              LIKELIHOOD OF GENERALIZED LOGISTIC
# empirical likelihood sampler code for standard logistic regression
##############################################################################################
log_like_std <- function(beta, y, X){
  xbeta  <- X%*%beta
  xbeta  <- ifelse(xbeta>10,10,xbeta)
  like   <- sum(dbinom(y,1,expit(xbeta),log=TRUE))
  -like
}
expit<-function(x){exp(x)/(1+exp(x))}

##############################################################################################
##                          NUMERICAL DERIVATIVE OF LIKELIHOOD FOR GENERALIZED
# numerically evaluate the derivative of likelihood function of standard logit regression
##############################################################################################
loglike_std_numDev <- function(param, y, X){
  GradVal <- matrix(NA, length(y), length(param))
  for(i in 1:length(y)){
    GradVal[i,] <- grad(log_like_std , param, y=y[i], X=X[i,])
  }
  GradVal
}


##############################################################################################
##                          GIBBS SAMPLING FOR STANDARD LOGIT REGRESSION
##############################################################################################
Bayes.std.logistic.el<-function(y, X, beta.mn=rep(0,ncol(X)), beta.sd=1, n.samples){
  # EL sampler code for the model:
  # y[i]        ~ Bern(p[i])
  # Logit(p[i]) = X%*%beta
  # beta[j]     ~ N(prior.mn,prior.sd)
  nc <- ncol(X)
  n  <- length(y)
  K  = n.samples; # number of posterior samples
  L  = n.samples; # number of candidates for one posterior sample
  # latin hypercube sampling
  randomCDF = randomLHS(L, nc)
  beta.prop = matrix(NA, L, nc)
  for(j in 1:nc){
    beta.prop[ , j] = qnorm(randomCDF[ , j], mean = beta.mn[j], sd = beta.sd)
  }
  
  prob = numeric(L)
  for (i in 1:L){
    el = try(el.test(loglike_std_numDev(c(beta.prop[i, ]), y, X), rep(0, nc)), silent = TRUE)
    if ('try-error' %in% class(el)) {
      prob[i] = 0
    }else{
      prob[i] = exp(-el$'-2LLR'/2)
    }
  }
  pos             = sample(L, K, replace = TRUE, prob)
  beta.pos.sample = beta.prop[pos, ]
  
  # Return the posterior sample beta and p, estimated parameters from bootstrap samples
  list(beta = beta.pos.sample, prop = prob, beta.prop = beta.prop)
}


##############################################################################################
##                              MAIN PART TO BE RUN
# Generate the different simulation experiments
##############################################################################################
for(kk in 1:20){
  set.seed(kk)
  p         <- 0.1               # initial value of p
  true.beta <- c(1,-1,-3,1,3)    # initial value of beta
  dim.x     <- length(true.beta) # Number of coefficinets
  param_mle = rep(0, dim.x+1)
  
  n         <- 500
  x         <- cbind(1,matrix(rnorm(n*4), n, 4))
  true.p    <- exp(x%*%true.beta)/(1+exp(x%*%true.beta))
  y         <- rbinom(n, 1, true.p)
  xbetatrue <- x%*%true.beta
  x1        <- rbeta(n,p,p)    # random beta(p,p)
  z12       <- log(x1/(1-x1))+ xbetatrue
  y         <- ifelse(z12<0,0,1)  # our new variable if z<0 then y=0 else y=1
    
  # Generate simulation from Binomial distribution
  True.p2 <- pbeta(exp(x%*%true.beta)/(1+exp(x%*%true.beta)), p, p)
  y2      <- rbinom(n, 1, true.p2)
    
  # fit the model
  n.samples <- 20000     # number of posterior samples
  beta.mn   <- true.beta # prior mean for beta
  beta.sd   <- 1         # prior std deviation for beta
  hyper.a   <- 1         # prior hyperparameters for p
  hyper.b   <- 1         # prior hyperparameters for p
    
  # fit the data with generalized logit
  t0  = Sys.time()  
  fit <- Bayes.logistic.el(y = y2, X = x, beta.mn = true.beta, beta.sd = 1,
                              hyper.a, hyper.b, n.samples)
  runtime <- Sys.time() - t0
    
  # fit the data with std logistic
  n2.samples <- 10000 # number of posterior samples
  fit.std.el <- Bayes.std.logistic.el(y = y2, X = x, beta.mn = true.beta, beta.sd = 1, n2.samples)
}


