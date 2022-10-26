##generate true values for delta(t)=S_1(t|S=0)-S_0(t|S=0) in simulations

rm(list=ls())
library(survival)

N <- 1000000
times <- c(0.1281846, 0.3312968, 0.7124819)
nsim = 1000
true_eff = matrix(0,nsim,ncol=3)
scenario = 1

for(v in 1:nsim){
  # Generate covariates
  set.seed(1000000+v)
  # Generate covariates
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  x3 <- rbinom(N,1,0.5)
  xpop <- cbind(x1,x2,x3)
  # generate failure time for treatment and control
  upop_trt <- runif(N)
  upop_crl <- runif(N)
  lambda_trt <- exp(xpop %*% c(0.2,-0.6,0.6))
  lambda_crl <- exp(xpop %*% c(0.2,0.6,-0.6))
  lambda_cpop <- exp(xpop %*% c(0.1,0.4,-0.6))
  tpop_trt <- -log(upop_trt)/(lambda_trt)
  tpop_crl <- -log(upop_crl)/(lambda_crl)	
  if(scenario==1){ #weak sampling and covariate dependent censoring
    cpop <- -log(runif(N))/(lambda_cpop)
    gampar <- c(-6.5,-0.2,-0.1,-0.5) 
  }else if(scenario==2){ #strong sampling and covariate dependent censoring
    cpop <- -log(runif(N))/(lambda_cpop)
    gampar <- c(-6.7,-0.8,-0.4,-1)
  }else if(scenario==3){
    cpop = runif(N,min = 0,max=2) #weak sampling and fully independent censoring
    gampar <- c(-6.5,-0.2,-0.1,-0.5) 
  }else{
    cpop = runif(N,min = 0,max=2) #strong sampling and fully independent censoring
    gampar <- c(-6.7,-0.8,-0.4,-1)
  }
  # Sampling mechanism
  samp <-function(par,x){
    exp(x %*% par)/(1+exp(x %*% par))
  }
  # selection probability for trial
  p_selS <- samp(gampar,cbind(1,xpop))
  S <- rbinom(N,1,p_selS)
  tpop_trt_S0 = tpop_trt[S==0]
  tpop_crl_S0 = tpop_crl[S==0]
  for(i in 1:3){
    true_eff[v,i] = mean(c(tpop_trt_S0>=times[i])-c(tpop_crl_S0>=times[i]))
  }
  print(v)
}

#when scenario = 1 or 3
#-0.07714378,-0.14102360,-0.17789499
#when scenario = 2 or 4
#-0.07711887,-0.14097153,-0.17781734