##simulation comparing weighted KM, two IPWs and two DRs for transporting
##under random censoring for weak and strong sampling 

rm(list=ls())
library(survival)
library(survey)

#Set directory to the "R Code" folder
# setwd(".../R Code") # change working directory to the current directory
#setwd("C:/Users/user/Dropbox/research/generability/R_code_generability_bai/submission")

N <- 1000000
nsim = 1000
cen.rate <- rep(0,nsim)
ipw.su.est.ATE <- matrix(0,ncol=3,nrow=nsim)
ipw.su.se.ATE <- matrix(0,ncol=3,nrow=nsim)
ipw2.su.est.ATE <- matrix(0,ncol=3,nrow=nsim)
ipw2.su.se.ATE <- matrix(0,ncol=3,nrow=nsim)

dr.cox.unst.est.ATE <- matrix(0,ncol=3,nrow=nsim)
dr.cox.unst.se.ATE <- matrix(0,ncol=3,nrow=nsim)
dr2.cox.unst.est.ATE <- matrix(0,ncol=3,nrow=nsim)
dr2.cox.unst.se.ATE <- matrix(0,ncol=3,nrow=nsim)

km.unst.est.ATE <- matrix(0,ncol=3,nrow=nsim)
km.unst.se.ATE <- matrix(0,ncol=3,nrow=nsim)

check.v <- matrix(0,ncol=3,nrow=nsim)
cover.su.IPCW <- matrix(0,ncol=3*2,nrow=nsim) 
cover.cox.DR <- matrix(0,ncol=3,nrow=nsim)
cover.cox.DR2 <- matrix(0,ncol=3,nrow=nsim)
cover.km <- matrix(0,ncol=3,nrow=nsim)

########################################################################
#sindex is the index(indices) for covariates used in the logistic model;        
#cindex is the index(indices) for covariates used in the Cox model for censoring time; 
########################################################################
#scenario=1,5,4,6 (or 2) is for IPW
#scenario=1,2,3,6 is for DR
#scenario=1(or 3, or 4), 2(or 5, or 6) is for WKM
scenario = 2
if(scenario == 1){ #ccc model   
  cindex = c(1,2,3); sindex = c(1,2,3); condtype = "cox";
}else if(scenario == 2){ #wwc model
  cindex = c(1); sindex = c(1); condtype = "cox";
}else if(scenario == 3){ #ccw model
  cindex = c(1,2,3); sindex = c(1,2,3); condtype = "lognorm";
}else if(scenario == 4){ #wcc model
  cindex = c(1); sindex = c(1,2,3); condtype = "cox";
}else if(scenario ==5){  #cwc model
  cindex = c(1,2,3); sindex = c(1); condtype = "cox";
}else{  #www model
  cindex = c(1); sindex = c(1); condtype = "lognorm";
}

# three target survival times
times <- c(0.1281846, 0.3312968, 0.7124819)
choice = 2					
if(choice==1){						
  # selection probability for trial						
  gampar <- c(-6.5,-0.2,-0.1,-0.5)  #weak sampling						
  true_val <- c(-0.07714378,-0.14102360,-0.17789499)  						
}else{						
  # selection probability for trial						
  gampar <- c(-6.7,-0.8,-0.4,-1)  #strong sampling						
  true_val <- c(-0.07711887,-0.14097153,-0.17781734) #true value						
}

###main simulation
for(v in 1:nsim){
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
  #cpop <- -log(runif(N))/(lambda_cpop)
  #random censoring case
  cpop = runif(N, min = 0, max=2.1)
  # Sampling mechanism
  samp <-function(par,x){
    exp(x %*% par)/(1+exp(x %*% par))
  }
  p_selS <- samp(gampar,cbind(1,xpop))
  S <- rbinom(N,1,p_selS)
  # population data
  pop_all <-cbind(tpop_trt,tpop_crl,cpop,xpop)
  # covariate for rct
  rctsamp<- pop_all[which(S==1),]
  n <- sum(S==1)
  # treatment variable 
  z <- rbinom(n,1,0.5) 
  # observed time and censoring status
  timerct <- z*rctsamp[,1]+(1-z)*rctsamp[,2]
  crct <- rctsamp[,3]
  obs <- pmin(timerct,crct)
  delta<- as.numeric(timerct <= crct)
  cen.rate[v] <- mean(delta==0) # censoring rate
  xx.rct <- rctsamp[,4:6]
  # Generate survey sample
  M <-sum(S==0) # M = N-n
  xpopS0 <- pop_all[which(S==0),4:6] # N-n population
  rownames(xpopS0) <- which(S==0)
  # logistic model for D
  etapar <- c(-7,0.3,0.4,0.2)
  p_selD <- samp(etapar,cbind(1,xpopS0))
  D <- rbinom(M,1,p_selD)
  # data from survey
  svx <- xpopS0[which(D==1),]
  # data in the rct
  datrct <- data.frame(obs=obs,delta=delta,z=z, x1rct = rctsamp[,4],x2rct = rctsamp[,5],x3rct = rctsamp[,6])
  
  # estimation of sampling score
  resp <- c(S[S==1],S[S==0][D==1])
  xall <- rbind(rctsamp[,4:6],svx)
  sweight <- c(rep(1,n),1/(p_selD[which(D==1),]))
  dat.cb<- data.frame(resp=resp,x1cb = xall[,1],x2cb = xall[,2],x3cb = xall[,3],sweight=sweight)
  id <- 1:(dim(dat.cb)[1])
  m <- dim(svx)[1]
  #participation model
  fit.score <- glm(resp ~ xall[,sindex], weights = sweight,data=dat.cb, family = binomial)
  pred.score.tmp <- predict(fit.score,type="response")
  
  # p(X) for RCT sample
  K <- length(times)
  svw <- sum(sweight[(n+1):(n+m)])
  pred.score <- pred.score.tmp/(1-pred.score.tmp)
  w.h <- pred.score[1:n]
  print(c(n,m))
  
  #weighted KM
  ###########################################
  km.ATE <- function(obs.rct, delta.rct, z.rct, weights, n, taut){
    dat_rct = data.frame(id=1:n, time=obs.rct, delta=delta.rct, z=z.rct, weights=weights)
    drct <- svydesign(id=~1, weights=~weights, data=dat_rct)
    km_res <- svykm(Surv(time,delta)~z, se=TRUE, design=drct)  
    km_res1 = km_res$`1`  #for treatment
    km_res0 = km_res$`0`  #for control
    s1_hat = km_res1$surv
    obs1 = km_res1$time
    varlog_s1 = km_res1$varlog #var(log(s))
    var_s1_hat = (exp(log(s1_hat)))^2*varlog_s1  #delta method for variance(s)
    s0_hat = km_res0$surv
    obs0 = km_res0$time
    varlog_s0 = km_res0$varlog 
    var_s0_hat = (exp(log(s0_hat)))^2*varlog_s0 
    # return results
    u1 = obs1[order(obs1)]
    u0 = obs0[order(obs0)]
    est <- se <- NULL
    for(tt in taut){
      id1 = sum(u1 <= tt)
      id0 = sum(u0 <= tt)
      est <- c(est, s1_hat[id1]-s0_hat[id0])
      se <- c(se, sqrt(varlog_s1[id1]+var_s0_hat[id0]))
    }
    return(list(est=est, se=se))
  }
  
  weights = 1/w.h #sampling score weights for S=1
  km.est.ATE <- km.ATE(obs, delta, z, weights, n, times)
  km.unst.est.ATE[v,] <- km.est.ATE[[1]]
  km.unst.se.ATE[v,] <- km.est.ATE[[2]]
  # Coverage
  cover.km[v,1:3] <- as.numeric(((km.est.ATE[[1]]-qnorm(0.975)*km.est.ATE[[2]]) <= true_val) & 
                                  (true_val <= (km.est.ATE[[1]] + qnorm(0.975)*km.est.ATE[[2]])))
  
  
  # IPCW
  ##################################################################
  # Estimate K.prob with Cox model and otice that this model is treatment-specific,
  # and hence compatible with the Bai approach
  ##################################################################
  Gfunc_cox = function(obs.rct, delta.rct, z.rct, xx.rct, cindex, taut){
    n <- length(obs.rct)
    # ordering the data according to u increasing:
    u.order <- order(obs.rct)
    u <- obs.rct[u.order]
    delta <- delta.rct[u.order]
    z <- z.rct[u.order]
    
    xc <- xx.rct[,cindex]  #covariates for censoring
    lcindex = length(cindex)
    if(lcindex==1) xc <- xc[u.order] else xc <- xc[u.order,]
    #########################################
    ### estimate the proportional hazard  model parameters 
    ### for both failure time and cesoring time 
    #########################################
    uz <- u[z!=0]   # using data only when z=1
    deltaz <- delta[z!=0]
    if(lcindex==1) { # fit cox ph model for censoring
      bc_est <- coxph(Surv(uz, 1-deltaz) ~ xc[z!=0])$coef
    }  else {
      bc_est <- coxph(Surv(uz, 1-deltaz) ~ xc[z!=0,])$coef
    }
    
    #######################################
    ### computing propotional hazards   ###
    ### model for both failure time     ###
    ### and censoring, denoted by       ###
    ### dls(j) and dlc(j) for time u(j) ###
    #######################################
    if(lcindex==1){
      s2 <- cumsum(z[n:1]*exp(as.matrix(xc[n:1])%*%bc_est))
    } else{
      s2 <- cumsum(z[n:1]*exp(as.matrix(xc[n:1,])%*%bc_est))
    }
    dlc <- z*(1-delta)/s2[n:1]
    dlc[(1:length(dlc))[is.nan(dlc)]] <- 0
    
    ###########################################
    ### for patient i with covariates x[i], ###
    ### compute estimated survival time and ###
    ### censoring at time point u[j]        ###
    ###########################################
    # sss(i, j) -> H(u_j, X_i); ssc -> K(u_j, X_i)
    ssc <- matrix(0, n, n)
    # each column is one time point, row an individual
    ssc[, 1] <- dlc[1]*exp(as.matrix(xc)%*%bc_est)
    otc <- outer(c(exp(as.matrix(xc)%*%bc_est)), dlc, "*")   
    otc[, 1] <- rep(0, n)
    otc_cum <- t(apply(otc, 1, cumsum))
    ssc_1 <- outer(ssc[, 1], rep(1, n))
    ssc <- ssc_1 + otc_cum
    ssc <- exp(-ssc)
    #ssc -> K(u_j, X_i)
    d.ssc <- diag(ssc)
    K.res <- rep(NA,n)
    #we need to return to original observed time 
    K.res[u.order] <- d.ssc #K(U|X)
    return(K.res)
  }
  
  ####IPW estimation
  ipw.ATE <- function(obs,delta,z,xx,cindex,taut){
    est.ipw.eff1u <- est.ipw.eff2u <- rep(NA, K)
    se.ipw.eff1u <- se.ipw.eff2u <- rep(NA, K)
    
    K.prob1u <- Gfunc_cox(obs,delta,z,xx,cindex,taut) #K(U|X), vector, length is n
    K.prob0u <- Gfunc_cox(obs,delta,1-z,xx,cindex,taut) #K(U|X), vector, length is n
    
    w.h <- pred.score[1:n]
    ######################################################
    # stabilizing term (should subset S=1)
    st_z1 <- sum(z/w.h)
    st_z0 <- sum((1-z)/w.h)
    
    # point estimation
    for(j in 1:K){
      IUt <- delta * as.numeric(obs >= taut[j]) #have delta
      ift1u <- c((IUt*z)/(K.prob1u*w.h),rep(0,m))
      ift0u <- c((IUt*(1-z))/(K.prob0u*w.h),rep(0,m))
      mu1u <- 1/svw*sum(2*ift1u)
      mu0u <- 1/svw*sum(2*ift0u)
      est.ipw.eff1u[j] <- mu1u-mu0u
      ift1s1u = 2*ift1u-c(rep(0,n),sweight[(n+1):(n+m)])*mu1u
      ift0s1u = 2*ift0u-c(rep(0,n),sweight[(n+1):(n+m)])*mu0u
      se.ipw.eff1u[j] <- sqrt(sum((ift1s1u-ift0s1u-est.ipw.eff1u[j])^2)/svw^2)
      
      mu1.su <- sum(ift1u)/st_z1
      mu0.su <- sum(ift0u)/st_z0
      ift1.su <- c((z*(IUt/K.prob1u-mu1.su)/w.h), rep(0,m))
      ift0.su <- c(((1-z)*(IUt/K.prob0u-mu0.su)/w.h), rep(0,m))
      est.ipw.eff2u[j] <- mu1.su-mu0.su
      se.ipw.eff2u[j] <- sqrt(sum(ift1.su^2)/(st_z1^2)+sum(ift0.su^2)/(st_z0^2))
    }
    ipw.res = list(est.ipw.eff1u, est.ipw.eff2u, se.ipw.eff1u, se.ipw.eff2u)
    return(ipw.res)	
  }
  
  ####estimation results
  ipw.est.ATE <- ipw.ATE(obs,delta,z,xx.rct,cindex,times)
  # IPCW
  ipw.su.est.ATE[v,] <- ipw.est.ATE[[1]] #using K(u|X)
  ipw2.su.est.ATE[v,] <- ipw.est.ATE[[2]] 
  
  ipw.su.se.ATE[v,] <- ipw.est.ATE[[3]] 
  ipw2.su.se.ATE[v,] <- ipw.est.ATE[[4]] 
  
  # Coverage
  cover.su.IPCW[v,1:3] <- as.numeric((ipw.est.ATE[[1]] - qnorm(0.975)*ipw.est.ATE[[3]] <= true_val) & (true_val <= ipw.est.ATE[[1]] + qnorm(0.975)*ipw.est.ATE[[3]]))
  cover.su.IPCW[v,4:6] <- as.numeric((ipw.est.ATE[[2]] - qnorm(0.975)*ipw.est.ATE[[4]] <= true_val) & (true_val <= ipw.est.ATE[[2]] + qnorm(0.975)*ipw.est.ATE[[4]]))
  
  #DR estimator
  xx.svs <- svx
  sweight <- as.numeric(sweight)
  source("sim_code_taste_final.R")
  dr.ATE <- function(obs.rct, delta.rct, z.rct, xx.rct, xx.svs, taut, cindex, condtype){
    # estimate fot z=1:
    r_z1 <- survres_z(obs.rct, delta.rct, z.rct, xx.rct, xx.svs, taut, cindex, condtype)
    S_hat_z1 <- r_z1$surv   # s is used to estimate the S_hat (survival function estimator)
    sif_z1 <- r_z1$s_if   # ssq is used to estimate the v_hat (variance of S_hat)
    
    # estimate fot z=0:
    r_z0 <- survres_z(obs.rct, delta.rct, 1-z.rct, xx.rct, xx.svs, taut, cindex, condtype)
    S_hat_z0 <- r_z0$surv   # s is used to estimate the S_hat (survival function estimator)
    sif_z0 <- r_z0$s_if
    
    # return results
    u.order <- order(obs.rct)
    u <- obs.rct[u.order]
    est <- se <- NULL
    for(tt in taut){
      id = sum(u<=tt)
      sdiff = S_hat_z1[id]-S_hat_z0[id]
      est <- c(est, sdiff)
      se <- c(se, sqrt(sum((sif_z1[,id] - sif_z0[,id]- sdiff)^2)/svw^2)) 
    }
    return(list(est=est, se=se))
  }
  dr.est.ATE.cox <- dr.ATE(obs,delta,z,xx.rct,xx.svs,times,cindex,condtype)
  
  # DR1
  dr.cox.unst.est.ATE[v,] <- dr.est.ATE.cox[[1]]
  dr.cox.unst.se.ATE[v,] <- dr.est.ATE.cox[[2]]
  # Coverage
  cover.cox.DR[v,1:3] <- as.numeric(((dr.est.ATE.cox[[1]]-qnorm(0.975)*dr.est.ATE.cox[[2]]) <= true_val) & 
                                      (true_val <= (dr.est.ATE.cox[[1]] + qnorm(0.975)*dr.est.ATE.cox[[2]])))
  
  dr2.ATE <- function(obs.rct, delta.rct, z.rct, xx.rct, xx.svs, taut, cindex, condtype){
    # estimate fot z=1:
    r_z1 <- survres_z2(obs.rct, delta.rct, z.rct, xx.rct, xx.svs, taut, cindex, condtype)
    S_hat_z1 <- r_z1$surv   # s is used to estimate the S_hat (survival function estimator)
    sif_z1 <- r_z1$s_if   # ssq is used to estimate the v_hat (variance of S_hat)
    
    # estimate fot z=0:
    r_z0 <- survres_z2(obs.rct, delta.rct, 1-z.rct, xx.rct, xx.svs, taut, cindex, condtype)
    S_hat_z0 <- r_z0$surv   # s is used to estimate the S_hat (survival function estimator)
    sif_z0 <- r_z0$s_if
    
    # return results
    u.order <- order(obs.rct)
    u <- obs.rct[u.order]
    est <- se <- NULL
    for(tt in taut){
      id = sum(u<=tt)
      est <- c(est, S_hat_z1[id]-S_hat_z0[id])
      se <- c(se, sqrt(sum((sif_z1[,id])^2) + sum((sif_z0[,id])^2)))
    }
    return(list(est=est, se=se))
  }
  
  dr2.est.ATE.cox <- dr2.ATE(obs,delta,z,xx.rct,xx.svs,times,cindex, condtype)
  
  # DR2
  dr2.cox.unst.est.ATE[v,] <- dr2.est.ATE.cox[[1]]
  dr2.cox.unst.se.ATE[v,] <- dr2.est.ATE.cox[[2]]
  # Coverage
  cover.cox.DR2[v,1:3] <- as.numeric(((dr2.est.ATE.cox[[1]]-qnorm(0.975)*dr2.est.ATE.cox[[2]]) <= true_val) & 
                                       (true_val <= (dr2.est.ATE.cox[[1]] + qnorm(0.975)*dr2.est.ATE.cox[[2]])))
  
  print(v)
}
#### summarize
colMeans(ipw.su.est.ATE) - true_val
colMeans(ipw2.su.est.ATE) - true_val
colMeans(ipw.su.se.ATE) 
colMeans(ipw2.su.se.ATE) 
apply(ipw.su.est.ATE, 2, sd)
apply(ipw2.su.est.ATE, 2, sd)
colMeans(cover.su.IPCW)

#dr
colMeans(dr.cox.unst.est.ATE) - true_val
colMeans(dr2.cox.unst.est.ATE) - true_val
colMeans(dr.cox.unst.se.ATE) 
colMeans(dr2.cox.unst.se.ATE) 

apply(dr.cox.unst.est.ATE,2,sd)
apply(dr2.cox.unst.est.ATE,2,sd)

colMeans(cover.cox.DR)
colMeans(cover.cox.DR2)

##km
colMeans(km.unst.est.ATE) - true_val
colMeans(km.unst.se.ATE)
apply(km.unst.est.ATE,2,sd)
colMeans(cover.km)


bias.ipcw1u <- colMeans(ipw.su.est.ATE) - true_val
bias.ipcw2u <- colMeans(ipw2.su.est.ATE) - true_val

se.ipcw1u <- colMeans(ipw.su.se.ATE) 
se.ipcw2u <- colMeans(ipw2.su.se.ATE) 

emp.ipcw1u <- apply(ipw.su.est.ATE, 2, sd)
emp.ipcw2u <- apply(ipw2.su.est.ATE, 2, sd)

cover.ipcw1u <- colMeans(cover.su.IPCW)[1:3]
cover.ipcw2u <- colMeans(cover.su.IPCW)[4:6]

#DR
bias.cox.dr1 <- colMeans(dr.cox.unst.est.ATE) - true_val
bias.cox.dr2 <- colMeans(dr2.cox.unst.est.ATE) - true_val
emp.cox.dr1 <- apply(dr.cox.unst.est.ATE,2,sd)
emp.cox.dr2 <- apply(dr2.cox.unst.est.ATE,2,sd)
se.cox.dr1 <- colMeans(dr.cox.unst.se.ATE)
se.cox.dr2 <- colMeans(dr2.cox.unst.se.ATE)
cover.cox.dr1 <- colMeans(cover.cox.DR)
cover.cox.dr2 <- colMeans(cover.cox.DR2)

#weighted KM
bias.wkm = colMeans(km.unst.est.ATE) - true_val
emp.wkm = apply(km.unst.est.ATE,2,sd)
se.wkm = colMeans(km.unst.se.ATE)
cover.wkm = colMeans(cover.km)


###summary results
ipcw1u = c(bias.ipcw1u[1],emp.ipcw1u[1],se.ipcw1u[1],cover.ipcw1u[1],bias.ipcw1u[1],emp.ipcw1u[2],
           se.ipcw1u[2],cover.ipcw1u[2],bias.ipcw1u[3],emp.ipcw1u[3],se.ipcw1u[3],cover.ipcw1u[3])
ipcw2u = c(bias.ipcw2u[1],emp.ipcw2u[1],se.ipcw2u[1],cover.ipcw2u[1],bias.ipcw2u[1],emp.ipcw2u[2],
           se.ipcw2u[2],cover.ipcw2u[2],bias.ipcw2u[3],emp.ipcw2u[3],se.ipcw2u[3],cover.ipcw2u[3])


dr1 = c(bias.cox.dr1[1],emp.cox.dr1[1],se.cox.dr1[1],cover.cox.dr1[1],bias.cox.dr1[2],emp.cox.dr1[2],
        se.cox.dr1[2],cover.cox.dr1[2],bias.cox.dr1[3],emp.cox.dr1[3],se.cox.dr1[3],cover.cox.dr1[3])

dr2 = c(bias.cox.dr2[1],emp.cox.dr2[1],se.cox.dr2[1],cover.cox.dr2[1],bias.cox.dr2[2],emp.cox.dr2[2],
        se.cox.dr2[2],cover.cox.dr2[2],bias.cox.dr2[3],emp.cox.dr2[3],se.cox.dr2[3],cover.cox.dr2[3])

wkm = c(bias.wkm[1],emp.wkm[1],se.wkm[1],cover.wkm[1],bias.wkm[2],emp.wkm[2],
        se.wkm[2],cover.wkm[2],bias.wkm[3],emp.wkm[3],se.wkm[3],cover.wkm[3])

resu_final = round(rbind(ipcw1u,ipcw2u,dr1,dr2,wkm),4)
colnames(resu_final) = c("bias1","em.se1","se1","cr1","bias2","em.se2","se2","cr2",
                         "bias3","em.se3","se3","cr3")

resu_final
mean(cen.rate)
write.csv(resu_final,file="resu1.csv")

#library(xtable)
#xtable(resu_final)
