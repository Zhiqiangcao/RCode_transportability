survres_z <- function(obs.rct, delta.rct, z.rct, xx.rct, xx.svs, taut, cindex, condtype){
  n <- length(obs)
  m <- nrow(xx.svs)
  
  # ordering the data according to u increasing:
  u.order <- order(obs.rct)
  u <- obs.rct[u.order]
  x <- xx.rct[u.order,]
  delta <- delta.rct[u.order]
  z <- z.rct[u.order]
  
  xc <- xx.rct[,cindex]  #covariates for censoring
  lcindex = length(cindex)
  if(lcindex==1) xc <- xc[u.order] else xc <- xc[u.order,] #may be only one covariate
  ##################################
  ### use logistic model to fit  ###
  ### the propensity probability ###
  ##################################
  #pred.score <- pred.score.tmp/(1-pred.score.tmp)
  pp_est <- as.numeric(pred.score[1:n][u.order])  #sampling score
  if(condtype == "lognorm"){
    #########################################
    ### estimate the proportional hazard  ###
    ### model parameters for both failure ###
    ### time and cesoring time            ###
    #########################################
    uz <- u[z!=0]   # using data only when z=1
    deltaz <- delta[z!=0]
    b_est_all <- survreg(Surv(uz, deltaz) ~ x[z!=0,], dist = "lognormal")
    b_est <- b_est_all$coef
    b_scale <- b_est_all$scale
    if(lcindex==1) { # fit cox ph model for censoring
      bc_est <- coxph(Surv(uz, 1-deltaz) ~ xc[z!=0])$coef
    }  else {
      bc_est <- coxph(Surv(uz, 1-deltaz) ~ xc[z!=0,])$coef
    }
  
    #######################################
    ### computing parametric AFT model  ###
    ### for failure time          ###
    ### and propotional hazards for censoring, denoted by       ###
    ### dls(j) and dlc(j) for time u(j) ###
    #######################################
    # s1 <- cumsum(z[n:1]*exp(as.matrix(x[n:1,])%*%b_est))
    # dls <- z*delta/s1[n:1]
    # dls[(1:length(dls))[is.nan(dls)]] <- 0
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
    sss <- matrix(0, n, n)
    res <- matrix(0, n, n)
    sss_svs <- matrix(0, m, n)
    res_svs <- matrix(0, m, n)
    ssc <- matrix(0, n, n)
    # each column is one time point, row an individual
    for(j in 1:n){
      res[,j] <- (log(u[j]) - as.matrix(cbind(1,x)) %*%b_est)/b_scale
      sss[,j] <- 1-pnorm(res[,j])		
    }
    for(j in 1:n){
      res_svs[,j] <- (log(u[j]) - as.matrix(cbind(1,xx.svs)) %*%b_est)/b_scale
      sss_svs[,j] <- 1-pnorm(res_svs[,j])		
    }
    ssc[, 1] <- dlc[1]*exp(as.matrix(xc)%*%bc_est)
    otc <- outer(c(exp(as.matrix(xc)%*%bc_est)), dlc, "*")   
    otc[, 1] <- rep(0, n)
    otc_cum <- t(apply(otc, 1, cumsum))
    ssc_1 <- outer(ssc[, 1], rep(1, n))
    ssc <- ssc_1 + otc_cum
    ssc <- exp(-ssc)
  }
  if(condtype == "cox"){
    #########################################
    ### estimate the proportional hazard  ###
    ### model parameters for both failure ###
    ### time and cesoring time            ###
    #########################################
    uz <- u[z!=0]   # using data only when z=1
    deltaz <- delta[z!=0]
    b_est <- coxph(Surv(uz, deltaz) ~ x[z!=0,])$coef   # fit cox ph model for surv time
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
    s1 <- cumsum(z[n:1]*exp(as.matrix(x[n:1,])%*%b_est))
    dls <- z*delta/s1[n:1]
    dls[(1:length(dls))[is.nan(dls)]] <- 0
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
    sss <- matrix(0, n, n)
    sss_svs <- matrix(0, m, n)
    ssc <- matrix(0, n, n)
    # each column is one time point, row an individual
    sss[, 1] <- dls[1]*exp(as.matrix(x)%*%b_est)
    sss_svs[, 1] <- dls[1]*exp(as.matrix(xx.svs)%*%b_est)
    ssc[, 1] <- dlc[1]*exp(as.matrix(xc)%*%bc_est)
    
    ots <- outer(c(exp(as.matrix(x)%*%b_est)), dls, "*")   
    ots[, 1] <- rep(0, n)
    ots_cum <- t(apply(ots, 1, cumsum))
    sss_1 <- outer(sss[, 1], rep(1, n))
    sss <- sss_1 + ots_cum
    
    otsvs <- outer(c(exp(as.matrix(xx.svs)%*%b_est)), dls, "*")   
    otsvs[, 1] <- rep(0, m)
    otsvs_cum <- t(apply(otsvs, 1, cumsum))
    sss_1_svs <- outer(sss_svs[, 1], rep(1, n))
    sss_svs <- sss_1_svs + otsvs_cum
    
    otc <- outer(c(exp(as.matrix(xc)%*%bc_est)), dlc, "*")   
    otc[, 1] <- rep(0, n)
    otc_cum <- t(apply(otc, 1, cumsum))
    ssc_1 <- outer(ssc[, 1], rep(1, n))
    ssc <- ssc_1 + otc_cum
    sss <- exp(-sss)
    sss_svs <- exp(-sss_svs)
    ssc <- exp(-ssc)
  }
  ################################################
  ### compute the martingale integral part,    ###
  ### denote the augmentation term for patient ###
  ### i at time u[j] by q[i, j]                ###
  ################################################
  q <- matrix(0,n,n)
  ld <- lower.tri(matrix(1, n, n), diag=TRUE)   # lower triangular matrix
  ud <- matrix(1, n, n)-ld   # upper triangular matrix
  otr <- outer(c(exp(as.matrix(xc)%*%bc_est)), dlc)/(sss*ssc)  
  q <- t(apply(otr, 1, cumsum)) # \int{lambda_c(r,X)/[K_c(r,X)*H(r,X)]}
  q <- q*ld+diag(q)*ud  # \int{Y(r)*lambda_c(r,X)/[K_c(r,X)*H(r,X)]}
  
  ###############################
  ### estimate survival dist  ###
  ### at each time point u[j] ###
  ###############################
  surv <- rep(0, n)
  s_if <- matrix(0, n+m, n)
  
  temp1 = z*(1-delta)*sss/(diag(sss)*diag(ssc)*c(0.5*pp_est))   # Z*H/pi*\int{dN_c/[K*H]}
  ou <- -z*sss/c(0.5*pp_est) - z*sss*q/c(0.5*pp_est) + temp1
  ol <- z/(c(0.5*pp_est)*ssc) - z*sss/c(0.5*pp_est) - z*sss*q/c(0.5*pp_est)

  temp <- ou*ud+ol*ld+diag(diag(temp1)) 
  for(j in 1:n){
    tempj = temp[,j]
    tempj[is.na(tempj)] = 0
    temp[,j] = tempj
  }
  
  temp_svs <- sweight[-(1:n)]*sss_svs
  temp_all <- rbind(temp,temp_svs)
  surv <- apply(temp_all, 2, sum)/svw #S\hat S_a^{DR1}(t) for t=t_1,t_2,...,t_n
  
  s_if_rct = temp   #n*n matrix
  s_if_svy = sweight[-(1:n)]*(sss_svs-matrix(rep(surv,m),byrow=T,m,n)) #m*n matrix
  s_if = rbind(s_if_rct,s_if_svy)
  
  list <- list(surv=surv, s_if=s_if)
  return(list)
}

survres_z2 <- function(obs.rct, delta.rct, z.rct, xx.rct, xx.svs, taut, cindex, condtype){
  n <- length(obs)
  m <- nrow(xx.svs)
  
  # ordering the data according to u increasing:
  u.order <- order(obs.rct)
  u <- obs.rct[u.order]
  x <- xx.rct[u.order,]
  delta <- delta.rct[u.order]
  z <- z.rct[u.order]
  
  xc <- xx.rct[,cindex]  #covariates for censoring
  lcindex = length(cindex)
  if(lcindex==1) xc <- xc[u.order] else xc <- xc[u.order,]
  ##################################
  ### use logistic model to fit  ###
  ### the propensity probability ###
  ##################################
  pp_est <- as.numeric(pred.score[1:n][u.order])
  if(condtype == "lognorm"){
    #########################################
    ### estimate the proportional hazard  ###
    ### model parameters for both failure ###
    ### time and cesoring time            ###
    #########################################
    uz <- u[z!=0]   # using data only when z=1
    deltaz <- delta[z!=0]
    b_est_all <- survreg(Surv(uz, deltaz) ~ x[z!=0,], dist = "lognormal")
    b_est <- b_est_all$coef
    b_scale <- b_est_all$scale
    if(lcindex==1) { # fit cox ph model for censoring
      bc_est <- coxph(Surv(uz, 1-deltaz) ~ xc[z!=0])$coef
    }  else {
      bc_est <- coxph(Surv(uz, 1-deltaz) ~ xc[z!=0,])$coef
    }
    
    #######################################
    ### computing parametric AFT model  ###
    ### for failure time          ###
    ### and propotional hazards for censoring, denoted by       ###
    ### dls(j) and dlc(j) for time u(j) ###
    #######################################
    # s1 <- cumsum(z[n:1]*exp(as.matrix(x[n:1,])%*%b_est))
    # dls <- z*delta/s1[n:1]
    # dls[(1:length(dls))[is.nan(dls)]] <- 0
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
    sss <- matrix(0, n, n)
    res <- matrix(0, n, n)
    sss_svs <- matrix(0, m, n)
    res_svs <- matrix(0, m, n)
    ssc <- matrix(0, n, n)
    # each column is one time point, row an individual
    for(j in 1:n){
      res[,j] <- (log(u[j]) - as.matrix(cbind(1,x)) %*%b_est)/b_scale
      sss[,j] <- 1-pnorm(res[,j])		
    }
    for(j in 1:n){
      res_svs[,j] <- (log(u[j]) - as.matrix(cbind(1,xx.svs)) %*%b_est)/b_scale
      sss_svs[,j] <- 1-pnorm(res_svs[,j])		
    }
    ssc[, 1] <- dlc[1]*exp(as.matrix(xc)%*%bc_est)
    
    otc <- outer(c(exp(as.matrix(xc)%*%bc_est)), dlc, "*")   
    otc[, 1] <- rep(0, n)
    otc_cum <- t(apply(otc, 1, cumsum))
    ssc_1 <- outer(ssc[, 1], rep(1, n))
    ssc <- ssc_1 + otc_cum
    ssc <- exp(-ssc)
  }
  if(condtype == "cox"){
    #########################################
    ### estimate the proportional hazard  ###
    ### model parameters for both failure ###
    ### time and cesoring time            ###
    #########################################
    uz <- u[z!=0]   # using data only when z=1
    deltaz <- delta[z!=0]
    b_est <- coxph(Surv(uz, deltaz) ~ x[z!=0,])$coef  
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
    s1 <- cumsum(z[n:1]*exp(as.matrix(x[n:1,])%*%b_est))
    dls <- z*delta/s1[n:1]
    dls[(1:length(dls))[is.nan(dls)]] <- 0
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
    sss <- matrix(0, n, n)
    sss_svs <- matrix(0, m, n)
    ssc <- matrix(0, n, n)
    # each column is one time point, row an individual
    sss[, 1] <- dls[1]*exp(as.matrix(x)%*%b_est)
    sss_svs[, 1] <- dls[1]*exp(as.matrix(xx.svs)%*%b_est)
    ssc[, 1] <- dlc[1]*exp(as.matrix(xc)%*%bc_est)
    
    ots <- outer(c(exp(as.matrix(x)%*%b_est)), dls, "*")   
    ots[, 1] <- rep(0, n)
    ots_cum <- t(apply(ots, 1, cumsum))
    sss_1 <- outer(sss[, 1], rep(1, n))
    sss <- sss_1 + ots_cum
    
    otsvs <- outer(c(exp(as.matrix(xx.svs)%*%b_est)), dls, "*")   
    otsvs[, 1] <- rep(0, m)
    otsvs_cum <- t(apply(otsvs, 1, cumsum))
    sss_1_svs <- outer(sss_svs[, 1], rep(1, n))
    sss_svs <- sss_1_svs + otsvs_cum
    
    otc <- outer(c(exp(as.matrix(xc)%*%bc_est)), dlc, "*")   
    otc[, 1] <- rep(0, n)
    otc_cum <- t(apply(otc, 1, cumsum))
    ssc_1 <- outer(ssc[, 1], rep(1, n))
    ssc <- ssc_1 + otc_cum
    sss <- exp(-sss)
    sss_svs <- exp(-sss_svs)
    ssc <- exp(-ssc)
  }
  ################################################
  ### compute the martingale integral part,    ###
  ### denote the augmentation term for patient ###
  ### i at time u[j] by q[i, j]                ###
  ################################################
  q <- matrix(0,n,n)
  ld <- lower.tri(matrix(1, n, n), diag=TRUE)   # lower triangular matrix
  ud <- matrix(1, n, n)-ld   # upper triangular matrix
  otr <- outer(c(exp(as.matrix(xc)%*%bc_est)), dlc)/(sss*ssc)  
  # lambda_c(r,X)/[K_c(r,X)*H(r,X)]
  q <- t(apply(otr, 1, cumsum)) # \int{lambda_c(r,X)/[K_c(r,X)*H(r,X)]}
  q <- q*ld+diag(q)*ud  # \int{Y(r)*lambda_c(r,X)/[K_c(r,X)*H(r,X)]}
  
  ###############################
  ### estimate survival dist  ###
  ### at each time point u[j] ###
  ###############################
  
  stfactor <- sum(z/pp_est)
  surv <- rep(0, n)
  s_if <- matrix(0, n+m, n)
  
  temp4 = (1/stfactor)*z*(1-delta)*sss/(diag(sss)*diag(ssc)*c(pp_est))   # Z*H/pi*\int{dN_c/[K*H]}     
  ou <- -(1/stfactor)*z*sss/pp_est - (1/stfactor)*z*sss*q/pp_est + temp4
  ol <- (1/stfactor)*z/(pp_est*ssc) - (1/stfactor)*z*sss/pp_est - (1/stfactor)*z*sss*q/pp_est
  
  temp <- ou*ud+ol*ld+diag(diag(temp4))
  for(j in 1:n){
    tempj = temp[,j]
    tempj[is.na(tempj)] = 0
    temp[,j] = tempj
  }
  temp_svs <- sweight[-(1:n)]*sss_svs/svw
  temp_all <- rbind(temp,temp_svs)
  surv <- apply(temp_all, 2, sum)
  
  v1 = apply(temp,2,sum)
  v2 = apply(temp_svs,2,sum)
  
  mv1 = matrix(rep(v1,n),byrow=T,n,n) #n*n matrix
  s_if_rct = temp - (1/stfactor)*(z/pp_est)*mv1
  mv2 = matrix(rep(v2,m),byrow=T,m,n)  #m*n matrix
  s_if_svy = sweight[-(1:n)]*(sss_svs - mv2)/svw
  s_if = rbind(s_if_rct,s_if_svy)
  
  list <- list(surv=surv, s_if=s_if)
  return(list)
}


survres_z2app <- function(obs.rct, delta.rct, z.rct, xx.rct.censor, xx.rct.out, xx.svs, taut, z, condtype){
  n <- length(obs.rct)
  m <- nrow(xx.svs)
  
  # ordering the data according to u increasing:
  u.order <- order(obs.rct)
  u <- obs.rct[u.order]
  x.censor <- xx.rct.censor[u.order, ]
  x.out <- xx.rct.out[u.order, ]
  delta <- delta.rct[u.order]
  z <- z[u.order]
  ##################################
  ### use logistic model to fit  ###
  ### the propensity probability ###
  ##################################
  pp_est <- as.numeric(pred.score[1:n][u.order])
  if(condtype == "lognorm"){
    #########################################
    ### estimate the proportional hazard  ###
    ### model parameters for both failure ###
    ### time and cesoring time            ###
    #########################################
    uz <- u[z!=0]   # using data only when z=1
    deltaz <- delta[z!=0]
    b_est_all <- survreg(Surv(uz, deltaz) ~ x.out[z!=0,], dist = "lognormal")
    b_est <- b_est_all$coef
    b_scale <- b_est_all$scale
    bc_est <- coxph(Surv(uz, 1-deltaz) ~ x.censor[z!=0,])$coef   # fit cox ph model for censoring
    
    #######################################
    ### computing parametric AFT model  ###
    ### for failure time          ###
    ### and propotional hazards for censoring, denoted by       ###
    ### dls(j) and dlc(j) for time u(j) ###
    #######################################
    s2 <- cumsum(z[n:1]*exp(as.matrix(x.censor[n:1,])%*%bc_est))
    dlc <- z*(1-delta)/s2[n:1]
    dlc[(1:length(dlc))[is.nan(dlc)]] <- 0
    
    ###########################################
    ### for patient i with covariates x[i], ###
    ### compute estimated survival time and ###
    ### censoring at time point u[j]        ###
    ###########################################
    # sss(i, j) -> H(u_j, X_i); ssc -> K(u_j, X_i)
    sss <- matrix(0, n, n)
    res <- matrix(0, n, n)
    sss_svs <- matrix(0, m, n)
    res_svs <- matrix(0, m, n)
    ssc <- matrix(0, n, n)
    for(j in 1:n){
      res[,j] <- (log(u[j]) - as.matrix(cbind(1,x.out)) %*%b_est)/b_scale
      sss[,j] <- 1-pnorm(res[,j])		
    }
    for(j in 1:n){
      res_svs[,j] <- (log(u[j]) - as.matrix(cbind(1,xx.svs)) %*%b_est)/b_scale
      sss_svs[,j] <- 1-pnorm(res_svs[,j])		
    }
    ssc[, 1] <- dlc[1]*exp(as.matrix(x.censor)%*%bc_est)
    
    otc <- outer(c(exp(as.matrix(x.censor)%*%bc_est)), dlc, "*")   
    otc[, 1] <- rep(0, n)
    otc_cum <- t(apply(otc, 1, cumsum))
    ssc_1 <- outer(ssc[, 1], rep(1, n))
    ssc <- ssc_1 + otc_cum
    ssc <- exp(-ssc)
  }
  if(condtype == "cox"){
    #########################################
    ### estimate the proportional hazard  ###
    ### model parameters for both failure ###
    ### time and cesoring time            ###
    #########################################
    uz <- u[z!=0]   # using data only when z=1
    deltaz <- delta[z!=0]
    b_est <- coxph(Surv(uz, deltaz) ~ x.out[z!=0,])$coef   # fit cox ph model for surv time
    bc_est <- coxph(Surv(uz, 1-deltaz) ~ x.censor[z!=0,])$coef
    
    #######################################
    ### computing propotional hazards   ###
    ### model for both failure time     ###
    ### and censoring, denoted by       ###
    ### dls(j) and dlc(j) for time u(j) ###
    #######################################
    s1 <- cumsum(z[n:1]*exp(as.matrix(x.out[n:1,])%*%b_est))
    dls <- z*delta/s1[n:1]
    dls[(1:length(dls))[is.nan(dls)]] <- 0
    s2 <- cumsum(z[n:1]*exp(as.matrix(x.censor[n:1,])%*%bc_est))
    dlc <- z*(1-delta)/s2[n:1]
    dlc[(1:length(dlc))[is.nan(dlc)]] <- 0
    
    ###########################################
    ### for patient i with covariates x[i], ###
    ### compute estimated survival time and ###
    ### censoring at time point u[j]        ###
    ###########################################
    # sss(i, j) -> H(u_j, X_i); ssc -> K(u_j, X_i)
    sss <- matrix(0, n, n)
    sss_svs <- matrix(0, m, n)
    ssc <- matrix(0, n, n)
    # each column is one time point, row an individual
    sss[, 1] <- dls[1]*exp(as.matrix(x.out)%*%b_est)
    sss_svs[, 1] <- dls[1]*exp(as.matrix(xx.svs)%*%b_est)
    ssc[, 1] <- dlc[1]*exp(as.matrix(x.censor)%*%bc_est)
    
    ots <- outer(c(exp(as.matrix(x.out)%*%b_est)), dls, "*")   
    ots[, 1] <- rep(0, n)
    ots_cum <- t(apply(ots, 1, cumsum))
    sss_1 <- outer(sss[, 1], rep(1, n))
    sss <- sss_1 + ots_cum
    
    otsvs <- outer(c(exp(as.matrix(xx.svs)%*%b_est)), dls, "*")   
    otsvs[, 1] <- rep(0, m)
    otsvs_cum <- t(apply(otsvs, 1, cumsum))
    sss_1_svs <- outer(sss_svs[, 1], rep(1, n))
    sss_svs <- sss_1_svs + otsvs_cum
    
    otc <- outer(c(exp(as.matrix(x.censor)%*%bc_est)), dlc, "*")   
    otc[, 1] <- rep(0, n)
    otc_cum <- t(apply(otc, 1, cumsum))
    ssc_1 <- outer(ssc[, 1], rep(1, n))
    ssc <- ssc_1 + otc_cum
    sss <- exp(-sss)
    sss_svs <- exp(-sss_svs)
    ssc <- exp(-ssc)
  }
  ################################################
  ### compute the martingale integral part,    ###
  ### denote the augmentation term for patient ###
  ### i at time u[j] by q[i, j]                ###
  ################################################
  q <- matrix(0,n,n)
  ld <- lower.tri(matrix(1, n, n), diag=TRUE)   # lower triangular matrix
  ud <- matrix(1, n, n)-ld   # upper triangular matrix
  otr <- outer(c(exp(as.matrix(x.censor)%*%bc_est)), dlc)/(sss*ssc)  
  # lambda_c(r,X)/[K_c(r,X)*H(r,X)]
  q <- t(apply(otr, 1, cumsum)) # \int{lambda_c(r,X)/[K_c(r,X)*H(r,X)]}
  q <- q*ld+diag(q)*ud  # \int{Y(r)*lambda_c(r,X)/[K_c(r,X)*H(r,X)]}
  
  ###############################
  ### estimate survival dist  ###
  ### at each time point u[j] ###
  ###############################
  
  stfactor <- sum(z/pp_est)
  surv <- rep(0, n)
  s_if <- matrix(0, n+m, n)
  
  temp4 = (1/stfactor)*z*(1-delta)*sss/(diag(sss)*diag(ssc)*c(pp_est))   # Z*H/pi*\int{dN_c/[K*H]}     
  ou <- -(1/stfactor)*z*sss/pp_est - (1/stfactor)*z*sss*q/pp_est + temp4
  ol <- (1/stfactor)*z/(pp_est*ssc) - (1/stfactor)*z*sss/pp_est - (1/stfactor)*z*sss*q/pp_est
  
  temp <- ou*ud+ol*ld+diag(diag(temp4))   
  for(j in 1:n){
    tempj = temp[,j]
    tempj[is.na(tempj)] = 0
    temp[,j] = tempj
  }
  temp_svs <- sweight[-(1:n)]*sss_svs/svw
  temp_all <- rbind(temp,temp_svs)
  surv <- apply(temp_all, 2, sum)
  
  v1 = apply(temp,2,sum)
  v2 = apply(temp_svs,2,sum)
  
  mv1 = matrix(rep(v1,n),byrow=T,n,n) #n*n matrix
  s_if_rct = temp - (1/stfactor)*(z/pp_est)*mv1
  mv2 = matrix(rep(v2,m),byrow=T,m,n)  #m*n matrix
  s_if_svy = sweight[-(1:n)]*(sss_svs - mv2)/svw
  s_if = rbind(s_if_rct,s_if_svy)
  
  list <- list(surv=surv, s_if=s_if)
  return(list)
}
