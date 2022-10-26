#visualization for weak and srong sampling in Figure 1
rm(list = ls())

##case 1(weak sampling): gampar <- c(-6.5,-0.2,-0.1,-0.5)
visualizedata <- function(gampar, title, choice){
  #population size
  N = 10^6
  # covariates
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  x3 <- rbinom(N,1,0.5)
  
  # Sampling mechanism
  samp <-function(par,x){
    exp(x %*% par)/(1+exp(x %*% par))
  }
  
  xpop <- cbind(x1,x2,x3)
  p <- samp(gampar,cbind(1,xpop))
  S <- rbinom(N,1,p)
  
  M <-sum(S==0) # M = N-n
  xpopS0 <- xpop[which(S==0),] # N-n population
  # logistic model for D
  etapar <- c(-7,0.3,0.4,0.2)
  p0 <- samp(etapar,cbind(1,xpopS0))
  D <- rbinom(M,1,p0)
  
  idsamp = which(D==1)
  
  # proportion of ps outside range
  prop1 <- mean(p < 0.01)
  # prop2 <- mean(p < 0.1 | p > 0.9)
  prop2 <- mean(p > 0.95)
  prop3 <- mean(p0 < 0.01)
  prop4 <- mean(p0 > 0.95)
  print(c(prop1,prop2,prop3,prop4))
  if(choice==1){
    # histograms
    hist(p[S==1],breaks=50,col="gray77",border="gray77",main=title,xlim=c(0,0.005),xlab="true sampling scores",freq=F,
         cex.lab = 0.8, cex.axis = 0.8 ,cex.main = 0.8, cex = 0.8, ylim=c(0,1500))
    hist(p[S==0][idsamp],breaks=50,col="NA",add=TRUE,freq=F)
    legend("topright",legend=c("Trial participants","Non-participants"),col=c("gray77","black"),lty=1,lwd=1,bty='n',cex=0.8)
  }else{
    # histograms
    hist(p[S==1],breaks=50,col="gray77",border="gray77",main=title,xlim=c(0,0.04),xlab="true sampling scores",freq=F,
         cex.lab = 0.8, cex.axis = 0.8,cex.main = 0.8, cex = 0.8, ylim=c(0,1300))
    hist(p[S==0][idsamp],breaks=50,col="NA",add=TRUE,freq=F)
    legend("topright",legend=c("Trial participants","Non-participants"),col=c("gray77","black"),lty=1,lwd=1,bty='n',cex=0.8)
  }
}

gampar <- c(-6.5,-0.2,-0.1,-0.5)
visualizedata(gampar, title="weak sampling",choice=1)

gampar <- c(-6.7,-0.8,-0.4,-1)
visualizedata(gampar, title="strong sampling",choice=2)

