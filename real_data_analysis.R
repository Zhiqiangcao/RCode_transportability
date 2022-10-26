#### Real data analysis in Section 6
rm(list=ls())

#Set directory to the "R Code" folder
# setwd(".../R Code") # change working directory to the current directory

### Accord study
#out.dir="D:/A1 Methods Research/Causal Inference/Generalizability Survival/Application/"
#load(paste0(out.dir,"accord9.RData"))

### further cleaning data
accord9$ins_cover[is.na(accord9$ins_cover)]=0
accord9$x4smoke[accord9$x4smoke==2]=0
accord9$hartfail[accord9$hartfail==2]=0
accord9$x2mi[accord9$x2mi==2]=0
accord9$x2stroke[accord9$x2stroke==2]=0

accord9<-subset(accord9, select=-c(MaskID,treatment,dk_unins,wt_kg,ht_cm,
                                   x4bmi,Composite,CompositeDays))
accord9 <- accord9[complete.cases(accord9), ]
dim(accord9) # 4458, 35


# confirm trial statistics
library(table1)
table1(~ age + factor(female) + factor(race) + factor(ins_cover) +
         factor(x4smoke) + factor(edu) +  
         factor(hartfail) + factor(x2mi) + factor(x2stroke) + 
         yrsdiab + bmi + sbp + dbp + hdl + ldl + 
         trig + fpg + hba1c + gfr + uacr, data=accord9)


# install.packages('survminer')
library(survminer)
library(survival)
library(Rcpp)
library(ggplot2)

accord9$censor_po = abs(accord9$censor_po - 1)

###############################################################################################
### need to include IPCW here so that the trial estimates are robust under censoring at random.
###############################################################################################
fit.trial <- survfit(Surv(fuyrs_po, censor_po) ~ trt,
                     data = accord9)

#Zhiqiang: modify ltimes from 0 to 6.5
ltimes <- seq(0,7,by=0.1) 
K <- length(ltimes)

pred.data = summary(fit.trial,times=ltimes)


B=1000
EST1=EST0=matrix(NA,B,length(ltimes))
n=nrow(accord9)
set.seed(20220209)
for(b in 1:B){
  id = sample(1:n, n, replace = TRUE)
  accord9.boot = accord9[id,]
  fit.boot = survfit(Surv(fuyrs_po, censor_po) ~ trt,
                     data = accord9.boot)
  pred.boot = summary(fit.boot,times=ltimes)
  EST1[b,]<-pred.boot$surv[1:length(ltimes)]
  EST0[b,]<-pred.boot$surv[(length(ltimes)+1):(2*length(ltimes))]
}

## plot trial sample results
est=pred.data$surv
ub=pred.data$upper
lb=pred.data$lower

#rename the following data.frame as dat_rct (original is dat)
dat_rct=data.frame(ltimes=ltimes,est=est, ub=ub, lb=lb,
                   group=rep(c("Intensive BP", "Standard BP"),each=length(ltimes)))

fig1=ggplot(dat_rct,aes(x = ltimes,y = est, colour=group)) +
  geom_smooth(aes(group=group,ymin = lb, ymax = ub),stat = "identity",fill=rep(c(2,5),each=length(ltimes))) +
  xlab("Time (t) in years") + ylab("Survival probability") + 
  ggtitle("Counterfactual survival functions (trial)") + 
  theme(plot.title = element_text(size = 15, face = "bold"),
        axis.text=element_text(size=15),
        axis.title = element_text(size = 15),
        legend.position=c(0.2,0.2),legend.text = element_text(size=12),
        legend.title = element_text(size=13,face="bold"),
        legend.background = element_rect(fill="transparent", 
                                         size=0.5, linetype="solid")) + 
  scale_x_continuous(breaks = round(seq(0, 7, by = 1),1))+ylim(0.8,1)

fig1

# plot difference
diff=pred.data$surv[1:length(ltimes)]-pred.data$surv[(length(ltimes)+1):(2*length(ltimes))]
ubdiff=apply(EST1-EST0,2,function(x){quantile(x,0.975)})
lbdiff=apply(EST1-EST0,2,function(x){quantile(x,0.025)})

#rename the following data.frame as diff_data_rct (original is diff_dat)
diff_dat_rct <- data.frame(ltimes=ltimes,diff=diff, ubdiff=ubdiff, lbdiff=lbdiff)

fig2=ggplot(diff_dat_rct,aes(x = ltimes,y = diff)) +
  geom_smooth(aes(ymin = ubdiff, ymax = lbdiff),stat = "identity") +
  xlab("Time (t) in years") + ylab("Survival difference") + ggtitle("Sample average survival treatment effect") + 
  theme(plot.title = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 15),
        axis.text=element_text(size=15)) + 
  scale_x_continuous(breaks = round(seq(0, 7, by = 1),1)) + ylim(-0.06, 0.0777)
fig2 

library(ggplot2)
library(reshape2)
library(plyr)
library(gridExtra)
library(grid)
library(lattice)
library(Matrix)

uni = grid.arrange(arrangeGrob(fig1, fig2, nrow=1, ncol=2))

setwd(out.dir)
ggsave("Trial_results_new.pdf", plot=uni, width= 28, height=13, units="cm", dpi=600)

###################################################################
#### now estimate the target population average treatment effect
###################################################################
# install.packages(c("haven", "sas7bdat"))
library(sas7bdat)

nhanes = read.sas7bdat(paste0(out.dir,"nhanes_analysis_10yr_revisions.sas7bdat"))
nhanes_sub1 = subset(nhanes, include==1,drop=T)

nhanes9<-data.frame(age=nhanes_sub1$RIDAGEYR,
                    female=nhanes_sub1$female,
                    race=nhanes_sub1$race,
                    ins_cover=nhanes_sub1$insured,
                    x4smoke=as.numeric(nhanes_sub1$smoker==1),
                    edu=nhanes_sub1$edu_4cat,
                    hartfail=nhanes_sub1$chf,
                    x2mi=nhanes_sub1$mi,
                    x2stroke=nhanes_sub1$stroke,
                    yrsdiab=nhanes_sub1$yrsdiab,
                    bmi=nhanes_sub1$bmi,
                    sbp=nhanes_sub1$sbp,
                    dbp=nhanes_sub1$dbp,
                    hdl=nhanes_sub1$hdl,
                    ldl=nhanes_sub1$ldl,
                    trig=nhanes_sub1$trig,
                    fpg=nhanes_sub1$fpg,
                    hba1c=nhanes_sub1$hba1c,
                    gfr=nhanes_sub1$gfr,
                    uacr=nhanes_sub1$uacr,
                    weight=nhanes_sub1$new_weight)
dim(nhanes9) # 4289   21
nhanes9 = nhanes9[complete.cases(nhanes9$weight), ]
nhanes9 = nhanes9[nhanes9$weight > 0,]
dim(nhanes9) # 1943 21


table1(~ age + factor(female) + factor(race) + factor(ins_cover) +
         factor(x4smoke) + factor(edu) +  
         factor(hartfail) + factor(x2mi) + factor(x2stroke) + 
         yrsdiab + bmi + sbp + dbp + hdl + ldl + 
         trig + fpg + hba1c + gfr + uacr + weight, data=nhanes9)

nhanes9 <- nhanes9[complete.cases(nhanes9), ]
dim(nhanes9) # 1647, 21

nhanes9 <- nhanes9[nhanes9$dbp >= 30,]
dim(nhanes9) # 1639, 21

nhanes9=nhanes9[nhanes9$uacr <= 2000,]
dim(nhanes9) # 1619   21


table1(~ age + factor(female) + factor(race) + factor(ins_cover) +
         factor(x4smoke) + factor(edu) +  
         factor(hartfail) + factor(x2mi) + factor(x2stroke) + 
         yrsdiab + bmi + sbp + dbp + hdl + ldl + 
         trig + fpg + hba1c + gfr + uacr + weight, data=nhanes9)

## combine data sets
library(survival)
datrct <- data.frame(obs=accord9$fuyrs_po,delta=accord9$censor_po,
                     z=as.numeric(accord9$trt=="Intensive BP"))
varlist <- c("age","female","race","ins_cover","x4smoke",
             "edu","hartfail","x2mi","x2stroke","yrsdiab",
             "bmi","sbp","dbp","hdl","ldl","trig","fpg",
             "hba1c","gfr","uacr")
datrct <- cbind(datrct, accord9[,varlist])


# estimation of sampling score
n = nrow(datrct) #4458
m = nrow(nhanes9) #1619
resp <- c(rep(1,n),rep(0,m))
datrct$race = as.numeric(datrct$race)
xall <- rbind(datrct[,-c(1:3)],nhanes9[,varlist])
sweight <- c(rep(1,n),nhanes9$weight)
dat.cb <- data.frame(resp=resp,xall,sweight=sweight)
id <- 1:(dim(dat.cb)[1])

dat.cb$age=(dat.cb$age - mean(dat.cb$age))/sd(dat.cb$age)
dat.cb$yrsdiab=(dat.cb$yrsdiab - mean(dat.cb$yrsdiab))/sd(dat.cb$yrsdiab)
dat.cb$bmi=(dat.cb$bmi - mean(dat.cb$bmi))/sd(dat.cb$bmi)
dat.cb$sbp=(dat.cb$sbp - mean(dat.cb$sbp))/sd(dat.cb$sbp)
dat.cb$dbp=(dat.cb$dbp - mean(dat.cb$dbp))/sd(dat.cb$dbp)
dat.cb$hdl=(dat.cb$hdl - mean(dat.cb$hdl))/sd(dat.cb$hdl)
dat.cb$ldl=(dat.cb$ldl - mean(dat.cb$ldl))/sd(dat.cb$ldl)
dat.cb$trig=(dat.cb$trig - mean(dat.cb$trig))/sd(dat.cb$trig)
dat.cb$fpg=(dat.cb$fpg - mean(dat.cb$fpg))/sd(dat.cb$fpg)
dat.cb$hba1c=(dat.cb$hba1c - mean(dat.cb$hba1c))/sd(dat.cb$hba1c)
dat.cb$gfr=(dat.cb$gfr - mean(dat.cb$gfr))/sd(dat.cb$gfr)
dat.cb$uacr=(dat.cb$uacr - mean(dat.cb$uacr))/sd(dat.cb$uacr)

formula.ss <- as.formula(resp ~ age + female + factor(race) + ins_cover +
                           x4smoke + factor(edu) +  
                           hartfail + x2mi + x2stroke + 
                           yrsdiab + bmi + sbp + dbp + hdl + ldl + 
                           poly(trig,2) + fpg + hba1c + gfr + uacr + 
                           factor(race):bmi + yrsdiab:sbp)

fit.score <- glm(formula.ss, weights = sweight,data=dat.cb, family = binomial)
pred.score.tmp <- predict(fit.score,type="response")
clp<-as.numeric(quantile(pred.score.tmp,c(0.1)))
cup<-as.numeric(quantile(pred.score.tmp,c(1))) 
pred.score.tmp[pred.score.tmp<=clp]=clp
pred.score.tmp[pred.score.tmp>=cup]=cup


#survey weights, PS weights
svw <- sum(sweight[(n+1):(n+m)])
pred.score <- pred.score.tmp/(1-pred.score.tmp)
w.h <- pred.score[1:n]

#we used all covariates in censoring model
censor.fit = coxph(Surv(datrct$obs, 1-datrct$delta) ~ age + female + factor(race) + ins_cover +
                     x4smoke + factor(edu) + hartfail + x2mi + x2stroke + 
                     yrsdiab + bmi + sbp + dbp + hdl + ldl + trig + fpg + hba1c + gfr + uacr, 
                   data=datrct)


## DR estimation
xx.svs <- xall[resp==0,]
sweight <- as.numeric(sweight)

# Estimation of conditional expectation
setwd(out.dir)
source("sim_code_taste_final.R")

#we include all covariates in the outcome model
out.fit = coxph(Surv(datrct$obs, datrct$delta) ~ age + female + factor(race) + ins_cover +
                  x4smoke + factor(edu) + hartfail + x2mi + x2stroke + 
                  yrsdiab + bmi + sbp + dbp + hdl + ldl + trig + fpg + hba1c + gfr + uacr, 
                data=datrct)

xx.rct.out = model.matrix(out.fit)
xx.svs = model.matrix(out.fit, data=xall[(n+1):(n+m),])

#DR2
dr2.ATE <- function(obs.rct, delta.rct, z.rct, xx.rct.censor, xx.rct.out, xx.svs, taut, condtype){
  # treatment-specific counterfactual
  
  # estimate fot z=1:
  r_z1 <- survres_z2app(obs.rct, delta.rct, z.rct, xx.rct.censor, xx.rct.out, xx.svs, taut, z.rct, condtype)
  S_hat_z1 <- r_z1$surv   # s is used to estimate the S_hat (survival function estimator)
  sif_z1 <- r_z1$s_if   # ssq is used to estimate the v_hat (variance of S_hat)
  
  # estimate fot z=0:
  r_z0 <- survres_z2app(obs.rct, delta.rct, z.rct, xx.rct.censor, xx.rct.out, xx.svs, taut, 1-z.rct, condtype)
  S_hat_z0 <- r_z0$surv   # s is used to estimate the S_hat (survival function estimator)
  sif_z0 <- r_z0$s_if
  
  # return results
  u.order <- order(obs.rct)
  u <- obs.rct[u.order]
  est <- se <- NULL
  S1 <- S0 <- S1.se <- S0.se <- NULL
  for(tt in taut){
    id = sum(u<=tt)
    if(id==0){
      est <- c(est,0)
      S1 <- c(S1, 1)
      S0 <- c(S0, 1)
    } else {
      est <- c(est, S_hat_z1[id]-S_hat_z0[id])
      S1 <- c(S1, S_hat_z1[id])
      S0 <- c(S0, S_hat_z0[id])
    }
    se <- c(se, sqrt(sum((sif_z1[,id])^2) + sum((sif_z0[,id])^2)))
    S1.se <- c(S1.se, sqrt(sum((sif_z1[,id])^2)))
    S0.se <- c(S0.se, sqrt(sum((sif_z0[,id])^2)))
  }
  return(list(est=est, se=se, S1=S1, S0=S0, S1.se=S1.se, S0.se=S0.se))
}

xx.rct.out = model.matrix(out.fit)
xx.svs = model.matrix(out.fit, data=xall[(n+1):(n+m),])

dr2.est.ATE.cox <- dr2.ATE(obs.rct=datrct$obs,delta.rct=datrct$delta,z.rct=datrct$z,
                           xx.rct.censor=model.matrix(censor.fit),
                           xx.rct.out = xx.rct.out,
                           xx.svs=xx.svs,
                           taut=ltimes,condtype="cox")

#############################################
# Plot DR2 Cox
#############################################
## plot trial sample results
est=c(dr2.est.ATE.cox$S1, dr2.est.ATE.cox$S0)
ub=pmin(est+qnorm(0.975)*c(dr2.est.ATE.cox$S1.se, dr2.est.ATE.cox$S0.se),1)
lb=pmax(est-qnorm(0.975)*c(dr2.est.ATE.cox$S1.se, dr2.est.ATE.cox$S0.se),0)
#Zhiqiang: rename the following data.frame as dat_dr2 (original is dat)
dat_dr2=data.frame(ltimes=ltimes,est=est, ub=ub, lb=lb,
                   group=rep(c("Intensive BP", "Standard BP"),each=length(ltimes)))


fig_pop_dr2=ggplot(dat_dr2,aes(x = ltimes,y = est, colour=group)) +
  geom_smooth(aes(group=group,ymin = lb, ymax = ub),stat = "identity",fill=rep(c(2,5),each=length(ltimes))) +
  xlab("Time (t) in years") + ylab("Survival probability") + 
  ggtitle("Counterfactual survival functions (population)") + 
  theme(plot.title = element_text(size = 15, face = "bold"),
        axis.text=element_text(size=15),
        axis.title = element_text(size = 15),
        legend.position=c(0.2,0.2),legend.text = element_text(size=12),
        legend.title = element_text(size=13,face="bold"),
        legend.background = element_rect(fill="transparent", 
                                         size=0.5, linetype="solid")) + 
  scale_x_continuous(breaks = round(seq(0, 7, by = 1),1)) + ylim(0.8,1)

fig_pop_dr2



##### initial plotting on survival differences
diff=dr2.est.ATE.cox[[1]]
ubdiff=dr2.est.ATE.cox[[1]]+qnorm(0.975)*dr2.est.ATE.cox[[2]]
lbdiff=dr2.est.ATE.cox[[1]]-qnorm(0.975)*dr2.est.ATE.cox[[2]]
#Zhiqiang: rename the following data.frame as diff_dat_dr2 (original is diff_dat)
diff_dat_dr2 <- data.frame(ltimes=ltimes,diff=diff, ubdiff=ubdiff, lbdiff=lbdiff)


fig_diff_dr2=ggplot(diff_dat_dr2,aes(x = ltimes,y = diff)) +
  geom_smooth(aes(ymin = ubdiff, ymax = lbdiff),stat = "identity",method="loess") +
  xlab("Time (t) in years") + ylab("Survival difference") + ggtitle("Target average survival treatment effect") + 
  theme(plot.title = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 15),
        axis.text=element_text(size=15)) + 
  scale_x_continuous(breaks = round(seq(0, 7, by = 1),1)) + ylim(-0.06, 0.0777)
fig_diff_dr2

uni.DR2 = grid.arrange(arrangeGrob(fig_pop_dr2, fig_diff_dr2, nrow=1, ncol=2))

setwd(out.dir)
ggsave("TASTE_DR2_results_new.pdf", plot=uni.DR2, width= 28, height=13, units="cm", dpi=600)


##calcualte results in Table 4 of manuscript
ltimes
index = seq(6,71,by=5) #interested target survival time
ltimes[index]
diff_surv_rct = diff_dat_rct[index,c(1,2,4,3)]
diff_surv_dr2 = diff_dat_dr2[index,c(1,2,4,3)]

#library(xtable)
#diff_surv_cb = cbind(diff_surv_rct,diff_surv_dr2)
#diff_surv_cb = diff_surv_cb[,-5]
#xtable(diff_surv_cb,digits=5)
