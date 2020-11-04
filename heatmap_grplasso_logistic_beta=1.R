library(glmnet)
library(dummies)
library(grplasso)
setwd(getwd())
source("functions_refit_bic.R")
set.seed(12345)
outputPath="./Heatmap_grplasso_logistic_beta=1"
if(!dir.exists(outputPath)){
  dir.create(outputPath)
}

p=15  ### number of covariates
n=20000
true.mean=0
true.sd=1
beta=1
tau=0.8
index=c(NA,rep(1,3),rep(2,4),rep(3,3),rep(4,3),rep(5,2))
par=c(0,rep(beta,3),rep(beta,4),rep(0,3),rep(beta,3),rep(0,2)) 


data=data_glm_gr(n,p,true.mean,true.sd,par)

basis=read.csv(paste(outputPath,"/beta_hat.csv",sep=""))
not_zero=apply(basis,2,function(x){sum(x!=0)})/nrow(basis)>0.8
sd0=apply(basis,2,sd)
c0=apply(basis,2,function(x){quantile(x,0.975)-quantile(x,0.025)})

##################################################################################
#################################################
r_vec=c(5,10,30,50,100,200,500)
s_vec=c(3,5,10,25,50,100)
sd_rel=matrix(0,length(r_vec),length(s_vec))
ci_rel=matrix(0,length(r_vec),length(s_vec))

for(i in 1:length(r_vec)){
  for(j in 6:6){
    temp=BLB_grp(data,0.7,index,r_vec[i],tau,s_vec[j],family="binomial")
    sd_rel[i,j]=rel.err0(temp$std_avg,sd0,not_zero)[s_vec[j]]
    ci_rel[i,j]=rel.err0(temp$c_avg,c0,not_zero)[s_vec[j]]
    }
}


write.csv(sd_rel,paste(outputPath,"/heat.mapSD.csv",sep=""),row.names=F)
write.csv(ci_rel,paste(outputPath,"/heat.mapCI.csv",sep=""),row.names=F)






