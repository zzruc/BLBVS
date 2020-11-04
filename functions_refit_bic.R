quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}


########################################################################################
######################### Lasso ################################################


data_lm_la<-function(n,p,true.mean,true.sd,par){
  continue=matrix(0,n,p)
  for(i in 1:p){
    continue[,i]=rnorm(n,true.mean,true.sd)
  }
  scale_x=scale(continue)
  scale_x=cbind(1,scale_x)
  colnames(scale_x)=c("Intercept",paste("X",1:p,sep=""))
  
  epsilon=rnorm(n,mean=0,sd=1)
  y=scale_x%*%par+epsilon
  
  test.data=as.matrix(data.frame(scale_x[,-1],y))
  return(test.data)
}

data_glm_la<-function(n,p,true.mean,true.sd,par){
  continue=matrix(0,n,p)
  for(i in 1:p){
    continue[,i]=rnorm(n,true.mean,true.sd)
  }
  scale_x=scale(continue)
  scale_x=cbind(1,scale_x)
  colnames(scale_x)=c("Intercept",paste("X",1:p,sep=""))
  
  prob=round(1/(1+as.vector(exp(-scale_x%*%par))),digits=3)
  y=c()
  for(i in 1:n){
    y[i]=rbinom(1,size=1,prob=prob[i])
  }
  
  test.data=as.matrix(data.frame(scale_x[,-1],y))
  return(test.data)
}

appr_sd_la<-function(n,iter.num,p,true.mean,true.sd,par,family){
  coef.simu=matrix(0,nrow=iter.num,ncol=p) 
  for(i in 1:iter.num){
    set.seed(i)
    print(i)
    if(family=="gaussian"){
      data=data_lm_la(n,p,true.mean,true.sd,par)
    }
    if(family=="binomial"){
      data=data_glm_la(n,p,true.mean,true.sd,par)
    }
    
    fit=glmnet(x=data[,-(p+1)],y=data[,(p+1)],family=family)
    lambda.max=max(fit$lambda)  
    lambda=c(lambda.max*0.96^seq(0,100,1))
    fit_train=quiet(glmnet(x=data[,-(p+1)],y=data[,(p+1)]
                           ,family=family,lambda=lambda))
    
    bic=rep(Inf,length(lambda))
    for(ii in 1:length(lambda)){
      temp.coef=as.vector(coef(fit_train)[,ii])[-1]
      remove.ind=which(temp.coef==0)
      if(length(remove.ind)<(p-1)){
        fit=glmnet(data[,-c(remove.ind,p+1)],data[,(p+1)],
                   family=family,alpha=0,lambda=10^(-5))
        tLL=fit$nulldev-deviance(fit)
        k=fit$df
        nn=fit$nobs
        bic[ii]=log(nn)*k-tLL
      }
    }
    min.ind=which.min(bic)
    
    coef.simu[i,]=as.vector(coef(fit_train)[,min.ind])[-1]
    remove.ind=which(coef.simu[i,]==0)
    fit=glmnet(data[,-c(remove.ind,p+1)],data[,(p+1)],
            family=family,alpha=0,lambda=10^(-5))
    
    if(length(remove.ind)==0){
      coef.simu[i,]=fit$beta[,1]
    }else{
      coef.simu[i,-remove.ind]=fit$beta[,1]
    }
  }
  beta.sd=apply(coef.simu,2,sd)
  c0=apply(coef.simu,2,function(x){quantile(x,0.975)-quantile(x,0.025)})
  res=list("beta"=coef.simu,"beta.sd"=beta.sd,"c0"=c0)
  return(res)  
}


BLB_la<-function(data,gamma,iter.num,tau,s,family){
    
  p=ncol(data)-1
  n=nrow(data)
  b=floor(n^gamma)

  std.bag=matrix(0,nrow=s,ncol=p)
  CI.bag=matrix(0,nrow=s,ncol=p)
  time=vector()
  var.instances=matrix(0,nrow=s*iter.num,ncol=p)
  colnames(var.instances)=c(paste("X",1:p,sep=""))
  #no sd for the first time
  var.instances[1,]=0
  beta_hat=NULL
  coef=matrix(0,nrow=iter.num,ncol=p)
  
  # loop for every bag
  for(i in 1:s) {  
    time[i]=system.time({
      bag.data=data[sample(c(1:n),b),]     
      for(j in 1:iter.num) {
        print(paste("i=",i,",j=",j,sep=""))
        boot.ind=sample(c(1:nrow(bag.data)),nrow(data),replace=T)
        boot.ind.uni=sort(unique(boot.ind))
        boot.weights=as.vector(table(boot.ind))
        
        fit=glmnet(bag.data[boot.ind.uni,-(p+1)],
                   bag.data[boot.ind.uni,(p+1)],
                   weights=boot.weights,family=family)
        lambda.max=max(fit$lambda)  
        lambda=c(lambda.max*0.96^seq(0,100,1))
        fit_train=quiet(glmnet(bag.data[boot.ind.uni,-(p+1)],
                               bag.data[boot.ind.uni,(p+1)],
                               weights=boot.weights,
                               family=family,lambda=lambda))
        
        bic=rep(Inf,length(lambda))
        for(ii in 1:length(lambda)){
          temp.coef=as.vector(coef(fit_train)[,ii])[-1]
          remove.ind=which(temp.coef==0)
          if(length(remove.ind)<(p-1)){
            fit=glmnet(bag.data[boot.ind.uni,-c(remove.ind,p+1)],
                       bag.data[boot.ind.uni,p+1],
                       weights=boot.weights,family=family,
                       alpha=0,lambda=10^(-5))
            tLL=fit$nulldev-deviance(fit)
            k=fit$df
            nn=fit$nobs
            bic[ii]=log(nn)*k-tLL
          }
        }
        min.ind=which.min(bic)
        
        coef[j,]=as.vector(coef(fit_train)[,min.ind])[-1]
        remove.ind=which(coef[j,]==0)
        fit=glmnet(bag.data[boot.ind.uni,-c(remove.ind,p+1)],
                   bag.data[boot.ind.uni,p+1],
                   weights=boot.weights,family=family,
                   alpha=0,lambda=10^(-5))
        
        if(length(remove.ind)==0){
          coef[j,]=fit$beta[,1]
        }else{
          coef[j,-remove.ind]=fit$beta[,1]
        }
        inclusion.index=which(coef[j,]!=0)
        # add 0 to turn bool value into int
        var.instances[j+(i-1)*iter.num,]=0+(colnames(data[,-(p+1)])%in%colnames(bag.data[,inclusion.index]))    
      }
      
      beta_hat=rbind(beta_hat,coef)
      
      std.bag[i,]=apply(coef,2,sd)
      CI.bag[i,]=apply(coef,2,function(x){quantile(x,0.975)-quantile(x,0.025)})
    })[1]
  } 
  # slected variable 
  sele.var=apply(var.instances,2,sum)
  sele.ratio=sele.var/(s*iter.num)  
  selected=which(sele.ratio>tau)   
  
  time_cum=cumsum(time)  
  
  ###### Standard Devision
  std_cum=apply(std.bag,2,cumsum)
  std_avg=std_cum
  for(i in 1:nrow(std_cum)){
    std_avg[i,]=std_cum[i,]/rep(i,ncol(std_cum))
  }
  
  ###### length of CI
  CI_cum=apply(CI.bag,2,cumsum)
  CI_avg=CI_cum
  for(i in 1:nrow(std_cum)){
    CI_avg[i,]=CI_cum[i,]/rep(i,ncol(CI_cum))
  }
  ## output
  res=list("sele.var"=sele.var,"sele.ratio"=sele.ratio,"selected"=selected,
           "std"=std.bag,"std_avg"=std_avg,
           "c"=CI.bag,"c_avg"=CI_avg,
           "time"=time,"tim_cum"=time_cum,"beta_hat"=beta_hat)
  return(res)
}


BOOT_la<-function(data,iter.num,tau,family){
  p=ncol(data)-1  
  coef.simu=matrix(0,iter.num,p);colnames(coef.simu)=paste("X",1:p,sep="")
  var.instances=matrix(0,nrow=iter.num,ncol=p);colnames(var.instances)=c(paste("X",1:p,sep=""))
  time=vector()
  for(i in 1:iter.num) {
    print(i)  
    boot.data=data[sample(1:nrow(data),nrow(data),replace=T),]
    time[i]=system.time({
      fit=glmnet(x=boot.data[,-(p+1)],y=boot.data[,(p+1)],
                 family=family)  
      lambda.max=max(fit$lambda)
      lambda=c(lambda.max*0.96^seq(0,100,1))
      fit_train=quiet(glmnet(x=boot.data[,-(p+1)],
                             y=boot.data[,(p+1)],
                             family=family,lambda=lambda))
        
      bic=rep(Inf,length(lambda))
      for(ii in 1:length(lambda)){
        temp.coef=as.vector(coef(fit_train)[,ii])[-1]
        remove.ind=which(temp.coef==0)
        if(length(remove.ind)<(p-1)){
          fit=glmnet(boot.data[,-c(remove.ind,p+1)],
                     boot.data[,(p+1)],
                     family=family,alpha=0,lambda=10^(-5))
          tLL=fit$nulldev-deviance(fit)
          k=fit$df
          nn=fit$nobs
          bic[ii]=log(nn)*k-tLL
        }
      }
      min.ind=which.min(bic)
      
      coef.simu[i,]=as.vector(coef(fit_train)[,min.ind])[-1]
      remove.ind=which(coef.simu[i,]==0)
      fit=glmnet(boot.data[,-c(remove.ind,p+1)],
                 boot.data[,(p+1)],
                 family=family,alpha=0,lambda=10^(-5))
      	
      if(length(remove.ind)==0){
          coef.simu[i,]=fit$beta[,1]
        }else{
          coef.simu[i,-remove.ind]=fit$beta[,1]
        }
    })[1]
    inclusion.index=which(coef.simu[i,]!=0)
    var.instances[i,]=0+(colnames(data[,-(p+1)])%in%colnames(data[,inclusion.index]))
  }
  # slected variable : vote
  sele.var=apply(var.instances,2,sum)
  std.bag=matrix(0,iter.num,p)
  for(i in 2:nrow(coef.simu)){
    for(j in 1:ncol(coef.simu)){
      std.bag[i,j]=sd(coef.simu[1:i,j])
    }
  }
  
  CI.bag=matrix(0,iter.num,p)
  for(i in 2:nrow(coef.simu)){
    for(j in 1:ncol(coef.simu)){
      CI.bag[i,j]=quantile(coef.simu[1:i,j],0.975)-quantile(coef.simu[1:i,j],0.025)
    }
  }
  
  sele.ratio=sele.var/iter.num  
  selected=which(sele.ratio>tau)   

  time_cum=cumsum(time)   
  
  
  
  res=list("sele.var"=sele.var,"sele.ratio"=sele.ratio,"selected"=selected,
            "std_avg"=std.bag,"c_avg"=CI.bag,
            "time"=time,"tim_cum"=time_cum,"beta_hat"=coef.simu)

  return(res)
}


#########################################################
###########Grouplasso#########################
data_lm_gr<-function(n,p,true.mean,true.sd,par){
  continue=matrix(0,n,(p-5))
  cat_var1=sample(c("A","B","C","D"),size=n,replace=TRUE,prob=rep(1/4,4))
  cat_var0=sample(c("E","F","G"),size=n,replace=TRUE,prob=rep(1/3,3))
  for(i in 1:(p-5)){
    continue[,i]=rnorm(n,true.mean,true.sd)
  }
  x=cbind(continue,dummy(cat_var1)[,-1],dummy(cat_var0)[,-1])
  scale_x=scale(x)
  scale_x=cbind(1,scale_x)
  colnames(scale_x)=c("Intercept",paste("X",1:p,sep=""))
  
  epsilon=rnorm(n,mean=0,sd=1)
  y=scale_x%*%par+epsilon
  
  test.data=as.matrix(data.frame(scale_x,y))
  return(test.data)
}


data_glm_gr<-function(n,p,true.mean,true.sd,par){
  continue=matrix(0,n,(p-5))
  cat_var1=sample(c("A","B","C","D"),size=n,replace=TRUE,prob=rep(1/4,4))
  cat_var0=sample(c("E","F","G"),size=n,replace=TRUE,prob=rep(1/3,3))
  for(i in 1:(p-5)){
    continue[,i]=rnorm(n,true.mean,true.sd)
  }
  x=cbind(continue,dummy(cat_var1)[,-1],dummy(cat_var0)[,-1])
  scale_x=scale(x)
  scale_x=cbind(1,scale_x)
  colnames(scale_x)=c("Intercept",paste("X",1:p,sep=""))
  
  prob=round(1/(1+as.vector(exp(-scale_x%*%par))),digits=3)
  y=c()
  for(i in 1:n){
    y[i]=rbinom(1,size=1,prob=prob[i])
  }
  
  test.data=as.matrix(data.frame(scale_x,y))
  return(test.data)
}


appr_sd_gr<-function(n,iter.num,p,true.mean,true.sd,par,index,family){
  
  coef.simu=matrix(0,nrow=iter.num,ncol=p)
  for(i in 1:iter.num){
    print(i)
    set.seed(i)
    if(family=="gaussian"){
      data=data_lm_gr(n,p,true.mean,true.sd,par)
      lambda.max=lambdamax(x=as.matrix(data[,-(p+2)]),
                           y=as.matrix(data[,(p+2)]),
                           index=index,penscale=sqrt,
                           model=LinReg())
      lambda=c(lambda.max*0.96^seq(0,100,1))
      fit_train=quiet(grplasso(x=as.matrix(data[,-(p+2)]),
                               y=as.matrix(data[,(p+2)]),
                               index=index,lambda=lambda,
                               model=LinReg(),penscale=sqrt))
    }
    if(family=="binomial"){
      data=data_glm_gr(n,p,true.mean,true.sd,par)
      lambda.max=lambdamax(x=as.matrix(data[,-(p+2)]),
                           y=as.matrix(data[,(p+2)]),
                           index=index,penscale=sqrt,
                           model=LogReg())
      lambda=c(lambda.max*0.96^seq(0,100,1))
      fit_train=quiet(grplasso(x=as.matrix(data[,-(p+2)]),
                               y=as.matrix(data[,(p+2)]),
                               index=index,lambda=lambda,
                               model=LogReg(),penscale=sqrt))
    }
    
    bic=rep(Inf,length(lambda))
    for(ii in 1:length(lambda)){
      temp.coef=as.vector(coef(fit_train)[,ii])[-1]
      remove.ind=(which(temp.coef==0)+1)
      if(length(remove.ind)<(p-1)){
        fit=glmnet(data[,-c(1,remove.ind,(p+2))],data[,(p+2)],
                   family=family,lambda=10^(-5))
        tLL=fit$nulldev-deviance(fit)
        k=fit$df
        nn=fit$nobs
        bic[ii]=log(nn)*k-tLL
      }
    }
    min.ind=which.min(bic)

    coef.simu[i,]=as.vector(coef(fit_train)[,min.ind])[-1]
    remove.ind=(which(coef.simu[i,]==0)+1)
    fit=glmnet(data[,-c(1,remove.ind,(p+2))],data[,(p+2)],
               family=family,lambda=10^(-5))
    if(length(remove.ind)==0){
      coef.simu[i,]=fit$beta[,1]
    }else{
      coef.simu[i,-(remove.ind-1)]=fit$beta[,1]
    }
  }
  beta.sd=apply(coef.simu,2,sd)
  c0=apply(coef.simu,2,function(x){quantile(x,0.975)-quantile(x,0.025)})
  res=list("beta"=coef.simu,"beta.sd"=beta.sd,"c0"=c0)
  return(res)
}


BLB_grp<-function(data,gamma,index,iter.num,tau,s,family){
  
  
  beta_hat=NULL
  p=ncol(data)-2
  n=nrow(data)
  b=floor(n^gamma)
  
  std.bag=matrix(0,nrow=s,ncol=p)
  CI.bag=matrix(0,nrow=s,ncol=p)
  time=vector()
  var.instances=matrix(0,nrow=s*iter.num,ncol=p)
  
  coef=matrix(0,nrow=iter.num,ncol=p)
  
  for(i in 1:s) {
    time[i]=system.time({
      bag.data=data[sample(c(1:n),b),]
      for(j in 1:iter.num) {
        print(paste("i=",i,",j=",j,sep=""))
        boot.ind=sample(c(1:nrow(bag.data)),nrow(data),replace=T)
        boot.ind.uni=sort(unique(boot.ind))
        boot.weight=as.vector(table(boot.ind))
        
        if(family=="gaussian"){
          lambda.max=lambdamax(x=as.matrix(bag.data[boot.ind.uni,-(p+2)]),
                               y=as.matrix(bag.data[boot.ind.uni,(p+2)]),
                               weights=boot.weight,
                               index=index,penscale=sqrt,
                               model=LinReg())
          lambda=c(lambda.max*0.96^seq(0,100,1))
          fit_train=quiet(grplasso(x=as.matrix(bag.data[boot.ind.uni,-(p+2)]),
                                   y=as.matrix(bag.data[boot.ind.uni,(p+2)]),
                                   weights=boot.weight,
                                   index=index,lambda=lambda,
                                   model=LinReg(),penscale=sqrt)) 
        }
        if(family=="binomial"){
          lambda.max=lambdamax(x=as.matrix(bag.data[boot.ind.uni,-(p+2)]),
                               y=as.matrix(bag.data[boot.ind.uni,(p+2)]),
                               weights=boot.weight,
                               index=index,penscale=sqrt,
                               model=LogReg())
          lambda=c(lambda.max*0.96^seq(0,100,1))
          fit_train=quiet(grplasso(x=as.matrix(bag.data[boot.ind.uni,-(p+2)]),
                                   y=as.matrix(bag.data[boot.ind.uni,(p+2)]),
                                   weights=boot.weight,
                                   index=index,lambda=lambda,
                                   model=LogReg(),penscale=sqrt))
        }
        
        bic=rep(Inf,length(lambda))
        for(ii in 1:length(lambda)){
          temp.coef=as.vector(coef(fit_train)[,ii])[-1]
          remove.ind=(which(temp.coef==0)+1)
          if(length(remove.ind)<(p-1)){
            fit=glmnet(bag.data[boot.ind.uni,-c(1,remove.ind,(p+2))],
                       bag.data[boot.ind.uni,(p+2)],
                       family=family,weights=boot.weight,
                       lambda=10^(-5))
            tLL=fit$nulldev-deviance(fit)
            k=fit$df
            nn=fit$nobs
            bic[ii]=log(nn)*k-tLL
          }
        }
        min.ind=which.min(bic)
        
        coef[j,]=as.vector(coef(fit_train)[,min.ind])[-1]
        remove.ind=(which(coef[j,]==0)+1)
        if(length(remove.ind)<(p-1)){
          fit=glmnet(bag.data[boot.ind.uni,-c(1,remove.ind,(p+2))],
                     bag.data[boot.ind.uni,(p+2)],
                     family=family,weights=boot.weight,
                     lambda=10^(-5))
          if(length(remove.ind)==0){
            coef[j,]=fit$beta[,1]
          }else{
            coef[j,-(remove.ind-1)]=fit$beta[,1]
          }
        }
        
        inclusion.index=which(coef[j,]!=0)
        # turn bool value into int
        var.instances[(j+(i-1)*iter.num),]=0+(colnames(data[,-c(1,p+2)])%in%colnames(data[,-1])[inclusion.index])
      }
      beta_hat=rbind(beta_hat,coef)
      
      std.bag[i,]=apply(coef,2,sd)
      CI.bag[i,]=apply(coef,2,function(x){quantile(x,0.975)-quantile(x,0.025)})
    })[1]
  }
  
  sele.var=apply(var.instances,2,sum)
  sele.ratio=sele.var/(s*iter.num)  
  selected=which(sele.ratio>tau)   
  
  time_cum=cumsum(time)   
  
  
  std_cum=apply(std.bag,2,cumsum)
  std_avg=std_cum
  for(i in 1:nrow(std_cum)){
    std_avg[i,]=std_cum[i,]/rep(i,ncol(std_cum))
  }
  
  
  CI_cum=apply(CI.bag,2,cumsum)
  CI_avg=CI_cum
  for(i in 1:nrow(std_cum)){
    CI_avg[i,]=CI_cum[i,]/rep(i,ncol(CI_cum))
  }
  
  res=list("sele.var"=sele.var,"sele.ratio"=sele.ratio,"selected"=selected,
           "std"=std.bag,"std_avg"=std_avg,
           "c"=CI.bag,"c_avg"=CI_avg,
           "time"=time,"tim_cum"=time_cum,"beta_hat"=beta_hat)
  return(res)
}



BOOT_grp<-function(data,index,iter.num,tau,family){
  p=ncol(data)-2
  coef=matrix(0,iter.num,p)
  var.instances=matrix(0,nrow=iter.num,ncol=p)
  time=vector()
  for(i in 1:iter.num){
    print(i)
    time[i]=system.time({
      boot.data=data[sample(1:nrow(data),nrow(data),replace=T),]
      
      if(family=="gaussian"){
        lambda.max=lambdamax(x=as.matrix(boot.data[,-(p+2)]),
                             y=as.matrix(boot.data[,(p+2)]),
                             index=index,penscale=sqrt,
                             model=LinReg())
        lambda=c(lambda.max*0.96^seq(0,100,1))
        fit_train=quiet(grplasso(x=as.matrix(boot.data[,-(p+2)]),
                                 y=as.matrix(boot.data[,(p+2)]),
                                 index=index,lambda=lambda,
                                 model=LinReg(),penscale=sqrt)) 
      }
      if(family=="binomial"){
        lambda.max=lambdamax(x=as.matrix(boot.data[,-(p+2)]),
                             y=as.matrix(boot.data[,(p+2)]),
                             index=index,penscale=sqrt,
                             model=LogReg())
        lambda=c(lambda.max*0.96^seq(0,100,1))
        fit_train=quiet(grplasso(x=as.matrix(boot.data[,-(p+2)]),
                                 y=as.matrix(boot.data[,(p+2)]),
                                 index=index,lambda=lambda,
                                 model=LogReg(),penscale=sqrt))
      }
      
      bic=rep(Inf,length(lambda))
      for(ii in 1:length(lambda)){
        temp.coef=as.vector(coef(fit_train)[,ii])[-1]
        remove.ind=(which(temp.coef==0)+1)
        if(length(remove.ind)<(p-1)){
          fit=glmnet(boot.data[,-c(1,remove.ind,(p+2))],
                     boot.data[,(p+2)],
                     family=family,lambda=10^(-5))
          tLL=fit$nulldev-deviance(fit)
          k=fit$df
          nn=fit$nobs
          bic[ii]=log(nn)*k-tLL
        }
      }
      min.ind=which.min(bic)
      
      coef[i,]=as.vector(coef(fit_train)[-1,min.ind])
      remove.ind=(which(coef[i,]==0)+1)
      
      if(length(remove.ind)<(p-1)){
        fit=glmnet(boot.data[,-c(1,remove.ind,(p+2))],
                   boot.data[,(p+2)],
                   family=family,lambda=10^(-5))
        if(length(remove.ind)==0){
          coef[i,]=fit$beta[,1]
        }else{
          coef[i,-(remove.ind-1)]=fit$beta[,1]
        }
      }
    })[1]
    inclusion.index=which(coef[i,]!=0)
    var.instances[i,]=0+(colnames(data[,-c(1,p+2)])%in%colnames(data[,-1])[inclusion.index])
  }
  
  sele.var=apply(var.instances,2,sum)
  std.bag=matrix(0,iter.num,p)
  for(i in 2:nrow(coef)){
    for(j in 1:ncol(coef)){
      std.bag[i,j]=sd(coef[1:i,j])
    }
  }
  
  CI.bag=matrix(0,iter.num,p)
  for(i in 2:nrow(coef)){
    for(j in 1:ncol(coef)){
      CI.bag[i,j]=quantile(coef[1:i,j],0.975)-quantile(coef[1:i,j],0.025)
    }
  }
  
  sele.ratio=sele.var/iter.num  
  selected=which(sele.ratio>tau)   
  
  time_cum=cumsum(time)   
  
  res=list("sele.var"=sele.var,"sele.ratio"=sele.ratio,"selected"=selected,
           "std_avg"=std.bag,"c_avg"=CI.bag,
           "time"=time,"tim_cum"=time_cum,"beta_hat"=coef)
  
  return(res)
}


##################################################
###################relative error#################
rel.err0<-function(x,basis,not_zero){
  res=rep(0,nrow(x))
  for(i in 1:nrow(x)){
    res[i]=mean((abs(x[i,not_zero]-basis[not_zero]))/basis[not_zero])
  }
  return(res)
}


####################################################
################### plot ###########################
plotfun<-function(x,y,filename,my_ylab){
  png(filename=filename, width = 25,height = 18, units ="cm", bg = "white", res = 600)
  par(mar=c(4,5,2,2),mfrow=c(1,1),cex.axis=1.8,cex.lab=1.8,lty=1,lwd=2)
  xlim.vec=range(x)
  
  plot(x=x[[1]],y=y[[1]],type="l",
       xlab="Time (sec)",ylab=my_ylab,xlim = xlim.vec,ylim=c(0,1),col="red")
  lines(x=x[[2]],y=y[[2]],type="l",
        col="blue")
  lines(x=x[[3]],y=y[[3]],type="l",
        col="green3")
  lines(x=x[[4]],y=y[[4]],type="l",
        col="orange")
  lines(x=x[[5]],y=y[[5]],type="l",
        col="black",lty=2)
  
  dev.off()
}

seleplot<-function(x,filename){
  png(filename=filename, width = 25,height = 18, units ="cm", bg = "white", res = 600)
  par(mar=c(3,3,2,2),mfrow=c(2,3),cex.axis=2,cex.lab=2,
      cex.main=2)
  mycol=c(rep("red",5),rep("green",4),
          rep("blue",6),rep("yellow",5),
          rep("gray",5),rep("black",5),
          rep("cyan",3),rep("magenta",2))
  df.bar<-barplot(x[[1]],main="BLBVS-0.6",
                  ylab="Proportion",col=mycol)
  lines(x=df.bar,y=rep(0.5,p))
  
  df.bar<-barplot(x[[2]],main="BLBVS-0.7",
                  ylab="Proportion",col=mycol)
  lines(x=df.bar,y=rep(0.5,p))
  
  df.bar<-barplot(x[[3]],main="BLBVS-0.8",
                  ylab="Proportion",col=mycol)
  lines(x=df.bar,y=rep(0.5,p))
  
  df.bar<-barplot(x[[4]],main="BLBVS-0.9",
                  ylab="Proportion",col=mycol)
  lines(x=df.bar,y=rep(0.5,p))
  
  df.bar<-barplot(x[[5]],main="BootVS",
                  ylab="Proportion",col=mycol)
  lines(x=df.bar,y=rep(0.5,p))
  dev.off()
}