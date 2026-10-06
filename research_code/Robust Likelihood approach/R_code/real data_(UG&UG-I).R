#==============================================================================
#-----------real data : body fat percentage data (UG & UG-I model)-------------
#==============================================================================

#read data
bfp_data=read.csv("C:/Users/User/Desktop/R/mountain bro/real data(body fat percentage)/data_set_mglmm.csv")
#setting covariates X & response Y
bfpY=data.frame(bfp_data$ARMS,bfp_data$LEGS,bfp_data$TRUNK,bfp_data$ANDROID,bfp_data$GYNOID)/100
bfpX=data.frame(bfp_data$AGE,bfp_data$BMI,(bfp_data$SEX-1),as.numeric(bfp_data$IPAQ==1),as.numeric(bfp_data$IPAQ==2))
colnames(bfpY)=c("ARMS","LEGS","TRUNK","ANDROID","GYNOID")
colnames(bfpX)=c("AGE","BMI","Male","IA","A")

#package-----------------------------------------------------------------------
library(numDeriv)
library(copula)
library(MASS)
library(pracma)
library(mvtnorm)
#------------------------------------------------------------------------------

#cdf & log likelihood function-------------------------------------------------
unit_gamma_cdf=function(phi,mu,y){
    gf=function(t){
        t^(phi-1)*exp(-t)/factorial(phi-1)
    }
    ugbeta=(mu^(1/phi))/(1-mu^(1/phi))
    ans1=integrate(gf,0,ugbeta*(-log(y)))$value
    ans2=1-ans1
    return(ans2)
}
unit_gamma_pdf=function(phi,mu,y){
    ugbeta=(mu^(1/phi))/(1-mu^(1/phi))
    ans=(ugbeta^phi)*(y^(ugbeta-1))*((-log(y))^(phi-1))/factorial(phi-1)
    return(ans)
}
log_unit_gamma_pdf=function(phi,mu,y){
    term1=log(mu)-phi*log(1-mu^(1/phi))-lgamma(phi)
    term2=((mu^(1/phi))/(1-mu^(1/phi))-1)*log(y)
    term3=(phi-1)*log(-log(y))
    ans=term1+term2+term3
    return(ans)
}
log_dbeta=function(y,a,b){
    ans=log(factorial(a+b-1))-log(factorial(a-1))-log(factorial(b-1))+
        (a-1)*log(y)+(b-1)*log(1-y)
    return(ans)
}

#link function
link=function(x){
    1-1/(1+exp(x))
}
dlink=function(x){
    exp(x)/((1+exp(x))^2)
}
#------------------------------------------------------------------------------

#setting-----------------------------------------------------------------------
dm=5
p=5
n=length(bfpY[,1])
df=50
YYY_ug=as.matrix(bfpY)
XXX=as.matrix(bfpX)
#------------------------------------------------------------------------------


#-----------------------------univariate model ug------------------------------
#(not MGLMM, just regression)
#(take the estimates of univariate model as initial point for UG-I model)
#log likelihood function
loglike_ug=function(para,d){
    logprob=0
    itavalue=XXX%*%matrix(para[3:7],p,1)+para[2]
    muug=sapply(itavalue,link)
    for(i in c(1:n)){
        logprob=logprob+log_unit_gamma_pdf(exp(para[1]),muug[i],YYY_ug[i,d])
    }
    ans=logprob
    return(ans)
}

#arms
reg_ug1=optim(c(0,0,0,0,0,0,0),function(para){-loglike_ug(para,1)},
              method="Nelder-Mead",hessian = T)
updatereg=0
if(reg_ug1$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        temppar=reg_ug1$par
        reg_ug1=optim(temppar,function(para){-loglike_ug(para,1)},
                      method="Nelder-Mead",hessian = T)
        if(reg_ug1$counts[1]<501 | updatereg==100)break
    }
}
reg_ug1par=reg_ug1$par
reg_ug1varm=solve(reg_ug1$hessian)
reg_ug1var=diag(reg_ug1varm)
reg_ug1llh=-reg_ug1$value
testug1=(reg_ug1par^2)/(reg_ug1var)
reg_ug1p=1-pchisq(testug1,1)
#likelihood
reg_ug1llh
#estimates of parameters
reg_ug1par


#legs
reg_ug2=optim(c(0,0,0,0,0,0,0),function(para){-loglike_ug(para,2)},
              method="Nelder-Mead",hessian = T)
updatereg=0
if(reg_ug2$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        temppar=reg_ug2$par
        reg_ug2=optim(temppar,function(para){-loglike_ug(para,2)},
                      method="Nelder-Mead",hessian = T)
        if(reg_ug2$counts[1]<501 | updatereg==100)break
    }
}
reg_ug2par=reg_ug2$par
reg_ug2varm=solve(reg_ug2$hessian)
reg_ug2var=diag(reg_ug2varm)
reg_ug2llh=-reg_ug2$value
testug2=(reg_ug2par^2)/(reg_ug2var)
reg_ug2p=1-pchisq(testug2,1)
#likelihood
reg_ug2llh
#estimates of parameters
reg_ug2par


#trunk
reg_ug3=optim(c(0,0,0,0,0,0,0),function(para){-loglike_ug(para,3)},
              method="Nelder-Mead",hessian = T)
updatereg=0
if(reg_ug3$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        temppar=reg_ug3$par
        reg_ug3=optim(temppar,function(para){-loglike_ug(para,3)},
                      method="Nelder-Mead",hessian = T)
        if(reg_ug3$counts[1]<501 | updatereg==100)break
    }
}
reg_ug3par=reg_ug3$par
reg_ug3varm=solve(reg_ug3$hessian)
reg_ug3var=diag(reg_ug3varm)
reg_ug3llh=-reg_ug3$value
testug3=(reg_ug3par^2)/(reg_ug3var)
reg_ug3p=1-pchisq(testug3,1)
#likelihood
reg_ug3llh
#estimates of parameters
reg_ug3par


#android
reg_ug4=optim(c(0,0,0,0,0,0,0),function(para){-loglike_ug(para,4)},
              method="Nelder-Mead",hessian = T)
updatereg=0
if(reg_ug4$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        temppar=reg_ug4$par
        reg_ug4=optim(temppar,function(para){-loglike_ug(para,4)},
                      method="Nelder-Mead",hessian = T)
        if(reg_ug4$counts[1]<501 | updatereg==100)break
    }
}
reg_ug4par=reg_ug4$par
reg_ug4varm=solve(reg_ug4$hessian)
reg_ug4var=diag(reg_ug4varm)
reg_ug4llh=-reg_ug4$value
testug4=(reg_ug4par^2)/(reg_ug4var)
reg_ug4p=1-pchisq(testug4,1)
#likelihood
reg_ug4llh
#estimates of parameters
reg_ug4par


#gynoid
reg_ug5=optim(c(0,0,0,0,0,0,0),function(para){-loglike_ug(para,5)},
              method="Nelder-Mead",hessian = T)
updatereg=0
if(reg_ug5$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        temppar=reg_ug5$par
        reg_ug5=optim(temppar,function(para){-loglike_ug(para,5)},
                      method="Nelder-Mead",hessian = T)
        if(reg_ug5$counts[1]<501 | updatereg==100)break
    }
}
reg_ug5par=reg_ug5$par
reg_ug5varm=solve(reg_ug5$hessian)
reg_ug5var=diag(reg_ug5varm)
reg_ug5llh=-reg_ug5$value
testug5=(reg_ug5par^2)/(reg_ug5var)
reg_ug5p=1-pchisq(testug5,1)
#likelihood
reg_ug5llh
#estimates of parameters
reg_ug5par

#------------------------------------------------------------------------------



#----------------------univariate MGLMM ug (UG-I model)------------------------
#initial point (from the result of univariate model)
phi=c(reg_ug1par[1],reg_ug2par[1],reg_ug3par[1],reg_ug4par[1],reg_ug5par[1])
phi
bt=c(reg_ug1par[-1],reg_ug2par[-1],reg_ug3par[-1],reg_ug4par[-1],reg_ug5par[-1])
bt
upt=0

#laplace approximation for likelihood function(UG-I)
loglike_ug_1d_laplace=function(para,d){
    upt=upt+1
    print(c("iterating"))
    print(para)
    G=exp(para[8])
    invG=solve(G)
    marg_loglike=0
    itavalue=XXX%*%matrix(para[3:7],p,1)+para[2]
    for(i in c(1:n)){
        f1=function(u){
            muug=sapply((itavalue[i,1]+u[1]),link)
            if(muug>(1-10^-14)){
                muug=(1-10^-14)
            }else if(muug<(10^-14)){
                muug=(10^-14)
            }
            if(para[1]>5){
                para[1]=4.99
            }
            uu=rbind(u[1])
            logprob=log_unit_gamma_pdf(exp(para[1]),muug,YYY_ug[i,d])
            logprobjoint=logprob+log(((2*pi)^(-1/2))*((G)^(-1/2)))+
                (-(1/2)*diag(t(uu)%*%invG%*%uu))
            return(logprobjoint)
        }
        optf1=optim(c(0),function(u){-f1(u)},method="Brent",lower=-5,upper=5,hessian = T)
        hessf1=optf1$hessian
        term1=(-1/2)*log(det(hessf1))
        marg_loglike=marg_loglike+(term1+f1(optf1$par))+1/2*log(2*pi)
    }
    ans=marg_loglike
    return(ans)
}

#====================================arms======================================
#optimize once
optinit=c(phi[1],bt[1:6],log(0.1))
opt_ug=optim(optinit,function(para){-loglike_ug_1d_laplace(para,1)},
             method="Nelder-Mead",hessian = T)

#Scattering points to check it is good place to find MLE
set.seed(52210)
count=0
if(opt_ug$counts[1]>=501){
    temppar=opt_ug$par
    value1d=opt_ug$value
    for(k in c(1:1000)){
        temptemppar=c(temppar+rnorm(8,0,0.0005))
        tempvalue=-loglike_ug_1d_laplace(temptemppar,1)
        if(tempvalue<value1d){
            temppar=temptemppar
            value1d=tempvalue
            count=count+1
        }
    }
    print(count)
}
HM1=hessian(function(para){-loglike_ug_1d_laplace(para,1)},temppar)
egvalue_ug=eigen(HM1)$values
detect=length(egvalue_ug[egvalue_ug<=0])
if(count>0){
    detect=detect+1
}

#check the hessian matrix is positive-define
#if the matrix is not positive-define, repeat optim()
if(detect>0){
    update_times=0
    repeat{
        inpar=temppar*0.75+optinit*0.25
        g=rnorm(8,0,0.0005)
        init=inpar+g
        opt_ug_best=optim(init,function(para){-loglike_ug_1d_laplace(para,1)},
                          method="Nelder-Mead",hessian = T)
        egvalue_ug=eigen(opt_ug_best$hessian)$values
        detect=length(egvalue_ug[egvalue_ug<=0])
        update_times=update_times+1
        temppar=opt_ug_best$par
        if(detect==0 | update_times==15)break
    }
    HM1=opt_ug_best$hessian
}
#(double check)
#To ensure it is maximum, repeat optim() until convergence
#if not convergence, run #-----*** again (optim$count<501)
#-----***
update_times=0
repeat{
    inpar=temppar
    init=inpar
    opt_ug_best=optim(init,function(para){-loglike_ug_1d_laplace(para,1)},
                      method="Nelder-Mead",hessian = T)
    egvalue_ug=eigen(opt_ug_best$hessian)$values
    detect=length(egvalue_ug[egvalue_ug<=0])
    update_times=update_times+1
    temppar=opt_ug_best$par
    if(detect==0 | update_times==15)break
}
#-----***

HM1=opt_ug_best$hessian
opt_ug_best$counts
opt_ug_uv1par=temppar
opt_ug_uv1varm=solve(HM1)
opt_ug_uv1var=diag(opt_ug_uv1varm)
testug_uv1=(opt_ug_uv1par^2)/opt_ug_uv1var
testug_uv1p=1-pchisq(testug_uv1,1)
CIug_uv1L=c(opt_ug_uv1par-qnorm(0.975)*sqrt(opt_ug_uv1var))
CIug_uv1U=c(opt_ug_uv1par+qnorm(0.975)*sqrt(opt_ug_uv1var))
loglikug_uv1=loglike_ug_1d_laplace(opt_ug_uv1par,1)
AICug_uv1=(-2)*loglikug_uv1+2*8
BICug_uv1=(-2)*loglikug_uv1+log(n)*8

#result
#estimates of parameters
opt_ug_uv1par
#estimates of variance & standard deviation
opt_ug_uv1var
sqrt(opt_ug_uv1var)
#Wald test p-value
testug_uv1p
#Confident interval of Wald test
CIug_uv1L
CIug_uv1U
#likelihood, AIC & BIC
loglikug_uv1
AICug_uv1
BICug_uv1
#CI length
2*qnorm(0.975)*sqrt(opt_ug_uv1var)
#==============================================================================


#====================================legs======================================
#optimize once
optinit=c(phi[2],bt[7:12],log(0.1))
opt_ug=optim(optinit,function(para){-loglike_ug_1d_laplace(para,2)},
             method="Nelder-Mead",hessian = T)

#Scattering points to check it is good place to find MLE
set.seed(52210)
count=0
if(opt_ug$counts[1]>=501){
    temppar=opt_ug$par
    value1d=opt_ug$value
    for(k in c(1:1000)){
        temptemppar=c(temppar+rnorm(8,0,0.0005))
        tempvalue=-loglike_ug_1d_laplace(temptemppar,2)
        if(tempvalue<value1d){
            temppar=temptemppar
            value1d=tempvalue
            count=count+1
        }
    }
    print(count)
}
HM2=hessian(function(para){-loglike_ug_1d_laplace(para,2)},temppar)
egvalue_ug=eigen(HM2)$values
detect=length(egvalue_ug[egvalue_ug<=0])
if(count>0){
    detect=detect+1
}

#check the hessian matrix is positive-define
#if the matrix is not positive-define, repeat optim()
if(detect>0){
    update_times=0
    repeat{
        inpar=temppar*0.75+optinit*0.25
        g=rnorm(8,0,0.0005)
        init=inpar+g
        opt_ug_best=optim(init,function(para){-loglike_ug_1d_laplace(para,2)},
                          method="Nelder-Mead",hessian = T)
        egvalue_ug=eigen(opt_ug_best$hessian)$values
        detect=length(egvalue_ug[egvalue_ug<=0])
        update_times=update_times+1
        temppar=opt_ug_best$par
        if(detect==0 | update_times==15)break
    }
    HM2=opt_ug_best$hessian
}
#(double check)
#To ensure it is maximum, repeat optim() until convergence
#if not convergence, run #-----*** again (optim$count<501)
#-----***
update_times=0
repeat{
    inpar=temppar
    init=inpar
    opt_ug_best=optim(init,function(para){-loglike_ug_1d_laplace(para,2)},
                      method="Nelder-Mead",hessian = T)
    egvalue_ug=eigen(opt_ug_best$hessian)$values
    detect=length(egvalue_ug[egvalue_ug<=0])
    update_times=update_times+1
    temppar=opt_ug_best$par
    if(detect==0 | update_times==15)break
}
#-----***

HM2=opt_ug_best$hessian
opt_ug_best$counts
opt_ug_uv2par=temppar
opt_ug_uv2varm=solve(HM2)
opt_ug_uv2var=diag(opt_ug_uv2varm)
testug_uv2=(opt_ug_uv2par^2)/opt_ug_uv2var
testug_uv2p=1-pchisq(testug_uv2,1)
CIug_uv2L=c(opt_ug_uv2par-qnorm(0.975)*sqrt(opt_ug_uv2var))
CIug_uv2U=c(opt_ug_uv2par+qnorm(0.975)*sqrt(opt_ug_uv2var))
loglikug_uv2=loglike_ug_1d_laplace(opt_ug_uv2par,2)
AICug_uv2=(-2)*loglikug_uv2+2*8
BICug_uv2=(-2)*loglikug_uv2+log(n)*8

#result
#estimates of parameters
opt_ug_uv2par
#estimates of variance & standard deviation
opt_ug_uv2var
sqrt(opt_ug_uv2var)
#Wald test p-value
testug_uv2p
#Confident interval of Wald test
CIug_uv2L
CIug_uv2U
#likelihood, AIC & BIC
loglikug_uv2
AICug_uv2
BICug_uv2
#CI length
2*qnorm(0.975)*sqrt(opt_ug_uv2var)
#==============================================================================


#====================================trunk=====================================
#optimize once
optinit=c(phi[3],bt[13:18],log(0.1))
opt_ug=optim(optinit,function(para){-loglike_ug_1d_laplace(para,3)},
             method="Nelder-Mead",hessian = T)

#Scattering points to check it is good place to find MLE
set.seed(52210)
count=0
if(opt_ug$counts[1]>=501){
    temppar=opt_ug$par
    value1d=opt_ug$value
    for(k in c(1:1000)){
        temptemppar=c(temppar+rnorm(8,0,0.0005))
        tempvalue=-loglike_ug_1d_laplace(temptemppar,3)
        if(tempvalue<value1d){
            temppar=temptemppar
            value1d=tempvalue
            count=count+1
        }
    }
    print(count)
}
HM3=hessian(function(para){-loglike_ug_1d_laplace(para,3)},temppar)
egvalue_ug=eigen(HM3)$values
detect=length(egvalue_ug[egvalue_ug<=0])
if(count>0){
    detect=detect+1
}

#check the hessian matrix is positive-define
#if the matrix is not positive-define, repeat optim()
if(detect>0){
    update_times=0
    repeat{
        inpar=temppar*0.75+optinit*0.25
        g=rnorm(8,0,0.0005)
        init=inpar+g
        opt_ug_best=optim(init,function(para){-loglike_ug_1d_laplace(para,3)},
                          method="Nelder-Mead",hessian = T)
        egvalue_ug=eigen(opt_ug_best$hessian)$values
        detect=length(egvalue_ug[egvalue_ug<=0])
        update_times=update_times+1
        temppar=opt_ug_best$par
        if(detect==0 | update_times==15)break
    }
    HM3=opt_ug_best$hessian
}
#(double check)
#To ensure it is maximum, repeat optim() until convergence
#if not convergence, run #-----*** again (optim$count<501)
#-----***
if(opt_ug_best$value+reg_ug3llh>0){
    update_times=0
    repeat{
        inpar=temppar
        init=inpar
        opt_ug_best=optim(init,function(para){-loglike_ug_1d_laplace(para,3)},
                          method="Nelder-Mead",hessian = T)
        egvalue_ug=eigen(opt_ug_best$hessian)$values
        detect=length(egvalue_ug[egvalue_ug<=0])
        update_times=update_times+1
        temppar=opt_ug_best$par
        if(detect==0 | update_times==15)break
    }
    HM3=opt_ug_best$hessian
}
#-----***

opt_ug_best$value+reg_ug3llh
opt_ug_best$counts
opt_ug_uv3par=temppar
opt_ug_uv3varm=solve(HM3)
opt_ug_uv3var=diag(opt_ug_uv3varm)
testug_uv3=(opt_ug_uv3par^2)/opt_ug_uv3var
testug_uv3p=1-pchisq(testug_uv3,1)
CIug_uv3L=c(opt_ug_uv3par-qnorm(0.975)*sqrt(opt_ug_uv3var))
CIug_uv3U=c(opt_ug_uv3par+qnorm(0.975)*sqrt(opt_ug_uv3var))
loglikug_uv3=loglike_ug_1d_laplace(opt_ug_uv3par,3)
AICug_uv3=(-2)*loglikug_uv3+2*8
BICug_uv3=(-2)*loglikug_uv3+log(n)*8

#result
#estimates of parameters
opt_ug_uv3par
#estimates of variance & standard deviation
opt_ug_uv3var
sqrt(opt_ug_uv3var)
#Wald test p-value
testug_uv3p
#Confident interval of Wald test
CIug_uv3L
CIug_uv3U
#likelihood, AIC & BIC
loglikug_uv3
AICug_uv3
BICug_uv3
#CI length
2*qnorm(0.975)*sqrt(opt_ug_uv3var)
#==============================================================================


#===================================android====================================
#optimize once
optinit=c(phi[4],bt[19:24],log(0.1))
opt_ug=optim(optinit,function(para){-loglike_ug_1d_laplace(para,4)},
             method="Nelder-Mead",hessian = T)

#Scattering points to check it is good place to find MLE
set.seed(52210)
count=0
if(opt_ug$counts[1]>=501){
    temppar=opt_ug$par
    value1d=opt_ug$value
    for(k in c(1:1000)){
        temptemppar=c(temppar+rnorm(8,0,0.0005))
        tempvalue=-loglike_ug_1d_laplace(temptemppar,4)
        if(tempvalue<value1d){
            temppar=temptemppar
            value1d=tempvalue
            count=count+1
        }
    }
    print(count)
}
HM4=hessian(function(para){-loglike_ug_1d_laplace(para,4)},temppar)
egvalue_ug=eigen(HM4)$values
detect=length(egvalue_ug[egvalue_ug<=0])
if(count>0){
    detect=detect+1
}

#check the hessian matrix is positive-define
#if the matrix is not positive-define, repeat optim()
if(detect>0){
    update_times=0
    repeat{
        inpar=temppar*0.75+optinit*0.25
        g=rnorm(8,0,0.0005)
        init=inpar+g
        opt_ug_best=optim(init,function(para){-loglike_ug_1d_laplace(para,4)},
                          method="Nelder-Mead",hessian = T)
        egvalue_ug=eigen(opt_ug_best$hessian)$values
        detect=length(egvalue_ug[egvalue_ug<=0])
        update_times=update_times+1
        temppar=opt_ug_best$par
        if(detect==0 | update_times==15)break
    }
    HM4=opt_ug_best$hessian
}
#(double check)
#To ensure it is maximum, repeat optim() until convergence
#if not convergence, run #-----*** again (optim$count<501)
#-----***
update_times=0
repeat{
    inpar=temppar
    init=inpar
    opt_ug_best=optim(init,function(para){-loglike_ug_1d_laplace(para,4)},
                      method="Nelder-Mead",hessian = T)
    egvalue_ug=eigen(opt_ug_best$hessian)$values
    detect=length(egvalue_ug[egvalue_ug<=0])
    update_times=update_times+1
    temppar=opt_ug_best$par
    if(detect==0 | update_times==15)break
}
#-----***

HM4=opt_ug_best$hessian
opt_ug_best$counts
opt_ug_uv4par=temppar
opt_ug_uv4varm=solve(HM4)
opt_ug_uv4var=diag(opt_ug_uv4varm)
testug_uv4=(opt_ug_uv4par^2)/opt_ug_uv4var
testug_uv4p=1-pchisq(testug_uv4,1)
CIug_uv4L=c(opt_ug_uv4par-qnorm(0.975)*sqrt(opt_ug_uv4var))
CIug_uv4U=c(opt_ug_uv4par+qnorm(0.975)*sqrt(opt_ug_uv4var))
loglikug_uv4=loglike_ug_1d_laplace(opt_ug_uv4par,4)
AICug_uv4=(-2)*loglikug_uv4+2*8
BICug_uv4=(-2)*loglikug_uv4+log(n)*8

#result
#estimates of parameters
opt_ug_uv4par
#estimates of variance & standard deviation
opt_ug_uv4var
sqrt(opt_ug_uv4var)
#Wald test p-value
testug_uv4p
#Confident interval of Wald test
CIug_uv4L
CIug_uv4U
#likelihood, AIC & BIC
loglikug_uv4
AICug_uv4
BICug_uv4
#CI length
2*qnorm(0.975)*sqrt(opt_ug_uv4var)
#==============================================================================


#===================================gynoid=====================================
#optimize once
optinit=c(phi[5],bt[25:30],log(0.1))
opt_ug=optim(optinit,function(para){-loglike_ug_1d_laplace(para,5)},
             method="Nelder-Mead",hessian = T)

#Scattering points to check it is good place to find MLE
set.seed(52210)
count=0
if(opt_ug$counts[1]>=501){
    temppar=opt_ug$par
    value1d=opt_ug$value
    for(k in c(1:1000)){
        temptemppar=c(temppar+rnorm(8,0,0.0005))
        tempvalue=-loglike_ug_1d_laplace(temptemppar,5)
        if(tempvalue<value1d){
            temppar=temptemppar
            value1d=tempvalue
            count=count+1
        }
    }
    print(count)
}
HM5=hessian(function(para){-loglike_ug_1d_laplace(para,5)},temppar)
egvalue_ug=eigen(HM5)$values
detect=length(egvalue_ug[egvalue_ug<=0])
if(count>0){
    detect=detect+1
}

#check the hessian matrix is positive-define
#if the matrix is not positive-define, repeat optim()
if(detect>0){
    update_times=0
    repeat{
        inpar=temppar*0.75+optinit*0.25
        g=rnorm(8,0,0.0005)
        init=inpar+g
        opt_ug_best=optim(init,function(para){-loglike_ug_1d_laplace(para,5)},
                          method="Nelder-Mead",hessian = T)
        egvalue_ug=eigen(opt_ug_best$hessian)$values
        detect=length(egvalue_ug[egvalue_ug<=0])
        update_times=update_times+1
        temppar=opt_ug_best$par
        if(detect==0 | update_times==15)break
    }
    HM5=opt_ug_best$hessian
}
#(double check)
#To ensure it is maximum, repeat optim() until convergence
#if not convergence, run #-----*** again (optim$count<501)
#-----***
update_times=0
repeat{
    inpar=temppar
    init=inpar
    opt_ug_best=optim(init,function(para){-loglike_ug_1d_laplace(para,5)},
                      method="Nelder-Mead",hessian = T)
    egvalue_ug=eigen(opt_ug_best$hessian)$values
    detect=length(egvalue_ug[egvalue_ug<=0])
    update_times=update_times+1
    temppar=opt_ug_best$par
    if(detect==0 | update_times==15)break
}
#-----***

HM5=opt_ug_best$hessian
opt_ug_best$counts
opt_ug_uv5par=temppar
opt_ug_uv5varm=solve(HM5)
opt_ug_uv5var=diag(opt_ug_uv5varm)
testug_uv5=(opt_ug_uv5par^2)/opt_ug_uv5var
testug_uv5p=1-pchisq(testug_uv5,1)
CIug_uv5L=c(opt_ug_uv5par-qnorm(0.975)*sqrt(opt_ug_uv5var))
CIug_uv5U=c(opt_ug_uv5par+qnorm(0.975)*sqrt(opt_ug_uv5var))
loglikug_uv5=loglike_ug_1d_laplace(opt_ug_uv5par,5)
AICug_uv5=(-2)*loglikug_uv5+2*8
BICug_uv5=(-2)*loglikug_uv5+log(n)*8

#result
#estimates of parameters
opt_ug_uv5par
#estimates of variance & standard deviation
opt_ug_uv5var
sqrt(opt_ug_uv5var)
#Wald test p-value
testug_uv5p
#Confident interval of Wald test
CIug_uv5L
CIug_uv5U
#likelihood, AIC & BIC
loglikug_uv5
AICug_uv5
BICug_uv5
#CI length
2*qnorm(0.975)*sqrt(opt_ug_uv5var)
#==============================================================================
#------------------------------------------------------------------------------


#---------------------------------MGLMM ug(2d)---------------------------------
#(arms & trunk) => d=1,t=3 ; (legs & android) => d=2,t=4
#initial point (from the result of UG-I model)
opt_ug_2dtable=matrix(0,9,10)
phi=c(opt_ug_uv1par[1],opt_ug_uv2par[1],opt_ug_uv3par[1],opt_ug_uv4par[1],opt_ug_uv5par[1])
phi
bt=c(opt_ug_uv1par[c(2:7)],opt_ug_uv2par[c(2:7)],opt_ug_uv3par[c(2:7)],opt_ug_uv4par[c(2:7)],opt_ug_uv5par[c(2:7)])
bt
sg=(c(opt_ug_uv1par[8],opt_ug_uv2par[8],opt_ug_uv3par[8],opt_ug_uv4par[8],opt_ug_uv5par[8]))/2
sg
d=2
t=4
optinit=c(phi[c(d,t)],bt[c(((p+1)*d-5):((p+1)*d),((p+1)*t-5):((p+1)*t))],sg[c(d,t)],0)

#laplace approximation for likelihood function
loglike_ug_2d_laplace=function(para,d,t){
    print(c("iterating"))
    print(para)
    hG=diag(1,2)
    hG[lower.tri(hG)]=para[17]
    hG=hG%*%t(hG)
    dghG=diag((diag(hG))^(-1/2),2)
    hG=dghG%*%hG%*%dghG
    G=diag(exp(para[15:16]),2)%*%hG%*%diag(exp(para[15:16]),2)
    invG=solve(G)
    marg_loglike=rep(0,n)
    itavalue1=XXX%*%matrix(para[4:8],p,1)+para[3]
    itavalue2=XXX%*%matrix(para[10:14],p,1)+para[9]
    for(i in c(1:n)){
        f1=function(u){
            muug1=sapply((itavalue1[i,1]+u[1]),link)
            muug2=sapply((itavalue2[i,1]+u[2]),link)
            if(muug1>(1-10^-10)){
                muug1=(1-10^-10)
            }else if(muug1<(10^-10)){
                muug1=(10^-10)
            }
            if(muug2>(1-10^-10)){
                muug2=(1-10^-10)
            }else if(muug2<(10^-10)){
                muug2=(10^-10)
            }
            uu=rbind(u[1],u[2])
            a=para[1]
            b=para[2]
            logprob=log_unit_gamma_pdf(exp(a),muug1,YYY_ug[i,d])+
                log_unit_gamma_pdf(exp(b),muug2,YYY_ug[i,t])
            logprobjoint=logprob+log(((2*pi)^(-2/2))*(det(G)^(-1/2)))+
                (-(1/2)*diag(t(uu)%*%invG%*%uu))
            return(logprobjoint)
        }
        optf1=optim(c(0,0),function(u){-f1(u)},hessian = T)
        hessf1=optf1$hessian
        term1=(-1/2)*log(det(hessf1))
        marg_loglike[i]=(term1+f1(optf1$par))+2/2*log(2*pi)
    }
    ans=sum(marg_loglike)
    return(ans)
}
optinit

#optimize
opt_ug=optim(c(optinit[1:17]),function(para){-loglike_ug_2d_laplace(c(para),d,t)},
             method="Nelder-Mead",hessian = T)

#To ensure it is maximum, repeat optim() until convergence
opt_ug_best=opt_ug
egvalue_ug=eigen(opt_ug$hessian)$values
detect=length(egvalue_ug[egvalue_ug<0])

#d=2,t=4 run this code=========================================================
ug_update=0
if(detect>0){
    ug_update=ug_update+1
    update_times=0
    repeat{
        inpar=c(opt_ug_best$par)
        init=inpar
        opt_ug_best=optim(init,function(para){-loglike_ug_2d_laplace(c(para),d,t)},
                          method="Nelder-Mead",hessian = T)
        egvalue_ug=eigen(opt_ug_best$hessian)$values
        detect=length(egvalue_ug[egvalue_ug<0])
        update_times=update_times+1
        print(c("ug_regression_updating",update_times,detect,opt_ug_best$value))
        if(detect==0 | update_times==25)break
    }
}
#==============================================================================

#d=1,t=3 run this code=========================================================
if(detect>0){
    ug_update=ug_update+1
    update_times=0
    inpar=c(opt_ug_best$par[1:16],log(opt_ug_best$par[17]))
    repeat{
        init=inpar
        opt_ug_best=optim(init,function(para){-loglike_ug_2d_laplace(c(para[1:16],exp(para[17])),d,t)},
                          method="Nelder-Mead",hessian = T)
        inpar=c(opt_ug_best$par)
        egvalue_ug=eigen(opt_ug_best$hessian)$values
        detect=length(egvalue_ug[egvalue_ug<0])
        update_times=update_times+1
        print(c("ug_regression_updating",update_times,detect,opt_ug_best$value))
        if(detect==0 | update_times==25)break
    }
}
#==============================================================================

opt_ug_uvdtpar=opt_ug_best$par
opt_ug_uvdtvarm=solve(opt_ug_best$hessian)
opt_ug_uvdtvar=diag(opt_ug_uvdtvarm)
loglikug_uvdt=-opt_ug_best$value
testug_uvdt=(opt_ug_uvdtpar^2)/opt_ug_uvdtvar
testug_uvdtp=1-pchisq(testug_uvdt,1)
CIug_uvdtL=c(opt_ug_uvdtpar-qnorm(0.975)*sqrt(opt_ug_uvdtvar))
CIug_uvdtU=c(opt_ug_uvdtpar+qnorm(0.975)*sqrt(opt_ug_uvdtvar))

#result
#estimates of parameters
opt_ug_uvdtpar
#estimates of variance & standard deviation
opt_ug_uvdtvar
sqrt(opt_ug_uvdtvar)
#likelihood
loglikug_uvdt
#Wald test p-value
testug_uvdtp
#Confident interval of Wald test
CIug_uvdtL
CIug_uvdtU

#variance matrix of random effect(d=2,t=4)
hG=diag(1,2)
hG[lower.tri(hG)]=opt_ug_uvdtpar[17]
hG=hG%*%t(hG)
dghG=diag((diag(hG))^(-1/2),2)
hG=dghG%*%hG%*%dghG
G=diag(exp(opt_ug_uvdtpar[15:16]),2)%*%hG%*%diag(exp(opt_ug_uvdtpar[15:16]),2)
hG
G

#variance matrix of random effect(d=1,t=3)
hG=diag(1,2)
hG[lower.tri(hG)]=exp(opt_ug_uvdtpar[17])
hG=hG%*%t(hG)
dghG=diag((diag(hG))^(-1/2),2)
hG=dghG%*%hG%*%dghG
G=diag(exp(opt_ug_uvdtpar[15:16]),2)%*%hG%*%diag(exp(opt_ug_uvdtpar[15:16]),2)
hG
G
#------------------------------------------------------------------------------

