#==============================================================================
#-------------simulation part 4_UG (independent unit gamma model)--------------
#==============================================================================

#package-----------------------------------------------------------------------
library(numDeriv)
library(copula)
library(MASS)
library(pracma)
library(mvtnorm)
library(tictoc)
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
    ans=lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(y)+(b-1)*log(1-y)
    return(ans)
}
link=function(x){
    1-1/(1+exp(x))
}
dlink=function(x){
    exp(x)/((1+exp(x))^2)
}
#------------------------------------------------------------------------------

#start-------------------------------------------------------------------------
#simulation times N=500, sample size n=100
#beta_t=(-1,1,-0.5,1) & XXX1~N(0,0.3^2), XXX2~Ber(0.5)
#unit gamma model precision parameters (10,5)
#Normal copula omega=0.6
#------------------------------------------------------------------------------
#setting
N=500
n=100
dm=2
beta_t=c(-1,1,-0.5,1)
p=length(beta_t)/dm-1
YYY_mix=matrix(0,n,dm)
aphb=10
aphug=5

#mu setting
set.seed(210)
XXX1=matrix(rnorm(n*dm,0,0.3),n*dm,1)
XXX2=matrix(rbinom(n*dm,1,0.5),n*dm,1)
XXX=XXX1
it1=XXX[1:n,]%*%matrix(beta_t[c(2)],p,1)+beta_t[1]
it2=XXX[(n+1):(n*dm),]%*%matrix(beta_t[c(4)],p,1)+beta_t[3]
it=rep(0,dm*n)
it[c(c(1:n)*dm-(dm-1))]=it1
it[c(c(1:n)*dm-(dm-2))]=it2
mufix=sapply(it,link)
matrix(mufix,n,dm,byrow=T)
mu=mufix
XXXtp=matrix(0,n*dm,p)
XXXtp[c(c(1:n)*dm-(dm-1)),]=XXX[1:n,]
XXXtp[c(c(1:n)*dm-(dm-2)),]=XXX[(n+1):(n*dm),]
XXXtp

#Normal copula setting
normal_cupcor_mix=0.6

#record the result
phi_sim_mix=matrix(0,dm,N)
phi_avar_mix=matrix(0,dm,N)
beta_sim_mix=matrix(0,(p+1)*dm,N)
beta_dllhvalue_mix=matrix(0,N)
beta_avar_mix=matrix(0,(p+1)*dm,N)
walds_beta_mix=matrix(0,(p+1)*dm,N)
sg_sim_mix=matrix(0,3,N)
sg_avar_mix=matrix(0,3,N)

#unit gamma regression---------------------------------------------------------
mix_update=0
a=1
upt=rep(0,500)
tic.clearlog()
detect_sim1=rep(0,500)
detect_sim2=rep(0,500)
for(M in c(1:N)){
    
    #time start
    tic(M)
    
    #generating data
    update_times=0
    set.seed(829+M*89)
    for(i in c(1:n)){
        mixcopula=normalCopula(normal_cupcor_mix,dim=dm,dispstr = "un")
        mvmix=mvdc(mixcopula,margins=c("gamma","gamma"),
                   paramMargins = list(list(shape=aphb,scale=((mu[dm*i-dm+1])^(-1/aphb)-1)),
                                       list(shape=aphug,scale=((mu[dm*i-dm+2])^(-1/aphug)-1))))
        temp=c(rMvdc(1,mvmix))
        YYY_mix[i,]=c(exp(-temp[1]),exp(-temp[2]))
    }
    
    ###dimension1--------------------------------------------------------------
    #laplace approximation for likelihood function
    loglike_mix_1d_laplace1=function(para){
        print(para)
        print(c(1,update_times))
        G=matrix(exp(para[4])^2,1,1)
        invG=solve(G)
        marg_loglike=rep(0,n)
        itavalue1=matrix(XXX[1:n,]%*%matrix(para[3],p,1)+para[2],n,1)
        itavalue=itavalue1
        for(i in c(1:n)){
            f1=function(u){
                mub1=sapply((itavalue[i,1]+u[1]),link)
                if(mub1<(10^(-14))){
                    mub1=(10^(-14))
                }else if(mub1>(1-10^(-14))){
                    mub1=(1-10^(-14))
                }
                uu=rbind(u[1])
                logprob=log_unit_gamma_pdf(exp(para[1]),mub1,YYY_mix[i,1])
                logprobjoint=logprob+log(((2*pi)^(-1/2))*(det(G)^(-1/2)))+
                    (-(1/2)*diag(t(uu)%*%invG%*%uu))
                return(logprobjoint)
            }
            optf1=optim(c(0),function(u){-f1(u)},hessian = T)
            hessf1=optf1$hessian
            term1=(-1/2)*log(det(hessf1))
            marg_loglike[i]=(term1+f1(optf1$par))+(1/2)*log(2*pi)
        }
        ans=sum(marg_loglike)
        return(ans)
    }
    
    #initial point
    optinit1=c(log(aphb),beta_t[1:2],-5)
    
    #optimize
    opt_mix1=optim(optinit1,function(para){-loglike_mix_1d_laplace1(c(para))},
                   method="Nelder-Mead",hessian = T)
    
    #To ensure it is maximum, repeat optim() until convergence
    opt_mix_best1=opt_mix1
    egvalue_mix1=eigen(opt_mix1$hessian)$values
    detect1=length(egvalue_mix1[egvalue_mix1<0])
    if(detect1>0){
        mix_update=mix_update+1
        inpar1=opt_mix1$par
        loglike1=opt_mix1$value
        repeat{
            a=1
            init1=inpar1+rnorm(4,0,0.0001)
            opt_mix_best1=optim(init1,function(para){-loglike_mix_1d_laplace1(c(para))},
                                method="Nelder-Mead",hessian = T)
            egvalue_mix1=eigen(opt_mix_best1$hessian)$values
            detect1=length(egvalue_mix1[egvalue_mix1<0]) 
            update_times=update_times+1
            print(c(M,"mix_regression1_updating",update_times,detect1))
            #first-check the hessian matrix is positive-define
            #if the matrix is positive-define, stop optim()
            if(detect1==0)break
            #second-reach convergence tolerance
            #if absolute error<=0.0000000001,also stop optim()
            if(opt_mix_best1$value<loglike1){
                a=abs(opt_mix_best1$value-loglike1)
                inpar1=opt_mix_best1$par
                loglike1=opt_mix_best1$value
            }
            if(a<=0.0000000001)break
        }
    }
    mix_varhat_best1=solve(opt_mix_best1$hessian) 
    beta_sim_mix[1:2,M]=opt_mix_best1$par[2:3]
    beta_avar_mix[1:2,M]=diag(mix_varhat_best1)[2:3]
    
    #other parameters
    phi_sim_mix[1,M]=opt_mix_best1$par[1]
    phi_avar_mix[1,M]=diag(mix_varhat_best1)[1]
    sg_sim_mix[1,M]=opt_mix_best1$par[4]
    sg_avar_mix[1,M]=diag(mix_varhat_best1)[4]
    ###------------------------------------------------------------------------
    
    ###dimension2--------------------------------------------------------------
    #laplace approximation for likelihood function
    loglike_mix_1d_laplace2=function(para){
        print(para)
        print(c(2,update_times))
        G=matrix(exp(para[4])^2,1,1)
        invG=solve(G)
        marg_loglike=rep(0,n)
        itavalue2=matrix(XXX[(n+1):(n*dm),]%*%matrix(para[3],p,1)+para[2],n,1)
        itavalue=itavalue2
        for(i in c(1:n)){
            f1=function(u){
                muug2=sapply((itavalue[i,1]+u[1]),link)
                if(muug2<(10^(-14))){
                    muug2=(10^(-14))
                }else if(muug2>(1-10^(-14))){
                    muug2=(1-10^(-14))
                }
                uu=rbind(u[1])
                logprob=log_unit_gamma_pdf(exp(para[1]),muug2,YYY_mix[i,2])
                logprobjoint=logprob+log(((2*pi)^(-1/2))*(det(G)^(-1/2)))+
                    (-(1/2)*diag(t(uu)%*%invG%*%uu))
                return(logprobjoint)
            }
            optf1=optim(c(0),function(u){-f1(u)},hessian = T)
            hessf1=optf1$hessian
            term1=(-1/2)*log(det(hessf1))
            marg_loglike[i]=(term1+f1(optf1$par))+(1/2)*log(2*pi)
        }
        ans=sum(marg_loglike)
        return(ans)
    }
    
    #initial point
    optinit2=c(log(aphug),beta_t[3:4],-5)
    
    #optimize
    opt_mix2=optim(optinit2,function(para){-loglike_mix_1d_laplace2(c(para))},
                   method="Nelder-Mead",hessian = T)
    
    #To ensure it is maximum, repeat optim() until convergence
    opt_mix_best2=opt_mix2
    egvalue_mix2=eigen(opt_mix2$hessian)$values
    detect2=length(egvalue_mix2[egvalue_mix2<0])
    if(detect2>0){
        mix_update=mix_update+1
        inpar2=opt_mix2$par
        loglike2=opt_mix2$value
        repeat{
            a=1
            init2=inpar2+rnorm(4,0,0.0001)
            opt_mix_best2=optim(init2,function(para){-loglike_mix_1d_laplace2(c(para))},
                                method="Nelder-Mead",hessian = T)
            egvalue_mix2=eigen(opt_mix_best2$hessian)$values
            detect2=length(egvalue_mix2[egvalue_mix2<0]) 
            update_times=update_times+1
            print(c(M,"mix_regression2_updating",update_times,detect2))
            #first-check the hessian matrix is positive-define
            #if the matrix is positive-define, stop optim()
            if(detect2==0)break
            #second-reach convergence tolerance
            #if absolute error<=0.0000000001,also stop optim()
            if(opt_mix_best2$value<loglike2){
                a=abs(opt_mix_best2$value-loglike2)
                inpar2=opt_mix_best2$par
                loglike2=opt_mix_best2$value
            }
            if(a<=0.0000000001)break
        }
    }
    mix_varhat_best2=solve(opt_mix_best2$hessian) 
    beta_sim_mix[3:4,M]=opt_mix_best2$par[2:3]
    beta_avar_mix[3:4,M]=diag(mix_varhat_best2)[2:3]
    
    #other parameters
    phi_sim_mix[2,M]=opt_mix_best2$par[1]
    phi_avar_mix[2,M]=diag(mix_varhat_best2)[1]
    sg_sim_mix[2,M]=opt_mix_best2$par[4]
    sg_avar_mix[2,M]=diag(mix_varhat_best2)[4]
    ###------------------------------------------------------------------------
    
    #wald test
    walds_beta_mix[1:4,M]=((beta_sim_mix[1:4,M]-beta_t)^2)/beta_avar_mix[1:4,M]
    print(c(M,"mix_regression complete"))
    
    #number of eigenvalues<=0 in each simulation
    detect_sim1[M]=detect1
    detect_sim2[M]=detect2
    upt[M]=update_times
    
    #time stop
    toc(log = TRUE, quiet = TRUE)
}
log.lst <- tic.log(format = FALSE)
timings <- unlist(lapply(log.lst, function(x){x$toc - x$tic}))
#average of taking time
mean(timings)
#convergence times
detectok1=length(detect_sim1[detect_sim1==0])
detectok1/500
detectok2=length(detect_sim2[detect_sim2==0])
detectok2/500
upt
#------------------------------------------------------------------------------

#unit gamma result-------------------------------------------------------------
beta_sim_mix
beta_dllhvalue_mix
beta_avar_mix
mix_update
#dimension1
#mean, sample variance, variance hat result
mean(beta_sim_mix[1,1:500])
mean(beta_sim_mix[2,1:500])
var(beta_sim_mix[1,1:500])
var(beta_sim_mix[2,1:500])
mean(beta_avar_mix[1,1:500])
mean(beta_avar_mix[2,1:500])
#wald test result
reject_W_mix=+(walds_beta_mix>qchisq(0.95,1))
mean(reject_W_mix[1,1:500])
mean(reject_W_mix[2,1:500])

#dimension2
#mean, sample variance, variance hat result
mean(beta_sim_mix[3,1:500])
mean(beta_sim_mix[4,1:500])
var(beta_sim_mix[3,1:500])
var(beta_sim_mix[4,1:500])
mean(beta_avar_mix[3,1:500])
mean(beta_avar_mix[4,1:500])
#wald test result
mean(reject_W_mix[3,1:500])
mean(reject_W_mix[4,1:500])

#other parameters
#phi
mean(exp(phi_sim_mix[1,1:500]))
mean(exp(phi_sim_mix[2,1:500]))
var(exp(phi_sim_mix[1,1:500]))
var(exp(phi_sim_mix[2,1:500]))
mean((exp(phi_sim_mix[1,1:500])^2)*phi_avar_mix[1,1:500])
mean((exp(phi_sim_mix[2,1:500])^2)*phi_avar_mix[2,1:500])
max(exp(phi_sim_mix[1,1:500]))
min(exp(phi_sim_mix[1,1:500]))
max(exp(phi_sim_mix[2,1:500]))
min(exp(phi_sim_mix[2,1:500]))

#sigma
mean(exp(sg_sim_mix[1,1:500]))
mean(exp(sg_sim_mix[2,1:500]))
var(exp(sg_sim_mix[1,1:500]))
var(exp(sg_sim_mix[2,1:500]))
mean((exp(sg_sim_mix[1,1:500])^2)*sg_avar_mix[1,1:500])
mean((exp(sg_sim_mix[2,1:500])^2)*sg_avar_mix[2,1:500])
max(exp(sg_sim_mix[1,1:500]))
min(exp(sg_sim_mix[1,1:500]))
max(exp(sg_sim_mix[2,1:500]))
min(exp(sg_sim_mix[2,1:500]))

#lo
mean(sg_sim_mix[3,1:500])
var(sg_sim_mix[3,1:500])
mean(sg_avar_mix[3,1:500])
mean((sg_sim_mix[3,1:500]/sqrt(1+sg_sim_mix[3,1:500]^2)))
var((sg_sim_mix[3,1:500]/sqrt(1+sg_sim_mix[3,1:500]^2)))
sqrt(var((sg_sim_mix[3,1:500]/sqrt(1+sg_sim_mix[3,1:500]^2))))
max(sg_sim_mix[3,1:500])
min(sg_sim_mix[3,1:500])
max((sg_sim_mix[3,1:500]/sqrt(1+sg_sim_mix[3,1:500]^2)))
min((sg_sim_mix[3,1:500]/sqrt(1+sg_sim_mix[3,1:500]^2)))
#------------------------------------------------------------------------------

