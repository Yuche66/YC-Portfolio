#==============================================================================
#-------------simulation part 2-1 (Beta&UG model misspecified)-----------------
#==============================================================================

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
    term1=log(mu)-phi*log(1-mu^(1/phi))-log(factorial(phi-1))
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
#------------------------------------------------------------------------------


#start-------------------------------------------------------------------------
#simulation times N=3000, sample size n=200,400
#beta_t=(0,1,-0.5) XXX~N(3,0.7^2) & XXX~Ber(0.5)
#beta model aphb=(25,15), unit gamma model aphug=(3,5)
#random effect variance matrix sg=0.1, rho=0.5
#------------------------------------------------------------------------------
#setting
N=3000
n=400
dm=2
beta_t=c(0,1,-0.5)
p=length(beta_t)-1
YYY_bt=matrix(0,n,dm)
YYY_ug=matrix(0,n,dm)
aphb=c(25,15)
aphug=c(3,5)

#random effect setting
rdefflo=c(0.5)
Crm=matrix(1,dm,dm)
Crm[upper.tri(Crm)]=rdefflo
Crm[lower.tri(Crm)]=rdefflo
G=diag(rep(0.1,dm),dm,dm)%*%Crm%*%diag(rep(0.1,dm),dm,dm)
L=t(chol(G))

#mu setting
set.seed(923)
XXX1=matrix(rnorm(n*dm,3,0.7),n*dm,1)
XXX2=matrix(rbinom(n*dm,1,0.5),n*dm,1)
XXX=cbind(XXX1,XXX2)
link=function(x){
    1-1/(1+exp(x))
}
dlink=function(x){
    exp(x)/((1+exp(x))^2)
}
it=XXX%*%matrix(beta_t[2:(p+1)],p,1)+beta_t[1]
mufix=XXX%*%matrix(beta_t[2:(p+1)],p,1)+beta_t[1]
mufix=sapply(mufix,link)
matrix(mufix,n,dm,byrow=T)

#scatter
gd1=seq(beta_t[1]-0,beta_t[1]+0,by=0.4)
gd2=seq(beta_t[2]-0,beta_t[2]+0,by=0.4)
gd3=seq(beta_t[3]-0,beta_t[3]+0,by=0.4)
gd=c()
gdn=length(gd1)
gdn3=gdn*gdn*gdn
gdtemp=as.matrix(expand.grid(gd1,gd2))
gdtemp=t(gdtemp)
for(w in c(1:gdn)){
    w1=rep(gd3[w],gdn)
    gd=cbind(gd,rbind(gdtemp,w1))
}

#record the result
phi_sim_bt=matrix(0,dm,N)
phi_avar_bt=matrix(0,dm,N)
beta_sim_bt=matrix(0,p+1,N)
beta_dllhvalue_bt=matrix(0,N)
beta_avar_bt=matrix(0,p+1,N)
walds_beta_bt=matrix(0,p+1,N)

phi_sim_ug=matrix(0,dm,N)
phi_avar_ug=matrix(0,dm,N)
beta_sim_ug=matrix(0,p+1,N)
beta_dllhvalue_ug=matrix(0,N)
beta_avar_ug=matrix(0,p+1,N)
walds_beta_ug=matrix(0,p+1,N)

#data come from unit gamma (working model is beta)-----------------------------
bt_update=0
for(M in c(1:N)){
    
    #generating data(UG)
    #n=200:Seed:619+210M, n=400:Seed:619M+210
    set.seed(619+210*M)
    for(i in c(1:n)){
        rdeff=rmvnorm(n,c(0,0),G)
        mu=sapply(it+c(t(rdeff)),link)
        mu_i=mu[(dm*i-dm+1):(dm*i)]
        cd_i=runif(dm,0,1)
        rug_i=rep(0,dm)
        for(ii in c(1:dm)){
            rug_i[ii]=uniroot(function(y){unit_gamma_cdf(aphug[ii],mu_i[ii],y)-cd_i[ii]},
                              interval = c(0,1))$root
            if(rug_i[ii]>=0.9999999){rug_i[ii]=0.9999999}
            if(rug_i[ii]<=0.0000001){rug_i[ii]=0.0000001}
        }
        YYY_ug[i,]=rug_i
    }
    YYY_bt=YYY_ug
    
    #laplace approximation for likelihood function(Beta)
    loglike_b_2d_laplace=function(para){
        hG=matrix(0,2,2)
        diag(hG)=para[6:7]
        hG[lower.tri(hG)]=para[8]
        G=hG%*%t(hG)
        invG=solve(G)
        marg_loglike=rep(0,n)
        itavalue=XXX%*%matrix(para[4:5],p,1)+para[3]
        itavalue=matrix(c(itavalue),length(YYY_bt[,1]),dm,byrow = T)
        for(i in c(1:n)){
            f1=function(u){
                mub1=sapply((itavalue[i,1]+u[1]),link)
                mub2=sapply((itavalue[i,2]+u[2]),link)
                a1=mub1*exp(para[1])
                a2=mub2*exp(para[2])
                b1=(1-mub1)*exp(para[1])
                b2=(1-mub2)*exp(para[2])
                uu=rbind(u[1],u[2])
                logprob=log_dbeta(YYY_bt[i,1],a1,b1)+log_dbeta(YYY_bt[i,2],a2,b2)
                logprobjoint=logprob+log(((2*pi)^(-dm/2))*(det(G)^(-1/2)))+
                    (-(1/2)*diag(t(uu)%*%invG%*%uu))
                return(logprobjoint)
            }
            optf1=optim(c(0,0),function(u){-f1(u)},hessian = T)
            hessf1=optf1$hessian
            term1=(-1/2)*log(det(hessf1))
            marg_loglike[i]=(term1+f1(optf1$par))+(dm/2)*log(2*pi)
        }
        ans=sum(marg_loglike)
        return(ans)
    }
    
    #initial point
    loglikevalue=rep(0,gdn3)
    for(w in c(1:gdn3)){
        loglikevalue[w]=-loglike_b_2d_laplace(c(log(aphb),as.numeric(gd[,w]),L[1,1],L[2,2],L[2,1]))
        print(c(M,"bt_regression_scatting",w))
    }
    optinit=gd[,order(loglikevalue)[1]]
    optinit=as.numeric(optinit)
    optinit=c(log(6),log(10),optinit)
    
    #optimize
    opt_bt=optim(optinit,function(para){-loglike_b_2d_laplace(c(para,L[1,1],L[2,2],L[2,1]))},
                 method="Nelder-Mead",hessian = T)
    opt_bt_best=opt_bt
    
    #To ensure it is not saddle
    #check all eigenvalues of the hessian matrix need to be positive(positive-define)
    egvalue_bt=eigen(opt_bt$hessian)$values
    detect=length(egvalue_bt[egvalue_bt<0])
    #if the matrix is not positive-define, repeat optim()
    if(detect>0){
        bt_update=bt_update+1
        update_times=0
        repeat{
            inpar=0.8*opt_bt$par+0.2*optinit
            r1=runif(1,inpar[1]-0.25,inpar[1]+0.25)
            r2=runif(1,inpar[2]-0.25,inpar[2]+0.25)
            r3=runif(1,inpar[3]-0.25,inpar[3]+0.25)
            r4=runif(1,inpar[4]-0.25,inpar[4]+0.25)
            r5=runif(1,inpar[5]-0.25,inpar[5]+0.25)
            init=c(r1,r2,r3,r4,r5)
            opt_bt_best=optim(init,function(para){-loglike_b_2d_laplace(c(para,L[1,1],L[2,2],L[2,1]))},
                              method="Nelder-Mead",hessian = T)
            egvalue_bt=eigen(opt_bt_best$hessian)$values
            detect=length(egvalue_bt[egvalue_bt<0]) 
            update_times=update_times+1
            print(c(M,"bt_regression_updating",update_times))
            if(detect==0 | update_times==15)break
        }
    }
    bt_varhat_best=solve(opt_bt_best$hessian) 
    beta_sim_bt[1:3,M]=opt_bt_best$par[3:5]
    beta_avar_bt[1:3,M]=diag(bt_varhat_best)[3:5]
    beta_dllhvalue_bt[M]=opt_bt$value-opt_bt_best$value
    phi_sim_bt[1:2,M]=exp(opt_bt_best$par[1:2])
    phi_avar_bt[1:2,M]=exp(opt_bt_best$par[1:2])^2*diag(bt_varhat_best)[1:2]
    
    #wald test
    walds_beta_bt[1:3,M]=((beta_sim_bt[1:3,M]-beta_t)^2)/beta_avar_bt[1:3,M]
    
    print(c(M,"bt_regression complete"))
}
#------------------------------------------------------------------------------


#data come from beta (working model is unit gamma)-----------------------------
YYY_bt=matrix(0,n,dm)
YYY_ug=matrix(0,n,dm)
ug_update=0
for(M in c(1:N)){
    
    #generating data
    #n=200:Seed:619+210M, n=400:Seed:619+210M
    set.seed(619+210*M)
    for(i in c(1:n)){
        rdeff=rmvnorm(n,c(0,0),G)
        mu=sapply(it+c(t(rdeff)),link)
        mu_i=mu[(dm*i-dm+1):(dm*i)]
        rbeta_i=rep(0,dm)
        for(ii in c(1:dm)){
            rbeta_i[ii]=rbeta(1,mu_i[ii]*aphb[ii],(1-mu_i[ii])*aphb[ii])
            if(rbeta_i[ii]>=0.9999999){rbeta_i[ii]=0.9999999}
            if(rbeta_i[ii]<=0.0000001){rbeta_i[ii]=0.0000001}
        }
        YYY_bt[i,]=rbeta_i
    }
    YYY_ug=YYY_bt
    
    #laplace approximation for likelihood function(UG)
    loglike_ug_2d_laplace=function(para){
        hG=matrix(0,2,2)
        diag(hG)=para[6:7]
        hG[lower.tri(hG)]=para[8]
        G=hG%*%t(hG)
        invG=solve(G)
        marg_loglike=0
        itavalue=XXX%*%matrix(para[4:5],p,1)+para[3]
        itavalue=matrix(c(itavalue),length(YYY_ug[,1]),dm,byrow = T)
        for(i in c(1:n)){
            f1=function(u){
                muug1=sapply((itavalue[i,1]+u[1]),link)
                muug2=sapply((itavalue[i,2]+u[2]),link)
                uu=rbind(u[1],u[2])
                logprob=log_unit_gamma_pdf(exp(para[1]),muug1,YYY_ug[i,1])+
                    log_unit_gamma_pdf(exp(para[2]),muug2,YYY_ug[i,2])
                logprobjoint=logprob+log(((2*pi)^(-dm/2))*(det(G)^(-1/2)))+
                    (-(1/2)*diag(t(uu)%*%invG%*%uu))
                return(logprobjoint)
            }
            optf1=optim(c(0,0),function(u){-f1(u)},hessian = T)
            hessf1=optf1$hessian
            term1=(-1/2)*log(det(hessf1))
            marg_loglike=marg_loglike+(term1+f1(optf1$par))+dm/2*log(2*pi)
        }
        ans=marg_loglike
        return(ans)
    }
    
    #initial point
    loglikevalue=rep(0,gdn3)
    for(w in c(1:gdn3)){
        loglikevalue[w]=-loglike_ug_2d_laplace(c(log(aphug),as.numeric(gd[,w]),L[1,1],L[2,2],L[2,1]))
        print(c(M,"ug_regression_scatting",w))
    }
    optinit=gd[,order(loglikevalue)[1]]
    optinit=as.numeric(optinit)
    optinit=c(log(5),log(2),optinit)
    
    #optimize
    opt_ug=optim(optinit,function(para){-loglike_ug_2d_laplace(c(para,L[1,1],L[2,2],L[2,1]))},
                 method="Nelder-Mead",hessian = T)
    opt_ug_best=opt_ug
    
    #To ensure it is not saddle
    #check all eigenvalues of the hessian matrix need to be positive(positive-define)
    egvalue_ug=eigen(opt_ug$hessian)$values
    detect=length(egvalue_ug[egvalue_ug<0])
    #if the matrix is not positive-define, repeat optim()
    if(detect>0){
        ug_update=ug_update+1
        update_times=0
        repeat{
            inpar=0.8*opt_ug$par+0.2*optinit
            g1=runif(1,inpar[1]-0.25,inpar[1]+0.25)
            g2=runif(1,inpar[2]-0.25,inpar[2]+0.25)
            g3=runif(1,inpar[3]-0.25,inpar[3]+0.25)
            g4=runif(1,inpar[4]-0.25,inpar[4]+0.25)
            g5=runif(1,inpar[5]-0.25,inpar[5]+0.25)
            init=c(g1,g2,g3,g4,g5)
            opt_ug_best=optim(init,function(para){-loglike_ug_2d_laplace(c(para,L[1,1],L[2,2],L[2,1]))},
                              method="Nelder-Mead",hessian = T)
            egvalue_ug=eigen(opt_ug_best$hessian)$values
            detect=length(egvalue_ug[egvalue_ug<0])
            update_times=update_times+1
            print(c(M,"ug_regression_updating",update_times))
            if(detect==0 | update_times==15)break
        }
    }
    ug_varhat_best=solve(opt_ug_best$hessian)
    beta_sim_ug[1:3,M]=opt_ug_best$par[3:5]
    beta_avar_ug[1:3,M]=diag(ug_varhat_best)[3:5]
    beta_dllhvalue_ug[M]=opt_ug$value-opt_ug_best$value
    phi_sim_ug[1:2,M]=exp(opt_ug_best$par[1:2])
    phi_avar_ug[1:2,M]=exp(opt_ug_best$par[1:2])^2*diag(ug_varhat_best)[1:2]
    
    #wald test
    walds_beta_ug[1:3,M]=((beta_sim_ug[1:3,M]-beta_t)^2)/beta_avar_ug[1:3,M]
    print(c(M,"ug_regression_complete"))
}
#------------------------------------------------------------------------------


#beta model result(data comes from UG)-----------------------------------------
beta_sim_bt
beta_dllhvalue_bt
beta_avar_bt
bt_update
#mean, sample variance, variance hat result
mean(beta_sim_bt[1,])
mean(beta_sim_bt[2,])
mean(beta_sim_bt[3,])
var(beta_sim_bt[1,])
var(beta_sim_bt[2,])
var(beta_sim_bt[3,])
mean(beta_avar_bt[1,])
mean(beta_avar_bt[2,])
mean(beta_avar_bt[3,])

#phi result
mean(phi_sim_bt[1,])
mean(phi_sim_bt[2,])
var(phi_sim_bt[1,])
var(phi_sim_bt[2,])
mean(phi_avar_bt[1,])
mean(phi_avar_bt[2,])

#wald test result
reject_W_bt=+(walds_beta_bt>qchisq(0.95,1))
mean(reject_W_bt[1,])
mean(reject_W_bt[2,])
mean(reject_W_bt[3,])
#------------------------------------------------------------------------------


#UG model result(data comes from beta)-----------------------------------------
beta_sim_ug
beta_dllhvalue_ug
beta_avar_ug
ug_update
#mean, sample variance, variance hat result
mean(beta_sim_ug[1,])
mean(beta_sim_ug[2,])
mean(beta_sim_ug[3,])
var(beta_sim_ug[1,])
var(beta_sim_ug[2,])
var(beta_sim_ug[3,])
mean(beta_avar_ug[1,])
mean(beta_avar_ug[2,])
mean(beta_avar_ug[3,])

#phi result
mean(phi_sim_ug[1,])
mean(phi_sim_ug[2,])
var(phi_sim_ug[1,])
var(phi_sim_ug[2,])
mean(phi_avar_ug[1,])
mean(phi_avar_ug[2,])

#wald test result
reject_W_ug=+(walds_beta_ug>qchisq(0.95,1))
mean(reject_W_ug[1,])
mean(reject_W_ug[2,])
mean(reject_W_ug[3,])
#------------------------------------------------------------------------------


