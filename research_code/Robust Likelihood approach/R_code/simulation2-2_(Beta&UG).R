#==============================================================================
#-------------------simulation part 2-2 (Beta& UG model)-----------------------
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
#beta_t=(-1,0.1,0.1), XXX~N(0,0.3^2) & XXX~Ber(0.5)
#beta model aphb=(9,4), unit gamma model aphug=(3,5)
#random effect(not follow normal dist.) variance matrix sg=0.3, rho=0.9 
#------------------------------------------------------------------------------
#setting
N=3000
n=200
dm=2
beta_t=c(-1,0.1,0.1)
p=length(beta_t)-1
YYY_bt=matrix(0,n,dm)
YYY_ug=matrix(0,n,dm)
aphb=c(9,4)
aphug=c(3,5)

#random effect setting
#generate correlated V~Beta(0.01,0.01) Var(V)=1/(4*(1+0.01+0.01))
#let r=sqrt(1/(100/9*Var(V))) => u=(V-0.5)*r, E(u)=0, Var(u)=0.09
#given correlation rho=0.9
rdefflo=c(0.9)
ap=c(0.01,0.01)
bt=c(0.01,0.01)
mn=ap/(ap+bt)
rr=1/(4*(1+0.01+0.01))
vr=ap*bt/(((ap+bt)^2)*(ap+bt+1))/(100/9*rr)
Crm=matrix(1,dm,dm)
Crm[upper.tri(Crm)]=rdefflo
Crm[lower.tri(Crm)]=rdefflo
G=diag(rep(sqrt(vr),dm),dm,dm)%*%Crm%*%diag(rep(sqrt(vr),dm),dm,dm)
L=t(chol(G))

#generate correlated Beta dist. by normal copula
#solve the parameter(omega) in normal copula-----------------------------------
normal_cupcor=rep(0,1)
count_cup=0
for(ii in c(2:dm)){
    for(jj in c(1:(ii-1))){
        count_cup=count_cup+1
        btlo_fun=function(lo){
            f1=function(x){
                f2=function(y){
                    qbeta(pnorm(x),ap[jj],bt[jj])*qbeta(pnorm(y),ap[ii],bt[ii])*
                        1/sqrt(1-lo^2)*1/(2*pi)*exp(-1/(2*(1-lo^2))*(x^2-2*lo*x*y+y^2))
                }
                ans=integrate(f2,-8,8)$value
                return(ans)
            }
            v=mn[ii]
            u=mn[jj]
            #calculate covariance by double integral
            #the first is in f1 & the second is Reimann integral in ans2
            ans2=sum(sapply(seq(-8,8,by=0.1),f1)*0.1)-u*v
            #Let ans3=0 and we will get the parameters in normal copula
            ans3=ans2/sqrt((100/9*rr)^2*vr[jj]*vr[ii])-rdefflo[count_cup]
            return(ans3)
        }
        #solve the root of ans3=0
        normal_cupcor[count_cup]=uniroot(btlo_fun,interval=c(-0.99,0.99))$root
    }
}

#check the geremating method is working
set.seed(619)
btcopula=normalCopula(normal_cupcor,dim=dm,dispstr = "un")
mvbt=mvdc(btcopula,margins=c("beta","beta"),paramMargins = list(list(ap[1],bt[1]),list(ap[2],bt[2])))
rdeff=rMvdc(n,mvbt)
for(k in c(1:dm)){
    rdeff[,k]=(rdeff[,k]-mn[k])/sqrt((100/9*rr))
}
mean(rdeff[,1])
mean(rdeff[,2])
var(rdeff[,1])
var(rdeff[,2])
vr
cor(rdeff[,1],rdeff[,2])
plot(rdeff[,1],rdeff[,2],xlim=c(-0.4,0.4),ylim=c(-0.4,0.4))
points(rdeff[,1],rdeff[,2],xlim=c(-0.4,0.4),ylim=c(-0.4,0.4),col="red",pch=4)
#------------------------------------------------------------------------------

#mu setting
set.seed(923)
XXX1=matrix(runif(n*dm,-1,1),n*dm,1)
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
gd1=seq(beta_t[1]-0,beta_t[1]+0,by=0.6)
gd2=seq(beta_t[2]-0,beta_t[2]+0,by=0.6)
gd3=seq(beta_t[3]-0,beta_t[3]+0,by=0.6)
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

#beta regression---------------------------------------------------------------
bt_update=0
for(M in c(1:N)){
    
    #generating data
    #n=200:Seed:619+210M, n=400:Seed:619+210M
    set.seed(619+210*M)
    for(i in c(1:n)){
        rdeff=rMvdc(n,mvbt)
        for(k in c(1:dm)){
            rdeff[,k]=(rdeff[,k]-mn[k])/sqrt((100/9*rr))
        }
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
    
    #laplace approximation for likelihood function
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
                a1=mub1*para[1]
                a2=mub2*para[2]
                b1=(1-mub1)*para[1]
                b2=(1-mub2)*para[2]
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
        loglikevalue[w]=-loglike_b_2d_laplace(c(aphb,as.numeric(gd[,w]),L[1,1],L[2,2],L[2,1]))
        print(c(M,"bt_regression_scatting",w))
    }
    optinit=gd[,order(loglikevalue)[1]]
    optinit=as.numeric(optinit)
    
    #optimize
    opt_bt=optim(optinit,function(para){-loglike_b_2d_laplace(c(aphb,para,L[1,1],L[2,2],L[2,1]))},
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
            r4=runif(1,inpar[1]-0.25,inpar[1]+0.25)
            r5=runif(1,inpar[2]-0.25,inpar[2]+0.25)
            r6=runif(1,inpar[3]-0.1,inpar[3]+0.1)
            init=c(r4,r5,r6)
            opt_bt_best=optim(init,function(para){-loglike_b_2d_laplace(c(aphb,para,L[1,1],L[2,2],L[2,1]))},
                              method="Nelder-Mead",hessian = T)
            egvalue_bt=eigen(opt_bt_best$hessian)$values
            detect=length(egvalue_bt[egvalue_bt<0]) 
            update_times=update_times+1
            print(c(M,"bt_regression_updating",update_times))
            if(detect==0 | update_times==15)break
        }
    }
    bt_varhat_best=solve(opt_bt_best$hessian) 
    beta_sim_bt[1:3,M]=opt_bt_best$par
    beta_avar_bt[1:3,M]=diag(bt_varhat_best)
    beta_dllhvalue_bt[M]=opt_bt$value-opt_bt_best$value
    
    #wald test
    walds_beta_bt[1:3,M]=((beta_sim_bt[1:3,M]-beta_t)^2)/beta_avar_bt[1:3,M]
    print(c(M,"bt_regression complete"))
}
#------------------------------------------------------------------------------


#unit gamma regression---------------------------------------------------------
ug_update=0
for(M in c(1:N)){
    
    #generating data
    #n=200:Seed:619+210M, n=400:Seed:619M+210
    set.seed(619+210*M)
    for(i in c(1:n)){
        rdeff=rMvdc(n,mvbt)
        for(k in c(1:dm)){
            rdeff[,k]=(rdeff[,k]-mn[k])/sqrt((100/9*rr))
        }
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
    
    #laplace approximation for likelihood function
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
                logprob=log_unit_gamma_pdf(para[1],muug1,YYY_ug[i,1])+
                    log_unit_gamma_pdf(para[2],muug2,YYY_ug[i,2])
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
        loglikevalue[w]=-loglike_ug_2d_laplace(c(aphug,as.numeric(gd[,w]),L[1,1],L[2,2],L[2,1]))
        print(c(M,"ug_regression_scatting",w))
    }
    optinit=gd[,order(loglikevalue)[1]]
    optinit=as.numeric(optinit)
    
    #optimize
    opt_ug=optim(optinit,function(para){-loglike_ug_2d_laplace(c(aphug,para,L[1,1],L[2,2],L[2,1]))},
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
            g4=runif(1,inpar[1]-0.25,inpar[1]+0.25)
            g5=runif(1,inpar[2]-0.25,inpar[2]+0.25)
            g6=runif(1,inpar[3]-0.1,inpar[3]+0.1)
            init=c(g4,g5,g6)
            opt_ug_best=optim(init,function(para){-loglike_ug_2d_laplace(c(aphug,para,L[1,1],L[2,2],L[2,1]))},
                              method="Nelder-Mead",hessian = T)
            egvalue_ug=eigen(opt_ug_best$hessian)$values
            detect=length(egvalue_ug[egvalue_ug<0])
            update_times=update_times+1
            print(c(M,"ug_regression_updating",update_times))
            if(detect==0 | update_times==15)break
        }
    }
    ug_varhat_best=solve(opt_ug_best$hessian)
    beta_sim_ug[1:3,M]=opt_ug_best$par
    beta_avar_ug[1:3,M]=diag(ug_varhat_best)
    beta_dllhvalue_ug[M]=opt_ug$value-opt_ug_best$value
    
    #wald test
    walds_beta_ug[1:3,M]=((beta_sim_ug[1:3,M]-beta_t)^2)/beta_avar_ug[1:3,M]
    print(c(M,"ug_regression_complete"))
}
#------------------------------------------------------------------------------


#beta result-------------------------------------------------------------------
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
#wald test result
reject_W_bt=+(walds_beta_bt>qchisq(0.95,1))
mean(reject_W_bt[1,])
mean(reject_W_bt[2,])
mean(reject_W_bt[3,])
#------------------------------------------------------------------------------

#unit gamma result-------------------------------------------------------------
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
#wald test result
reject_W_ug=+(walds_beta_ug>qchisq(0.95,1))
mean(reject_W_ug[1,])
mean(reject_W_ug[2,])
mean(reject_W_ug[3,])
#------------------------------------------------------------------------------


