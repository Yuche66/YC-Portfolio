#==============================================================================
#---------------------simulation part 5-2 (MNB model)--------------------------
#==============================================================================

#package-----------------------------------------------------------------------
library(numDeriv)
library(copula)
library(MASS)
library(pracma)
library(mvtnorm)
#------------------------------------------------------------------------------

#MNB log likelihood func., score func., I, V, A, B-----------------------------
mvnb_loglike_reg_rem=function(y_m,x_m,beta,link_f,phi){
    Y=y_m
    m=length(y_m[,1])
    p=length(x_m[1,])
    mn=length(x_m[,1])
    X=cbind(matrix(1,mn,1),x_m)
    ita=X%*%matrix(beta,p+1,1)
    mu=sapply(ita,link_f)
    detect= +(mu<=0)
    if(sum(detect)>0){
        loglike=-Inf
    }else{
        loglike=matrix(0,1,m)
        nonaY=c(na.omit(c(t(Y))))
        temp1=sum(nonaY*log(c(mu)))
        nn=0
        for(i in c(1:m)){
            YY=c(na.omit(c(Y[i,])))
            Yiplus=sum(YY)
            muiplus=0
            for(j in c((nn+1):(nn+length(YY)))){
                muiplus=muiplus+c(mu)[j]
            }
            loglike[,i]=-(1/phi+Yiplus)*log(1+muiplus*phi)
            nn=nn+length(YY)
        }
        loglike=sum(loglike)+temp1
    }
    return(c(loglike))
}
mvnb_score_beta_rem=function(y_m,x_m,beta,link_f,dlink_f,phi,ith){
    Y=y_m
    m=length(y_m[,1])
    p=length(x_m[1,])
    mn=length(x_m[,1])
    X=cbind(matrix(1,mn,1),x_m)
    ita=X%*%matrix(beta,p+1,1)
    mu=sapply(ita,link_f)
    gp=sapply(ita,dlink_f)
    nonaY=c(na.omit(c(t(Y))))
    temp1=sum((nonaY-c(mu))/c(mu)*c(X[,ith])*gp)
    nn=0
    score_f=matrix(0,1,m)
    for(i in c(1:m)){
        YY=c(na.omit(c(Y[i,])))
        Yiplus=sum(YY)
        muiplus=0
        giplus=0
        for(j in c((nn+1):(nn+length(YY)))){
            muiplus=muiplus+c(mu)[j]
            giplus=giplus+(c(gp)*c(X[,ith]))[j]
        }
        score_f[,i]=phi*(Yiplus-muiplus)/(1+phi*muiplus)*(giplus)
        nn=nn+length(YY)
    }
    score_f=-sum(score_f)+temp1
    return(score_f)
}
mvnb_fisher_reg_rem=function(y_m,x_m,beta,link_f,dlink_f,phi){
    Y=y_m
    m=length(y_m[,1])
    p=length(x_m[1,])
    mn=length(x_m[,1])
    X=cbind(matrix(1,mn,1),x_m)
    ita=X%*%matrix(beta,p+1,1)
    mu=sapply(ita,link_f)
    gp=sapply(ita,dlink_f)
    l0=t(X)%*%diag(c(gp),mn,mn)%*%diag(c(1/sqrt(mu)),mn,mn)
    l1=matrix(0,p+1,m)
    temp=t(X)%*%diag(c(gp),mn,mn)
    muip=rep(0,m)
    nn=0
    for(i in c(1:m)){
        YY=c(na.omit(c(Y[i,])))
        add=0
        muiplus=0
        for(j in c((nn+1):(nn+length(YY)))){
            muiplus=muiplus+c(mu)[j]
            add=add+temp[,j]
        }
        l1[,i]=add
        muip[i]=muiplus
        nn=nn+length(YY)
    }
    
    Fish_m=(l0%*%t(l0)-l1%*%diag(phi/(1+phi*muip),m,m)%*%t(l1))/m
    return(Fish_m)
}
mvnb_scorevar_reg_rem=function(y_m,x_m,beta,link_f,dlink_f,phi){
    Y=y_m
    m=length(y_m[,1])
    p=length(x_m[1,])
    mn=length(x_m[,1])
    X=cbind(matrix(1,mn,1),x_m)
    ita=X%*%matrix(beta,p+1,1)
    mu=sapply(ita,link_f)
    gp=sapply(ita,dlink_f)
    nonaY=c(na.omit(c(t(Y))))
    temp=nonaY/mu
    V0=matrix(0,m,p+1)
    yplus=matrix(0,m,1)
    muiplus=matrix(0,m,1)
    nn=0
    for(i in c(1:m)){
        YY=c(na.omit(c(Y[i,])))
        a=0
        b=0
        for(j in c((nn+1):(nn+length(YY)))){
            a=a+c(mu)[j]
            b=b+nonaY[j]
        }
        yplus[i,]=b
        muiplus[i,]=a
        nn=nn+length(YY)
    }
    nn=0
    for(i in c(1:m)){
        YY=c(na.omit(c(Y[i,])))
        add=0
        for(j in c((nn+1):(nn+length(YY)))){
            add=add+(temp[j]-(1+phi*yplus[i,])/(1+phi*muiplus[i,]))*c(gp)[j]*X[j,]
        }
        V0[i,]=add
        nn=nn+length(YY)
    }
    scoreV=var(V0)
    return(scoreV)
}
mvnb_A_reg_rem=function(y_m,x_m,beta,link_f,dlink_f,phi,ith){
    Fish_m=mvnb_fisher_reg_rem(y_m,x_m,beta,link_f,dlink_f,phi)
    I1=Fish_m[ith,ith]
    I2=Fish_m[-ith,ith]
    I3=Fish_m[-ith,-ith]
    A=I1-t(I2)%*%solve(I3)%*%I2
    return(A)
}
mvnb_B_reg_rem=function(y_m,x_m,beta,link_f,dlink_f,phi,ith){
    Fish_m=mvnb_fisher_reg_rem(y_m,x_m,beta,link_f,dlink_f,phi)
    scoreV=mvnb_scorevar_reg_rem(y_m,x_m,beta,link_f,dlink_f,phi)
    V1=scoreV[ith,ith]
    I2=Fish_m[-ith,ith]
    V2=scoreV[-ith,ith]
    I3=Fish_m[-ith,-ith]
    V3=scoreV[-ith,-ith]
    temp=t(I2)%*%solve(I3)
    B=V1-2*t(I2)%*%solve(I3)%*%V2+temp%*%V3%*%t(temp)
    return(B)
}
#beta & unit gamma pdf, cdf
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
#simulation times N=3000, sample size n=50,100,200,400
#beta_t=(-2,3,0) & XXX1~Unif(0,1), XXX2~Unif(0,1),
#beta model aphb=(9,4,1.5),
#random effect variance matrix sigma=1, rho12=0.6, rho13=0.4, rho23=0.1,
#------------------------------------------------------------------------------
#setting
N=3000
n=400
dm=3
beta_t=c(-2,3,0)
#the approximation of regression (Zeger)
c=16*sqrt(3)/(15*pi)
sg=1
beta_adj=(((c^2)*(sg^2)+1)^(-1/2))*beta_t
p=length(beta_t)-1
YYY_bt=matrix(0,n,dm)
YYY_ug=matrix(0,n,dm)
aphb=c(9,4,1.5)

#random effect setting
rdefflo=c(0.6,0.4,0.1)
Crm=matrix(1,dm,dm)
Crm[upper.tri(Crm)]=rdefflo
Crm[lower.tri(Crm)]=rdefflo
G=diag(rep(sg,dm),dm,dm)%*%Crm%*%diag(rep(sg,dm),dm,dm)
L=t(chol(G))

#mu setting
set.seed(923)
XXX1=matrix(runif(n*dm,0,1),n*dm,1)
XXX2=matrix(runif(n*dm,0,1),n*dm,1)
XXX3=matrix(rbinom(n*dm,1,0.5),n*dm,1)
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

#nuisance parameter phi=1
h=0.02*c(1:5000)
num=10
ph=1

#record the result
beta_sim=matrix(0,p+1,N)
betavar_nai=matrix(0,p+1,N)
betavar_rb=matrix(0,p+1,N)
walds_beta_rb=matrix(0,p+1,N)
walds_beta_nai=matrix(0,p+1,N)
score_beta_rb=matrix(0,p+1,N)
score_beta_nai=matrix(0,p+1,N)
LR_beta_rb=matrix(0,p+1,N)
LR_beta_nai=matrix(0,p+1,N)
walds_allpara=matrix(0,2,N)
score_allpara=matrix(0,2,N)

#data comes from beta----------------------------------------------------------
times=0
set.seed(621019)
for(M in c(1:N)){
    
    #generating data
    for(i in c(1:n)){
        rdeff=rmvnorm(n,c(0,0,0),G)
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
    
    #log likelihood function
    max_f=function(beta){
        -mvnb_loglike_reg_rem(YYY_bt,XXX,beta,link,ph)
    }
    
    #initial point
    v=rep(0,(p+1))
    for(j in c(1:(p+1))){
        v[j]=runif(1,beta_t[j]-0.5,beta_t[j]+0.5)
    }
    
    #optimize
    opt=optim(v,max_f)
    beta_hat=opt$par
    lr=opt$value
    beta_sim[,M]=beta_hat
    sf_all=matrix(0,(p+1),1)
    
    #I hat, V hat, naive variance hat, robust variance hat
    I=mvnb_fisher_reg_rem(YYY_bt,XXX,beta_hat,link,dlink,ph)
    V=mvnb_scorevar_reg_rem(YYY_bt,XXX,beta_hat,link,dlink,ph)
    var_na=solve(I)/n
    var_rb=solve(I)%*%V%*%solve(I)/n
    
    #three kinds of tests
    for(q in c(1:(p+1))){
        betavar_nai[q,M]=var_na[q,q]
        betavar_rb[q,M]=var_rb[q,q]
        A=mvnb_A_reg_rem(YYY_bt,XXX,beta_hat,link,dlink,ph,q)
        B=mvnb_B_reg_rem(YYY_bt,XXX,beta_hat,link,dlink,ph,q)
        
        #LR test 
        max_f0=function(beta0){
            beta=rep(beta_adj[q],(p+1))
            beta[-q]=beta0
            ans=-mvnb_loglike_reg_rem(YYY_bt,XXX,beta,link,ph)
            return(ans)
        }
        opt0=optim(v[-q],max_f0)
        beta_0_hat=opt0$par
        lr0=opt0$value
        LR_beta_nai[q,M]=2*(lr0-lr)
        LR_beta_rb[q,M]=2*(A/B)*(lr0-lr)
        
        #walds test
        w1=A%*%matrix(beta_hat[q]-beta_adj[q],1,1)
        walds_beta_rb[q,M]=n*t(w1)%*%solve(B)%*%w1
        walds_beta_nai[q,M]=n*t(w1)%*%solve(A)%*%w1
        
        #score test
        beta=rep(beta_adj[q],(p+1))
        beta[-q]=beta_0_hat
        sA=mvnb_A_reg_rem(YYY_bt,XXX,beta,link,dlink,ph,q)
        sB=mvnb_B_reg_rem(YYY_bt,XXX,beta,link,dlink,ph,q)
        s1=mvnb_score_beta_rem(YYY_bt,XXX,beta,link,dlink,ph,q)
        score_beta_rb[q,M]=1/n*t(s1)%*%solve(sB)%*%s1
        score_beta_nai[q,M]=1/n*t(s1)%*%solve(sA)%*%s1
        
        #score(all) test
        sf_all[q,]=mvnb_score_beta_rem(YYY_bt,XXX,beta_adj,link,dlink,ph,q)
    }
    
    #robust walds test(all)
    W=I%*%matrix(beta_hat-beta_adj,p+1,1)
    walds_allpara[1,M]=n*t(W)%*%solve(V)%*%W
    
    #naive walds test(all)
    W=I%*%matrix(beta_hat-beta_adj,p+1,1)
    walds_allpara[2,M]=n*t(W)%*%solve(I)%*%W
    
    #robust score test(all)
    S=sf_all
    score_allpara[1,M]=1/n*t(S)%*%solve(V)%*%S
    
    #naive score test(all)
    S=sf_all
    score_allpara[2,M]=1/n*t(S)%*%solve(I)%*%S
    
    times=times+1
    print(times)
}
#------------------------------------------------------------------------------

#data comes from beta result---------------------------------------------------
#result(mean, sample variance, variance hat)
table_reg_NB1=matrix(c(mean(beta_sim[1,]),mean(beta_sim[2,]),mean(beta_sim[3,]),
                       var(beta_sim[1,]),var(beta_sim[2,]),var(beta_sim[3,]),
                       mean(betavar_rb[1,]),mean(betavar_rb[2,]),mean(betavar_rb[3,]),
                       mean(betavar_nai[1,]),mean(betavar_nai[2,]),mean(betavar_nai[3,])),p+1,4)
colnames(table_reg_NB1)=c("m()","S^2()","var_rb","var_nai")
row.names(table_reg_NB1)=c("beta0","beta1","beta2")
table_reg_NB1

#robust wald test result
reject_W=matrix(0,3,N)
reject_W=+(walds_beta_rb>qchisq(0.95,1))
mean(reject_W[1,])
mean(reject_W[2,])
mean(reject_W[3,])
#naive wald test result
nareject_W=matrix(0,3,N)
nareject_W=+(walds_beta_nai>qchisq(0.95,1))
mean(nareject_W[1,])
mean(nareject_W[2,])
mean(nareject_W[3,])

#robust score test result
reject_S=matrix(0,3,N)
reject_S=+(score_beta_rb>qchisq(0.95,1))
mean(reject_S[1,])
mean(reject_S[2,])
mean(reject_S[3,])
#naive score test result
nareject_S=matrix(0,3,N)
nareject_S=+(score_beta_nai>qchisq(0.95,1))
mean(nareject_S[1,])
mean(nareject_S[2,])
mean(nareject_S[3,])

#robust LR test result
LR_beta_rb=LR_beta_rb*(+LR_beta_rb>0)
reject_LR=matrix(0,3,N)
reject_LR=+(LR_beta_rb>qchisq(0.95,1))
mean(reject_LR[1,])
mean(reject_LR[2,])
mean(reject_LR[3,])
#naive LR test result
LR_beta_nai=LR_beta_nai*(+LR_beta_nai>0)
nareject_LR=matrix(0,3,N)
nareject_LR=+(LR_beta_nai>qchisq(0.95,1))
mean(nareject_LR[1,])
mean(nareject_LR[2,])
mean(nareject_LR[3,])

#robust & naive wald test result
reject_W_all=matrix(0,2,N)
reject_W_all=+(walds_allpara>qchisq(0.95,3))
mean(reject_W_all[1,])
mean(reject_W_all[2,])

#robust & naive score test result
reject_S_all=matrix(0,2,N)
reject_S_all=+(score_allpara>qchisq(0.95,3))
mean(reject_S_all[1,])
mean(reject_S_all[2,])
#------------------------------------------------------------------------------
