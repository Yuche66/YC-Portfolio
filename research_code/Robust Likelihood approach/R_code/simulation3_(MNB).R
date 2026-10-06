#==============================================================================
#---------------------simulation part 3 (MNB model)----------------------------
#==============================================================================

#package-----------------------------------------------------------------------
library(numDeriv)
library(copula)
library(MASS)
library(pracma)
library(mvtnorm)
#------------------------------------------------------------------------------

#MNB log likelihood func., score func., I, V, A, B-----------------------------
mvnb_loglike_reg_rem_full_id.p=function(y_m,x_m,beta,link_f,phi){
    Y=y_m
    m=length(y_m[,1])
    p=length(x_m[1,])
    mn=length(x_m[,1])
    dm=length(y_m[1,])
    X=cbind(matrix(1,mn,1),x_m)
    bt=matrix(c(beta),mn,(p+1),byrow=T)
    itatemp=X*bt
    ita=rep(0,mn)
    for(k in c(1:mn)){
        ita[k]=sum(itatemp[k,])
    }
    mu=sapply(ita,link_f)
    detect= +(mu<=0)
    if(sum(detect)>0){
        loglike=-Inf
    }else{
        loglike=matrix(0,1,m)
        nonaY=c(t(Y))
        temp1=sum(nonaY*log(c(mu)))
        nn=0
        for(i in c(1:m)){
            YY=c(Y[i,])
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
mvnb_score_beta_rem_full_id.p=function(y_m,x_m,beta,link_f,dlink_f,phi,ith){
    Y=y_m
    m=length(y_m[,1])
    p=length(x_m[1,])
    mn=length(x_m[,1])
    dm=length(y_m[1,])
    X=cbind(matrix(1,mn,1),x_m)
    bt=matrix(c(beta),mn,(p+1),byrow=T)
    itatemp=X*bt
    ita=rep(0,mn)
    for(k in c(1:mn)){
        ita[k]=sum(itatemp[k,])
    }
    mu=sapply(ita,link_f)
    gp=sapply(ita,dlink_f)
    nonaY=c(t(Y))
    Xjr=c(t(X))[c(c(1:m)*dm*(p+1)-(dm*(p+1)-ith))]
    nn=0
    score_f=rep(0,m)
    for(i in c(1:m)){
        YY=c(Y[i,])
        Yiplus=sum(YY)
        muiplus=0
        for(j in c((nn+1):(nn+length(YY)))){
            muiplus=muiplus+c(mu)[j]
        }
        score_f[i]=(1+phi*Yiplus)/(1+phi*muiplus)
        nn=nn+length(YY)
    }
    temp1=c(nonaY/c(mu))[c(c(1:m)*dm-(dm-floor((ith-0.1)/(p+1))-1))]
    gp2=gp[c(c(1:m)*dm-(dm-floor((ith-0.1)/(p+1))-1))]
    score_f=(temp1-score_f)*gp2*Xjr
    score_f=sum(score_f)
    return(score_f)
}
mvnb_fisher_reg_rem_full_id.p=function(y_m,x_m,beta,link_f,dlink_f,phi){
    Y=y_m
    m=length(y_m[,1])
    p=length(x_m[1,])
    mn=length(x_m[,1])
    dm=length(y_m[1,])
    X=cbind(matrix(1,mn,1),x_m)
    bt=matrix(c(beta),mn,(p+1),byrow=T)
    itatemp=X*bt
    ita=rep(0,mn)
    for(k in c(1:mn)){
        ita[k]=sum(itatemp[k,])
    }
    mu=sapply(ita,link_f)
    gp=sapply(ita,dlink_f)
    nonaY=c(t(Y))
    I0=matrix(0,mn,dm*(p+1))
    I1=matrix(0,m,dm*(p+1))
    fish_f=matrix(0,m,m)
    mu_m=matrix(0,mn,mn)
    nn=0
    for(i in c(1:m)){
        YY=c(Y[i,])
        muiplus=0
        kk=c(c(1:dm)*m-(m-i))
        for(j in c((nn+1):(nn+length(YY)))){
            muiplus=muiplus+c(mu)[j]
            diag(mu_m)[kk[c(j-nn)]]=1/c(mu)[j]
        }
        fish_f[i,i]=phi/(1+phi*muiplus)
        nn=nn+length(YY)
    }
    for(ith in  c(1:(dm*(p+1)))){
        Xjr=c(t(X))[c(c(1:m)*dm*(p+1)-(dm*(p+1)-ith))]
        gp2=gp[c(c(1:m)*dm-(dm-floor((ith-0.1)/(p+1))-1))]
        I1[,ith]=gp2*Xjr
        I0[c(c(1:m)+(m*floor((ith-0.1)/(p+1)))),ith]=gp2*Xjr
    }
    Fish_m=(t(I0)%*%mu_m%*%I0-t(I1)%*%fish_f%*%I1)/m
    return(Fish_m)
}
mvnb_scorevar_reg_rem_full_id.p=function(y_m,x_m,beta,link_f,dlink_f,phi){
    Y=y_m
    m=length(y_m[,1])
    p=length(x_m[1,])
    mn=length(x_m[,1])
    dm=length(y_m[1,])
    X=cbind(matrix(1,mn,1),x_m)
    bt=matrix(c(beta),mn,(p+1),byrow=T)
    itatemp=X*bt
    ita=rep(0,mn)
    for(k in c(1:mn)){
        ita[k]=sum(itatemp[k,])
    }
    mu=sapply(ita,link_f)
    gp=sapply(ita,dlink_f)
    nonaY=c(t(Y))
    V0=matrix(0,m,dm*(p+1))
    score_f=rep(0,m)
    nn=0
    for(i in c(1:m)){
        YY=c(Y[i,])
        Yiplus=sum(YY)
        muiplus=0
        for(j in c((nn+1):(nn+length(YY)))){
            muiplus=muiplus+c(mu)[j]
        }
        score_f[i]=(1+phi*Yiplus)/(1+phi*muiplus)
        nn=nn+length(YY)
    }
    for(ith in  c(1:(dm*(p+1)))){
        Xjr=c(t(X))[c(c(1:m)*dm*(p+1)-(dm*(p+1)-ith))]
        temp1=c(nonaY/c(mu))[c(c(1:m)*dm-(dm-floor((ith-0.1)/(p+1))-1))]
        gp2=gp[c(c(1:m)*dm-(dm-floor((ith-0.1)/(p+1))-1))]
        ans=(temp1-score_f)*gp2*Xjr
        V0[,ith]=ans
    }
    scoreV=var(V0)
    return(scoreV)
}
mvnb_A_reg_rem_full_id.p=function(y_m,x_m,beta,link_f,dlink_f,phi,ith){
    Fish_m=mvnb_fisher_reg_rem_full_id.p(y_m,x_m,beta,link_f,dlink_f,phi)
    I1=Fish_m[ith,ith]
    I2=Fish_m[-ith,ith]
    I3=Fish_m[-ith,-ith]
    A=I1-t(I2)%*%solve(I3)%*%I2
    return(A)
}
mvnb_B_reg_rem_full_id.p=function(y_m,x_m,beta,link_f,dlink_f,phi,ith){
    Fish_m=mvnb_fisher_reg_rem_full_id.p(y_m,x_m,beta,link_f,dlink_f,phi)
    scoreV=mvnb_scorevar_reg_rem_full_id.p(y_m,x_m,beta,link_f,dlink_f,phi)
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
link=function(x){
    1-1/(1+exp(x))
}
dlink=function(x){
    exp(x)/((1+exp(x))^2)
}
#------------------------------------------------------------------------------


#start-------------------------------------------------------------------------
#simulation times N=3000, sample size n=100, 200
#beta_t=(0.5,1,-1.5,-0.5,1.5,1) & XXX1~N(0,0.3^2), XXX2~Ber(0.5)
#beta model aphb=(9,4), unit gamma model aphug=(3,5)
#random effect variance matrix sigma=0.1, 0.5, rho=0.5
#------------------------------------------------------------------------------
#setting
N=3000
n=100
dm=2
beta_t=c(0.5,1,-1.5,-0.5,1.5,1)
#the approximation of regression (Zeger)
c=16*sqrt(3)/(15*pi)
sg=0.5
beta_adj=(((c^2)*(sg^2)+1)^(-1/2))*beta_t
p=length(beta_t)/dm-1
YYY_bt=matrix(0,n,dm)
YYY_ug=matrix(0,n,dm)
aphb=c(9,4)
aphug=c(3,5)

#random effect setting
rdefflo=c(0.5)
Crm=matrix(1,dm,dm)
Crm[upper.tri(Crm)]=rdefflo
Crm[lower.tri(Crm)]=rdefflo
G=diag(rep(sg,dm),dm,dm)%*%Crm%*%diag(rep(sg,dm),dm,dm)
L=t(chol(G))

#mu setting
set.seed(923)
XXX1=matrix(rnorm(n*dm,0,0.3),n*dm,1)
XXX2=matrix(rbinom(n*dm,1,0.5),n*dm,1)
XXX=cbind(XXX1,XXX2)
it1=XXX[1:n,]%*%matrix(beta_t[c(2,3)],p,1)+beta_t[1]
it2=XXX[(n+1):(n*dm),]%*%matrix(beta_t[c(5,6)],p,1)+beta_t[4]
it=rep(0,dm*n)
it[c(c(1:n)*dm-(dm-1))]=it1
it[c(c(1:n)*dm-(dm-2))]=it2
mufix=sapply(it,link)
matrix(mufix,n,dm,byrow=T)
XXXtp=matrix(0,n*dm,p)
XXXtp[c(c(1:n)*dm-(dm-1)),]=XXX[1:n,]
XXXtp[c(c(1:n)*dm-(dm-2)),]=XXX[(n+1):(n*dm),]
XXX=XXXtp

#scatter
gd=c(beta_t)

#nuisance parameter phi=1
h=0.02*c(1:5000)
num=10
ph=1

#record the result
beta_sim=matrix(0,(p+1)*dm,N)
betavar_nai=matrix(0,(p+1)*dm,N)
betavar_rb=matrix(0,(p+1)*dm,N)
walds_beta_rb=matrix(0,(p+1)*dm,N)
walds_beta_nai=matrix(0,(p+1)*dm,N)
score_beta_rb=matrix(0,(p+1)*dm,N)
score_beta_nai=matrix(0,(p+1)*dm,N)
LR_beta_rb=matrix(0,(p+1)*dm,N)
LR_beta_nai=matrix(0,(p+1)*dm,N)

#data comes from beta----------------------------------------------------------
times=0
for(M in c(1:N)){
    
    #generating data
    set.seed(619+M*210)
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
    
    #log likelihood function
    max_f=function(beta){
        -mvnb_loglike_reg_rem_full_id.p(YYY_bt,XXX,beta,link,ph)
    }
    
    #initial point
    v=rep(0,(p+1)*dm)
    for(j in c(1:(p+1))){
        v[j]=runif(1,beta_t[j]-0.5,beta_t[j]+0.5)
    }
    
    #optimize
    opt=optim(v,max_f)
    beta_hat=opt$par
    lr=opt$value
    beta_sim[,M]=beta_hat
    
    #I hat, V hat, naive variance hat, robust variance hat
    I=mvnb_fisher_reg_rem_full_id.p(YYY_bt,XXX,beta_hat,link,dlink,ph)
    V=mvnb_scorevar_reg_rem_full_id.p(YYY_bt,XXX,beta_hat,link,dlink,ph)
    var_na=solve(I)/n
    var_rb=solve(I)%*%V%*%solve(I)/n
    
    #three kinds of tests
    for(q in c(1:((p+1)*dm))){
        betavar_nai[q,M]=var_na[q,q]
        betavar_rb[q,M]=var_rb[q,q]
        A=mvnb_A_reg_rem_full_id.p(YYY_bt,XXX,beta_hat,link,dlink,ph,q)
        B=mvnb_B_reg_rem_full_id.p(YYY_bt,XXX,beta_hat,link,dlink,ph,q)
        
        #LR test 
        max_f0=function(beta0){
            beta=rep(beta_adj[q],(p+1)*dm)
            beta[-q]=beta0
            ans=-mvnb_loglike_reg_rem_full_id.p(YYY_bt,XXX,beta,link,ph)
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
        beta=rep(beta_adj[q],(p+1)*dm)
        beta[-q]=beta_0_hat
        sA=mvnb_A_reg_rem_full_id.p(YYY_bt,XXX,beta,link,dlink,ph,q)
        sB=mvnb_B_reg_rem_full_id.p(YYY_bt,XXX,beta,link,dlink,ph,q)
        s1=mvnb_score_beta_rem_full_id.p(YYY_bt,XXX,beta,link,dlink,ph,q)
        score_beta_rb[q,M]=1/n*t(s1)%*%solve(sB)%*%s1
        score_beta_nai[q,M]=1/n*t(s1)%*%solve(sA)%*%s1
    }
    
    times=times+1
    print(times)
}
#------------------------------------------------------------------------------

#data comes from beta result---------------------------------------------------
#dimension1
#result(mean, sample variance, variance hat)
table_reg_NB1=matrix(c(mean(beta_sim[1,]),mean(beta_sim[2,]),mean(beta_sim[3,]),
                       var(beta_sim[1,]),var(beta_sim[2,]),var(beta_sim[3,]),
                       mean(betavar_rb[1,]),mean(betavar_rb[2,]),mean(betavar_rb[3,]),
                       mean(betavar_nai[1,]),mean(betavar_nai[2,]),mean(betavar_nai[3,])),p+1,4)
colnames(table_reg_NB1)=c("m()","S^2()","var_rb","var_nai")
row.names(table_reg_NB1)=c("beta01","beta11","beta21")
table_reg_NB1

#dimension2
#result(mean, sample variance, variance hat)
table_reg_NB2=matrix(c(mean(beta_sim[4,]),mean(beta_sim[5,]),mean(beta_sim[6,]),
                       var(beta_sim[4,]),var(beta_sim[5,]),var(beta_sim[6,]),
                       mean(betavar_rb[4,]),mean(betavar_rb[5,]),mean(betavar_rb[6,]),
                       mean(betavar_nai[4,]),mean(betavar_nai[5,]),mean(betavar_nai[6,])),p+1,4)
colnames(table_reg_NB2)=c("m()","S^2()","var_rb","var_nai")
row.names(table_reg_NB2)=c("beta02","beta12","beta22")
table_reg_NB2

#dimension1
#robust wald test result
reject_W=matrix(0,6,N)
reject_W=+(walds_beta_rb>qchisq(0.95,1))
mean(reject_W[1,])
mean(reject_W[2,])
mean(reject_W[3,])
#naive wald test result
nareject_W=matrix(0,6,N)
nareject_W=+(walds_beta_nai>qchisq(0.95,1))
mean(nareject_W[1,])
mean(nareject_W[2,])
mean(nareject_W[3,])

#robust score test result
reject_S=matrix(0,6,N)
reject_S=+(score_beta_rb>qchisq(0.95,1))
mean(reject_S[1,])
mean(reject_S[2,])
mean(reject_S[3,])
#naive score test result
nareject_S=matrix(0,6,N)
nareject_S=+(score_beta_nai>qchisq(0.95,1))
mean(nareject_S[1,])
mean(nareject_S[2,])
mean(nareject_S[3,])

#robust LR test result
LR_beta_rb=LR_beta_rb*(+LR_beta_rb>0)
reject_LR=matrix(0,6,N)
reject_LR=+(LR_beta_rb>qchisq(0.95,1))
mean(reject_LR[1,])
mean(reject_LR[2,])
mean(reject_LR[3,])
#naive LR test result
LR_beta_nai=LR_beta_nai*(+LR_beta_nai>0)
nareject_LR=matrix(0,6,N)
nareject_LR=+(LR_beta_nai>qchisq(0.95,1))
mean(nareject_LR[1,])
mean(nareject_LR[2,])
mean(nareject_LR[3,])

#dimension2
#robust wald test result
mean(reject_W[4,])
mean(reject_W[5,])
mean(reject_W[6,])
#naive wald test result
mean(nareject_W[4,])
mean(nareject_W[5,])
mean(nareject_W[6,])

#robust score test result
mean(reject_S[4,])
mean(reject_S[5,])
mean(reject_S[6,])
#naive score test result
mean(nareject_S[4,])
mean(nareject_S[5,])
mean(nareject_S[6,])

#robust LR test result
mean(reject_LR[4,])
mean(reject_LR[5,])
mean(reject_LR[6,])
#naive LR test result
mean(nareject_LR[4,])
mean(nareject_LR[5,])
mean(nareject_LR[6,])
#------------------------------------------------------------------------------


#data comes from unit gamma----------------------------------------------------
times=0
for(M in c(1:N)){
    
    #generating data
    set.seed(619*M+210)
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
    
    #loglikelihood function
    max_f=function(beta){
        -mvnb_loglike_reg_rem_full_id.p(YYY_ug,XXX,beta,link,ph)
    }
    
    #initial point
    v=rep(0,(p+1)*dm)
    for(j in c(1:((p+1)*dm))){
        v[j]=runif(1,beta_t[j]-0.5,beta_t[j]+0.5)
    }
    
    #optimize
    opt=optim(v,max_f)
    beta_hat=opt$par
    lr=opt$value
    beta_sim[,M]=beta_hat
    
    #I hat, V hat, naive variance hat, robust variance hat
    I=mvnb_fisher_reg_rem_full_id.p(YYY_ug,XXX,beta_hat,link,dlink,ph)
    V=mvnb_scorevar_reg_rem_full_id.p(YYY_ug,XXX,beta_hat,link,dlink,ph)
    var_na=solve(I)/n
    var_rb=solve(I)%*%V%*%solve(I)/n
    
    #three kinds of tests
    for(q in c(1:((p+1)*dm))){
        betavar_nai[q,M]=var_na[q,q]
        betavar_rb[q,M]=var_rb[q,q]
        A=mvnb_A_reg_rem_full_id.p(YYY_ug,XXX,beta_hat,link,dlink,ph,q)
        B=mvnb_B_reg_rem_full_id.p(YYY_ug,XXX,beta_hat,link,dlink,ph,q)
        
        #LR test
        max_f0=function(beta0){
            beta=rep(beta_adj[q],(p+1)*dm)
            beta[-q]=beta0
            ans=-mvnb_loglike_reg_rem_full_id.p(YYY_ug,XXX,beta,link,ph)
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
        beta=rep(beta_adj[q],(p+1)*dm)
        beta[-q]=beta_0_hat
        sA=mvnb_A_reg_rem_full_id.p(YYY_ug,XXX,beta,link,dlink,ph,q)
        sB=mvnb_B_reg_rem_full_id.p(YYY_ug,XXX,beta,link,dlink,ph,q)
        s1=mvnb_score_beta_rem_full_id.p(YYY_ug,XXX,beta,link,dlink,ph,q)
        score_beta_rb[q,M]=1/n*t(s1)%*%solve(sB)%*%s1
        score_beta_nai[q,M]=1/n*t(s1)%*%solve(sA)%*%s1
    }
    
    times=times+1
    print(times)
}
#------------------------------------------------------------------------------

#data comes from unit gamma result---------------------------------------------
#dimension1
#result(mean, sample variance, variance hat)
table_reg_NB1=matrix(c(mean(beta_sim[1,]),mean(beta_sim[2,]),mean(beta_sim[3,]),
                       var(beta_sim[1,]),var(beta_sim[2,]),var(beta_sim[3,]),
                       mean(betavar_rb[1,]),mean(betavar_rb[2,]),mean(betavar_rb[3,]),
                       mean(betavar_nai[1,]),mean(betavar_nai[2,]),mean(betavar_nai[3,])),p+1,4)
colnames(table_reg_NB1)=c("m()","S^2()","var_rb","var_nai")
row.names(table_reg_NB1)=c("beta01","beta11","beta21")
table_reg_NB1

#dimension2
#result(mean, sample variance, variance hat)
table_reg_NB2=matrix(c(mean(beta_sim[4,]),mean(beta_sim[5,]),mean(beta_sim[6,]),
                       var(beta_sim[4,]),var(beta_sim[5,]),var(beta_sim[6,]),
                       mean(betavar_rb[4,]),mean(betavar_rb[5,]),mean(betavar_rb[6,]),
                       mean(betavar_nai[4,]),mean(betavar_nai[5,]),mean(betavar_nai[6,])),p+1,4)
colnames(table_reg_NB2)=c("m()","S^2()","var_rb","var_nai")
row.names(table_reg_NB2)=c("beta02","beta12","beta22")
table_reg_NB2

#dimension1
#robust wald test result
reject_W=matrix(0,6,N)
reject_W=+(walds_beta_rb>qchisq(0.95,1))
mean(reject_W[1,])
mean(reject_W[2,])
mean(reject_W[3,])
#naive wald test result
nareject_W=matrix(0,6,N)
nareject_W=+(walds_beta_nai>qchisq(0.95,1))
mean(nareject_W[1,])
mean(nareject_W[2,])
mean(nareject_W[3,])

#robust score test result
reject_S=matrix(0,6,N)
reject_S=+(score_beta_rb>qchisq(0.95,1))
mean(reject_S[1,])
mean(reject_S[2,])
mean(reject_S[3,])
#naive score test result
nareject_S=matrix(0,6,N)
nareject_S=+(score_beta_nai>qchisq(0.95,1))
mean(nareject_S[1,])
mean(nareject_S[2,])
mean(nareject_S[3,])

#robust LR test result
LR_beta_rb=LR_beta_rb*(+LR_beta_rb>0)
reject_LR=matrix(0,6,N)
reject_LR=+(LR_beta_rb>qchisq(0.95,1))
mean(reject_LR[1,])
mean(reject_LR[2,])
mean(reject_LR[3,])
#naive LR test result
LR_beta_nai=LR_beta_nai*(+LR_beta_nai>0)
nareject_LR=matrix(0,6,N)
nareject_LR=+(LR_beta_nai>qchisq(0.95,1))
mean(nareject_LR[1,])
mean(nareject_LR[2,])
mean(nareject_LR[3,])

#dimension2
#robust wald test result
mean(reject_W[4,])
mean(reject_W[5,])
mean(reject_W[6,])
#naive wald test result
mean(nareject_W[4,])
mean(nareject_W[5,])
mean(nareject_W[6,])

#robust score test result
mean(reject_S[4,])
mean(reject_S[5,])
mean(reject_S[6,])
#naive score test result
mean(nareject_S[4,])
mean(nareject_S[5,])
mean(nareject_S[6,])

#robust LR test result
mean(reject_LR[4,])
mean(reject_LR[5,])
mean(reject_LR[6,])
#naive LR test result
mean(nareject_LR[4,])
mean(nareject_LR[5,])
mean(nareject_LR[6,])
#------------------------------------------------------------------------------

