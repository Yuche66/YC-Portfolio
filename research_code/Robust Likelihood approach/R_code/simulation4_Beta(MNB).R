#==============================================================================
#---------------------simulation part 4_Beta (MNB model)-----------------------
#==============================================================================

#package-----------------------------------------------------------------------
library(numDeriv)
library(copula)
library(MASS)
library(pracma)
library(mvtnorm)
library(tictoc)
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
unit_inverse_gaussian_pdf=function(x,u,phi){
    term=phi/(2*(u^2)*log(x))*((log(x)+u)^2)
    ans=sqrt(phi/(2*pi))*(1/(x*(-log(x))^(3/2)))*exp(term)
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
#beta model precision parameters (5,5)
#Normal copula omega=0.6
#------------------------------------------------------------------------------
#setting
N=500
n=100
dm=2
beta_t=c(-1,1,-0.5,1)
p=length(beta_t)/dm-1
YYY_mix=matrix(0,n,dm)
aphb=5
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
XXX=XXXtp

#Normal copula setting
normal_cupcor_mix=0.6

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

#data comes from beta distribution---------------------------------------------
times=0
tic.clearlog()
for(M in c(1:500)){
    
    #time start
    tic(M)
    
    #generating data
    set.seed(829+M*89)
    for(i in c(1:n)){
        mixcopula=normalCopula(normal_cupcor_mix,dim=dm,dispstr = "un")
        mvmix=mvdc(mixcopula,margins=c("beta","beta"),
                   paramMargins = list(list(aphb*mu[dm*i-dm+1],aphb*(1-mu[dm*i-dm+1])),
                                       list(aphug*mu[dm*i-dm+2],aphug*(1-mu[dm*i-dm+2]))))
        temp=c(rMvdc(1,mvmix))
        YYY_mix[i,]=temp
    }
    
    #log likelihood function
    max_f=function(beta){
        -mvnb_loglike_reg_rem_full_id.p(YYY_mix,XXX,beta,link,ph)
    }
    
    #optimize
    v=beta_t
    opt=optim(v,max_f,hessian=T)
    #To ensure it is maximum, repeat optim() until convergence
    if(opt$counts[1]>=501){
        repeat{
            print(c(opt$value,"no converge"))
            opt=optim(opt$par,max_f,hessian=T)
            if(opt$counts[1]<501)break
        }
    }
    print(c(opt$value,"converge"))
    beta_hat=opt$par
    lr=opt$value
    beta_sim[,M]=beta_hat
    
    #I hat, V hat, naive variance hat, robust variance hat
    I=mvnb_fisher_reg_rem_full_id.p(YYY_mix,XXX,beta_hat,link,dlink,ph)
    V=mvnb_scorevar_reg_rem_full_id.p(YYY_mix,XXX,beta_hat,link,dlink,ph)
    var_na=solve(I)/n
    var_rb=solve(I)%*%V%*%solve(I)/n
    
    #three kinds of tests
    for(q in c(1:((p+1)*dm))){
        print(c("test",q))
        betavar_nai[q,M]=var_na[q,q]
        betavar_rb[q,M]=var_rb[q,q]
        A=mvnb_A_reg_rem_full_id.p(YYY_mix,XXX,beta_hat,link,dlink,ph,q)
        B=mvnb_B_reg_rem_full_id.p(YYY_mix,XXX,beta_hat,link,dlink,ph,q)
        
        #LR test 
        max_f0=function(beta0){
            beta=rep(beta_t[q],(p+1)*dm)
            beta[-q]=beta0
            ans=-mvnb_loglike_reg_rem_full_id.p(YYY_mix,XXX,beta,link,ph)
            return(ans)
        }
        opt0=optim(v[-q],max_f0)
        beta_0_hat=opt0$par
        lr0=opt0$value
        LR_beta_nai[q,M]=2*(lr0-lr)
        LR_beta_rb[q,M]=2*(A/B)*(lr0-lr)
        
        #walds test
        w1=A%*%matrix(beta_hat[q]-beta_t[q],1,1)
        walds_beta_rb[q,M]=n*t(w1)%*%solve(B)%*%w1
        walds_beta_nai[q,M]=n*t(w1)%*%solve(A)%*%w1
        
        #score test
        beta=rep(beta_t[q],(p+1)*dm)
        beta[-q]=beta_0_hat
        sA=mvnb_A_reg_rem_full_id.p(YYY_mix,XXX,beta,link,dlink,ph,q)
        sB=mvnb_B_reg_rem_full_id.p(YYY_mix,XXX,beta,link,dlink,ph,q)
        s1=mvnb_score_beta_rem_full_id.p(YYY_mix,XXX,beta,link,dlink,ph,q)
        score_beta_rb[q,M]=1/n*t(s1)%*%solve(sB)%*%s1
        score_beta_nai[q,M]=1/n*t(s1)%*%solve(sA)%*%s1
    }
    times=times+1
    print(times)
    
    #time stop
    toc(log = TRUE, quiet = TRUE)
}
log.lst=tic.log(format=FALSE)
timings=unlist(lapply(log.lst,function(x){x$toc-x$tic}))
#average of taking time
mean(timings)
#------------------------------------------------------------------------------

#data comes from beta dist. result---------------------------------------------
#dimension1
#result(mean, sample variance, variance hat)
table_reg_NB1=matrix(c(mean(beta_sim[1,1:500]),mean(beta_sim[2,1:500]),
                       var(beta_sim[1,1:500]),var(beta_sim[2,1:500]),
                       mean(betavar_rb[1,1:500]),mean(betavar_rb[2,1:500]),
                       mean(betavar_nai[1,1:500]),mean(betavar_nai[2,1:500])),p+1,4)
colnames(table_reg_NB1)=c("m()","S^2()","var_rb","var_nai")
row.names(table_reg_NB1)=c("beta01","beta11")

#dimension2
#result(mean, sample variance, variance hat)
table_reg_NB2=matrix(c(mean(beta_sim[3,1:500]),mean(beta_sim[4,1:500]),
                       var(beta_sim[3,1:500]),var(beta_sim[4,1:500]),
                       mean(betavar_rb[3,1:500]),mean(betavar_rb[4,1:500]),
                       mean(betavar_nai[3,1:500]),mean(betavar_nai[4,1:500])),p+1,4)
colnames(table_reg_NB2)=c("m()","S^2()","var_rb","var_nai")
row.names(table_reg_NB2)=c("beta02","beta12")
table_reg_NB1
table_reg_NB2

#robust wald test result
reject_W=matrix(0,4,N)
reject_W=+(walds_beta_rb>qchisq(0.95,1))
mean(reject_W[1,1:500])
mean(reject_W[2,1:500])
mean(reject_W[3,1:500])
mean(reject_W[4,1:500])
#naive wald test result
nareject_W=matrix(0,4,N)
nareject_W=+(walds_beta_nai>qchisq(0.95,1))
mean(nareject_W[1,1:500])
mean(nareject_W[2,1:500])
mean(nareject_W[3,1:500])
mean(nareject_W[4,1:500])

#robust score test result
reject_S=matrix(0,4,N)
reject_S=+(walds_beta_rb>qchisq(0.95,1))
mean(reject_S[1,1:500])
mean(reject_S[2,1:500])
mean(reject_S[3,1:500])
mean(reject_S[4,1:500])
#naive score test result
nareject_S=matrix(0,4,N)
nareject_S=+(walds_beta_nai>qchisq(0.95,1))
mean(nareject_S[1,1:500])
mean(nareject_S[2,1:500])
mean(nareject_S[3,1:500])
mean(nareject_S[4,1:500])

#robust LR test result
LR_beta_rb=LR_beta_rb*(+LR_beta_rb>0)
reject_LR=matrix(0,4,N)
reject_LR=+(LR_beta_rb>qchisq(0.95,1))
mean(reject_LR[1,1:500])
mean(reject_LR[2,1:500])
mean(reject_LR[3,1:500])
mean(reject_LR[4,1:500])
#naive LR test result
LR_beta_nai=LR_beta_nai*(+LR_beta_nai>0)
nareject_LR=matrix(0,4,N)
nareject_LR=+(LR_beta_nai>qchisq(0.95,1))
mean(nareject_LR[1,1:500])
mean(nareject_LR[2,1:500])
mean(nareject_LR[3,1:500])
mean(nareject_LR[4,1:500])
#------------------------------------------------------------------------------


