#==============================================================================
#---------------------simulation part 5-3 (MNB model)--------------------------
#==============================================================================

#package-----------------------------------------------------------------------
library(numDeriv)
library(copula)
library(MASS)
library(pracma)
library(mvtnorm)
library(statmod)
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
#beta, unit gamma & unit inverse gaussian pdf, cdf
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
#simulation times N=3000, sample size n=50,100,200,400
#beta_t=(-1,0.5,0,1,0,0.5,0.5,0.5,-0.5) & XXX1~N(0,1), XXX2~Ber(0.5),
#beta model aphb=9, unit gamma model aphug=9
#Normal copula omega12=0, omega13=0.5, omega23=0.8
#------------------------------------------------------------------------------
#setting
N=3000
n=400
dm=3
beta_t=c(-1,0.5,0,1,0,0.5,0.5,0.5,-0.5)
p=length(beta_t)/dm-1
aphb=9
aphug=9
YYY_bui=matrix(0,n,dm)

#mu setting
set.seed(210)
XXX1=matrix(rnorm(dm*n,0,1),dm*n,1)
XXX2=matrix(rbinom(dm*n,1,0.5),dm*n,1)
XXX12=cbind(XXX1,XXX2)
XXXtp=matrix(0,n*dm,p)
XXXtp[c(c(1:n)*dm-(dm-1)),]=XXX12[1:n,]
XXXtp[c(c(1:n)*dm-(dm-2)),]=XXX12[c((n+1):(2*n)),]
XXXtp[c(c(1:n)*dm-(dm-3)),]=XXX12[c((2*n+1):(3*n)),]
XXX=XXXtp
it1=XXX12[1:n,]%*%matrix(beta_t[2:3],p,1)+beta_t[1]
it2=XXX12[c((n+1):(2*n)),]%*%matrix(beta_t[5:6],p,1)+beta_t[4]
it3=XXX12[c((2*n+1):(3*n)),]%*%matrix(beta_t[8:9],p,1)+beta_t[7]
it=rep(0,dm*n)
it[c(c(1:n)*dm-(dm-1))]=it1
it[c(c(1:n)*dm-(dm-2))]=it2
it[c(c(1:n)*dm-(dm-3))]=it3
mufix=sapply(it,link)
mu=mufix
matrix(mufix,n,dm,byrow=T)

#Normal copula setting
normal_cupcor_bui=c(0,0.5,0.8)

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

#data comes from different distributions each dimension------------------------
times=0
for(M in c(1:N)){
    
    #generating data
    set.seed(619+21*M)
    for(i in c(1:n)){
        #generate correlated data from beta, gamma & inverse gaussian dist.
        buicopula=normalCopula(normal_cupcor_bui,dim=dm,dispstr = "un")
        mvbui=mvdc(buicopula,margins=c("beta","gamma","invgauss"),
                   paramMargins = list(list(aphb*mu[dm*i-dm+1],aphb*(1-mu[dm*i-dm+1])),
                                       list(shape=aphug,scale=(((mu[dm*i-dm+2])^(-1/aphug))-1)),
                                       list(0.5/mu[dm*i-dm+3],((log(mu[dm*i-dm+3]))^2)/(2*(1+2*mu[dm*i-dm+3]*log(mu[dm*i-dm+3]))))))
        temp=c(rMvdc(1,mvbui))
        #exp(-X) transformation for gamma & inverse gaussian dist.
        YYY_bui[i,]=c(temp[1],exp(-temp[2:3]))
    }
    
    #log likelihood function
    max_f=function(beta){
        -mvnb_loglike_reg_rem_full_id.p(YYY_bui,XXX,beta,link,ph)
    }
    
    #optimize
    v=beta_t
    opt=optim(v,max_f,hessian=T)
    loglike=opt$value
    #To ensure it is maximum, repeat optim() until convergence
    repeat{
        print(c(opt$value,"no max"))
        opt=optim(opt$par,max_f,hessian=T)
        if((opt$value-loglike)<=0){
            a=abs(opt$value-loglike)
        }
        #convergence tolerance
        if(opt$counts[1]<501 & a<0.000001)break
        if(a>=0.000001){
            loglike=opt$value
        }
    }
    print(c(opt$value,"max"))
    beta_hat=opt$par
    lr=opt$value
    beta_sim[,M]=beta_hat
    
    #I hat, V hat, naive variance hat, robust variance hat
    I=opt$hessian/n
    V=mvnb_scorevar_reg_rem_full_id.p(YYY_bui,XXX,beta_hat,link,dlink,ph)
    var_na=solve(I)/n
    var_rb=solve(I)%*%V%*%solve(I)/n
    
    #three kinds of tests
    for(q in c(1:((p+1)*dm))){
        print(c("test",q))
        betavar_nai[q,M]=var_na[q,q]
        betavar_rb[q,M]=var_rb[q,q]
        A=(var_na[q,q]^-1)/n
        B=mvnb_B_reg_rem_full_id.p(YYY_bui,XXX,beta_hat,link,dlink,ph,q)
        
        #LR test
        max_f0=function(beta0){
            beta=rep(beta_t[q],(p+1)*dm)
            beta[-q]=beta0
            ans=-mvnb_loglike_reg_rem_full_id.p(YYY_bui,XXX,beta,link,ph)
            return(ans)
        }
        #given beta_rj and calculate other estimates of parameters
        opt0=optim(v[-q],max_f0)
        loglike0=opt0$value
        #To ensure it is maximum, repeat optim() until convergence
        repeat{
            opt0=optim(opt0$par,max_f0,hessian=T)
            if((opt0$value-loglike0)<=0){
                a0=abs(opt0$value-loglike0)
            }
            #convergence tolerance
            if(opt0$counts[1]<501 & a0<0.000001)break
            if(a0>=0.000001){
                loglike0=opt0$value
            }
            print(c("telda,no",a0))
        }
        print(c("telda,yes"))
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
        sA=mvnb_A_reg_rem_full_id.p(YYY_bui,XXX,beta,link,dlink,ph,q)
        sB=mvnb_B_reg_rem_full_id.p(YYY_bui,XXX,beta,link,dlink,ph,q)
        s1=mvnb_score_beta_rem_full_id.p(YYY_bui,XXX,beta,link,dlink,ph,q)
        score_beta_rb[q,M]=1/n*t(s1)%*%solve(sB)%*%s1
        score_beta_nai[q,M]=1/n*t(s1)%*%solve(sA)%*%s1
    }
    
    times=times+1
    print(times)
}
#------------------------------------------------------------------------------

#data comes from different dist. result----------------------------------------
#dimension1
#result(mean, sample variance, variance hat)
table_reg_NB1=matrix(c(mean(beta_sim[1,]),mean(beta_sim[2,]),mean(beta_sim[3,]),
                       var(beta_sim[1,]),var(beta_sim[2,]),var(beta_sim[3,]),
                       mean(betavar_rb[1,]),mean(betavar_rb[2,]),mean(betavar_rb[3,]),
                       mean(betavar_nai[1,]),mean(betavar_nai[2,]),mean(betavar_nai[3,])),p+1,4)
colnames(table_reg_NB1)=c("m()","S^2()","var_rb","var_nai")
row.names(table_reg_NB1)=c("beta01","beta11","beta21")

#dimension2
#result(mean, sample variance, variance hat)
table_reg_NB2=matrix(c(mean(beta_sim[4,]),mean(beta_sim[5,]),mean(beta_sim[6,]),
                       var(beta_sim[4,]),var(beta_sim[5,]),var(beta_sim[6,]),
                       mean(betavar_rb[4,]),mean(betavar_rb[5,]),mean(betavar_rb[6,]),
                       mean(betavar_nai[4,]),mean(betavar_nai[5,]),mean(betavar_nai[6,])),p+1,4)
colnames(table_reg_NB2)=c("m()","S^2()","var_rb","var_nai")
row.names(table_reg_NB2)=c("beta02","beta12","beta22")

#dimension3
#result(mean, sample variance, variance hat)
table_reg_NB3=matrix(c(mean(beta_sim[7,]),mean(beta_sim[8,]),mean(beta_sim[9,]),
                       var(beta_sim[7,]),var(beta_sim[8,]),var(beta_sim[9,]),
                       mean(betavar_rb[7,]),mean(betavar_rb[8,]),mean(betavar_rb[9,]),
                       mean(betavar_nai[7,]),mean(betavar_nai[8,]),mean(betavar_nai[9,])),p+1,4)
colnames(table_reg_NB3)=c("m()","S^2()","var_rb","var_nai")
row.names(table_reg_NB3)=c("beta03","beta13","beta23")
table_reg_NB1
table_reg_NB2
table_reg_NB3

#robust wald test result
reject_W=matrix(0,9,N)
reject_W=+(walds_beta_rb>qchisq(0.95,1))
mean(reject_W[1,])
mean(reject_W[2,])
mean(reject_W[3,])
mean(reject_W[4,])
mean(reject_W[5,])
mean(reject_W[6,])
mean(reject_W[7,])
mean(reject_W[8,])
mean(reject_W[9,])
#naive wald test result
nareject_W=matrix(0,9,N)
nareject_W=+(walds_beta_nai>qchisq(0.95,1))
mean(nareject_W[1,])
mean(nareject_W[2,])
mean(nareject_W[3,])
mean(nareject_W[4,])
mean(nareject_W[5,])
mean(nareject_W[6,])
mean(nareject_W[7,])
mean(nareject_W[8,])
mean(nareject_W[9,])

#robust score test result
reject_S=matrix(0,9,N)
reject_S=+(score_beta_rb>qchisq(0.95,1))
mean(reject_S[1,])
mean(reject_S[2,])
mean(reject_S[3,])
mean(reject_S[4,])
mean(reject_S[5,])
mean(reject_S[6,])
mean(reject_S[7,])
mean(reject_S[8,])
mean(reject_S[9,])
#naive score test result
nareject_S=matrix(0,9,N)
nareject_S=+(score_beta_nai>qchisq(0.95,1))
mean(nareject_S[1,])
mean(nareject_S[2,])
mean(nareject_S[3,])
mean(nareject_S[4,])
mean(nareject_S[5,])
mean(nareject_S[6,])
mean(nareject_S[7,])
mean(nareject_S[8,])
mean(nareject_S[9,])

#robust LR test result
LR_beta_rb=LR_beta_rb*(+LR_beta_rb>0)
reject_LR=matrix(0,9,N)
reject_LR=+(LR_beta_rb>qchisq(0.95,1))
mean(reject_LR[1,])
mean(reject_LR[2,])
mean(reject_LR[3,])
mean(reject_LR[4,])
mean(reject_LR[5,])
mean(reject_LR[6,])
mean(reject_LR[7,])
mean(reject_LR[8,])
mean(reject_LR[9,])
#naive LR test result
LR_beta_nai=LR_beta_nai*(+LR_beta_nai>0)
nareject_LR=matrix(0,9,N)
nareject_LR=+(LR_beta_nai>qchisq(0.95,1))
mean(nareject_LR[1,])
mean(nareject_LR[2,])
mean(nareject_LR[3,])
mean(nareject_LR[4,])
mean(nareject_LR[5,])
mean(nareject_LR[6,])
mean(nareject_LR[7,])
mean(nareject_LR[8,])
mean(nareject_LR[9,])
#------------------------------------------------------------------------------


#Let sample size of data = 50000 and calculates the correlation roughly--------
N=1
n=50000
dm=3
beta_t=c(-1,0.5,0,1,0,0.5,0.5,0.5,-0.5)
p=length(beta_t)/dm-1
aphb=9
aphug=9
YYY_bui=matrix(0,n,dm)

#mu setting
set.seed(210)
XXX1=matrix(rnorm(dm*n,0,1),dm*n,1)
XXX2=matrix(rbinom(dm*n,1,0.5),dm*n,1)
XXX12=cbind(XXX1,XXX2)
XXXtp=matrix(0,n*dm,p)
XXXtp[c(c(1:n)*dm-(dm-1)),]=XXX12[1:n,]
XXXtp[c(c(1:n)*dm-(dm-2)),]=XXX12[c((n+1):(2*n)),]
XXXtp[c(c(1:n)*dm-(dm-3)),]=XXX12[c((2*n+1):(3*n)),]
XXX=XXXtp
it1=XXX12[1:n,]%*%matrix(beta_t[2:3],p,1)+beta_t[1]
it2=XXX12[c((n+1):(2*n)),]%*%matrix(beta_t[5:6],p,1)+beta_t[4]
it3=XXX12[c((2*n+1):(3*n)),]%*%matrix(beta_t[8:9],p,1)+beta_t[7]
it=rep(0,dm*n)
it[c(c(1:n)*dm-(dm-1))]=it1
it[c(c(1:n)*dm-(dm-2))]=it2
it[c(c(1:n)*dm-(dm-3))]=it3
mufix=sapply(it,link)
mu=mufix
matrix(mufix,n,dm,byrow=T)

#Normal copula setting
normal_cupcor_bui=c(0,0.5,0.8)


set.seed(61921)
for(i in c(1:n)){
    
    #generating data
    buicopula=normalCopula(normal_cupcor_bui,dim=dm,dispstr = "un")
    mvbui=mvdc(buicopula,margins=c("beta","gamma","invgauss"),
               paramMargins = list(list(aphb*mu[dm*i-dm+1],aphb*(1-mu[dm*i-dm+1])),
                                   list(shape=aphug,scale=(((mu[dm*i-dm+2])^(-1/aphug))-1)),
                                   list(0.5/mu[dm*i-dm+3],((log(mu[dm*i-dm+3]))^2)/(2*(1+2*mu[dm*i-dm+3]*log(mu[dm*i-dm+3]))))))
    temp=c(rMvdc(1,mvbui))
    YYY_bui[i,]=c(temp[1],exp(-temp[2:3]))
}
#calculate the correlation of (Yi-mu)
cor(YYY_bui-matrix(mufix,n,dm,byrow=T))

#------------------------------------------------------------------------------


