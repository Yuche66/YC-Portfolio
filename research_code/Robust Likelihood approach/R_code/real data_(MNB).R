#==============================================================================
#-------------real data : body fat percentage data (MNB model)-----------------
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


#MNB log likelihood func., score func., I, V, A, B-----------------------------
#univariate model function
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

#5d model function
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

#link function
link=function(x){
    1-1/(1+exp(x))
}
dlink=function(x){
    exp(x)/((1+exp(x))^2)
}
#------------------------------------------------------------------------------


#----------------------------univariate model MNB------------------------------
#(take the estimates of univariate model as initial point)
#setting(nuisance parameter phi=1)
dm=5
p=5
n=length(bfpY[,1])
df=30
ph=1
YYY=as.matrix(bfpY)
XXX=as.matrix(bfpX)

#====================================arms======================================
d=1

#log likelihood function
max_f=function(beta){
    -mvnb_loglike_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,ph)
}

#optimize
reg_nb1=optim(c(0,0,0,0,0,0),max_f,method="Nelder-Mead",hessian = T)
updatereg=0
if(reg_nb1$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        temppar=reg_nb1$par
        reg_nb1=optim(temppar,max_f,method="Nelder-Mead",hessian = T)
        if(reg_nb1$counts[1]<501 | updatereg==100)break
    }
}
reg_nb1par=reg_nb1$par
lr=reg_nb1$value

#I hat, V hat, naive variance hat, robust variance hat
I=mvnb_fisher_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb1par,link,dlink,ph)
V=mvnb_scorevar_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb1par,link,dlink,ph)
reg_nb1var_nam=solve(I)/n
reg_nb1var_rbm=solve(I)%*%V%*%solve(I)/n

#three kinds of tests
LR_reg_nb1=c()
reg_nb1llh=c()
walds_reg_nb1=c()
score_reg_nb1=c()
for(q in c(1:(p+1))){
    A=mvnb_A_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb1par,link,dlink,ph,q)
    B=mvnb_B_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb1par,link,dlink,ph,q)
    
    #LR test 
    max_f0=function(beta0){
        beta=rep(0,(p+1))
        beta[-q]=beta0
        ans=-mvnb_loglike_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,ph)
        return(ans)
    }
    opt0=optim(c(0,0,0,0,0),max_f0,method="Nelder-Mead",hessian = T)
    updatereg=0
    if(opt0$counts[1]>=501){
        repeat{
            updatereg=updatereg+1
            temppar=opt0$par
            opt0=optim(temppar,max_f0,method="Nelder-Mead",hessian = T)
            if(opt0$counts[1]<501 | updatereg==100)break
        }
    }
    beta_0_hat=opt0$par
    lr0=opt0$value
    LR_reg_nb1[q]=2*(A/B)*(lr0-lr)
    reg_nb1llh[q]=(A/B)*(-lr)
    
    #walds test
    w1=A%*%matrix(reg_nb1par[q]-0,1,1)
    walds_reg_nb1[q]=n*t(w1)%*%solve(B)%*%w1
    
    #score test
    beta=rep(0,(p+1))
    beta[-q]=beta_0_hat
    sA=mvnb_A_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    sB=mvnb_B_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    s1=mvnb_score_beta_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    score_reg_nb1[q]=1/n*t(s1)%*%solve(sB)%*%s1
}
reg_nb1var_na=diag(reg_nb1var_nam)
reg_nb1var_rb=diag(reg_nb1var_rbm)
reg_nb1p_w=1-pchisq(walds_reg_nb1,1)
reg_nb1p_s=1-pchisq(score_reg_nb1,1)
reg_nb1p_lr=1-pchisq(LR_reg_nb1,1)
CInb1L=c(reg_nb1par-qnorm(0.975)*sqrt(reg_nb1var_rb))
CInb1U=c(reg_nb1par+qnorm(0.975)*sqrt(reg_nb1var_rb))

#result
#robust 1,LR test, 2,score test 1,Wald test, p-value 
reg_nb1p_lr
reg_nb1p_s
reg_nb1p_w
#estimates of parameters
reg_nb1par
#estimates of standard deviation(robust)
as.numeric(sqrt(reg_nb1var_rb))
as.numeric(sqrt(reg_nb1var_na))
#Confident interval of Wald test(robust)
as.numeric(CInb1L)
as.numeric(CInb1U)
#==============================================================================

#====================================legs======================================
d=2

#log likelihood function
max_f=function(beta){
    -mvnb_loglike_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,ph)
}

#optimize
reg_nb2=optim(c(0,0,0,0,0,0),max_f,method="Nelder-Mead",hessian = T)
updatereg=0
if(reg_nb2$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        temppar=reg_nb2$par
        reg_nb2=optim(temppar,max_f,method="Nelder-Mead",hessian = T)
        if(reg_nb2$counts[1]<501 | updatereg==100)break
    }
}
reg_nb2par=reg_nb2$par
lr=reg_nb2$value

#I hat, V hat, naive variance hat, robust variance hat
I=mvnb_fisher_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb2par,link,dlink,ph)
V=mvnb_scorevar_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb2par,link,dlink,ph)
reg_nb2var_nam=solve(I)/n
reg_nb2var_rbm=solve(I)%*%V%*%solve(I)/n

#three kinds of tests
LR_reg_nb2=c()
reg_nb2llh=c()
walds_reg_nb2=c()
score_reg_nb2=c()
for(q in c(1:(p+1))){
    A=mvnb_A_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb2par,link,dlink,ph,q)
    B=mvnb_B_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb2par,link,dlink,ph,q)
    
    #LR test 
    max_f0=function(beta0){
        beta=rep(0,(p+1))
        beta[-q]=beta0
        ans=-mvnb_loglike_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,ph)
        return(ans)
    }
    opt0=optim(c(0,0,0,0,0),max_f0,method="Nelder-Mead",hessian = T)
    updatereg=0
    if(opt0$counts[1]>=501){
        repeat{
            updatereg=updatereg+1
            temppar=opt0$par
            opt0=optim(temppar,max_f0,method="Nelder-Mead",hessian = T)
            if(opt0$counts[1]<501 | updatereg==100)break
        }
    }
    beta_0_hat=opt0$par
    lr0=opt0$value
    LR_reg_nb2[q]=2*(A/B)*(lr0-lr)
    reg_nb2llh[q]=(A/B)*(-lr)
    
    #walds test
    w1=A%*%matrix(reg_nb2par[q]-0,1,1)
    walds_reg_nb2[q]=n*t(w1)%*%solve(B)%*%w1
    
    #score test
    beta=rep(0,(p+1))
    beta[-q]=beta_0_hat
    sA=mvnb_A_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    sB=mvnb_B_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    s1=mvnb_score_beta_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    score_reg_nb2[q]=1/n*t(s1)%*%solve(sB)%*%s1
}
reg_nb2var_na=diag(reg_nb2var_nam)
reg_nb2var_rb=diag(reg_nb2var_rbm)
reg_nb2p_w=1-pchisq(walds_reg_nb2,1)
reg_nb2p_s=1-pchisq(score_reg_nb2,1)
reg_nb2p_lr=1-pchisq(LR_reg_nb2,1)
CInb2L=c(reg_nb2par-qnorm(0.975)*sqrt(reg_nb2var_rb))
CInb2U=c(reg_nb2par+qnorm(0.975)*sqrt(reg_nb2var_rb))

#result
#robust 1,LR test, 2,score test 1,Wald test, p-value 
reg_nb2p_lr
reg_nb2p_s
reg_nb2p_w
#estimates of parameters
reg_nb2par
#estimates of standard deviation(robust)
as.numeric(sqrt(reg_nb2var_rb))
as.numeric(sqrt(reg_nb2var_na))
#Confident interval of Wald test(robust)
as.numeric(CInb2L)
as.numeric(CInb2U)
#==============================================================================

#====================================trunk=====================================
d=3

#log likelihood function
max_f=function(beta){
    -mvnb_loglike_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,ph)
}

#optimize
reg_nb3=optim(c(0,0,0,0,0,0),max_f,method="Nelder-Mead",hessian = T)
updatereg=0
if(reg_nb3$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        temppar=reg_nb3$par
        reg_nb3=optim(temppar,max_f,method="Nelder-Mead",hessian = T)
        if(reg_nb3$counts[1]<501 | updatereg==100)break
    }
}
reg_nb3par=reg_nb3$par
lr=reg_nb3$value

#I hat, V hat, naive variance hat, robust variance hat
I=mvnb_fisher_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb3par,link,dlink,ph)
V=mvnb_scorevar_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb3par,link,dlink,ph)
reg_nb3var_nam=solve(I)/n
reg_nb3var_rbm=solve(I)%*%V%*%solve(I)/n

#three kinds of tests
LR_reg_nb3=c()
reg_nb3llh=c()
walds_reg_nb3=c()
score_reg_nb3=c()
for(q in c(1:(p+1))){
    A=mvnb_A_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb3par,link,dlink,ph,q)
    B=mvnb_B_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb3par,link,dlink,ph,q)
    
    #LR test 
    max_f0=function(beta0){
        beta=rep(0,(p+1))
        beta[-q]=beta0
        ans=-mvnb_loglike_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,ph)
        return(ans)
    }
    opt0=optim(c(0,0,0,0,0),max_f0,method="Nelder-Mead",hessian = T)
    updatereg=0
    if(opt0$counts[1]>=501){
        repeat{
            updatereg=updatereg+1
            temppar=opt0$par
            opt0=optim(temppar,max_f0,method="Nelder-Mead",hessian = T)
            if(opt0$counts[1]<501 | updatereg==100)break
        }
    }
    beta_0_hat=opt0$par
    lr0=opt0$value
    LR_reg_nb3[q]=2*(A/B)*(lr0-lr)
    reg_nb3llh[q]=(A/B)*(-lr)
    
    #walds test
    w1=A%*%matrix(reg_nb3par[q]-0,1,1)
    walds_reg_nb3[q]=n*t(w1)%*%solve(B)%*%w1
    
    #score test
    beta=rep(0,(p+1))
    beta[-q]=beta_0_hat
    sA=mvnb_A_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    sB=mvnb_B_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    s1=mvnb_score_beta_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    score_reg_nb3[q]=1/n*t(s1)%*%solve(sB)%*%s1
}
reg_nb3var_na=diag(reg_nb3var_nam)
reg_nb3var_rb=diag(reg_nb3var_rbm)
reg_nb3p_w=1-pchisq(walds_reg_nb3,1)
reg_nb3p_s=1-pchisq(score_reg_nb3,1)
reg_nb3p_lr=1-pchisq(LR_reg_nb3,1)
CInb3L=c(reg_nb3par-qnorm(0.975)*sqrt(reg_nb3var_rb))
CInb3U=c(reg_nb3par+qnorm(0.975)*sqrt(reg_nb3var_rb))

#result
#robust 1,LR test, 2,score test 1,Wald test, p-value 
reg_nb3p_lr
reg_nb3p_s
reg_nb3p_w
#estimates of parameters
reg_nb3par
#estimates of standard deviation(robust)
as.numeric(sqrt(reg_nb3var_rb))
as.numeric(sqrt(reg_nb3var_na))
#Confident interval of Wald test(robust)
as.numeric(CInb3L)
as.numeric(CInb3U)
#==============================================================================


#===================================android====================================
d=4

#log likelihood function
max_f=function(beta){
    -mvnb_loglike_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,ph)
}

#optimize
reg_nb4=optim(c(0,0,0,0,0,0),max_f,method="Nelder-Mead",hessian = T)
updatereg=0
if(reg_nb4$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        temppar=reg_nb4$par
        reg_nb4=optim(temppar,max_f,method="Nelder-Mead",hessian = T)
        if(reg_nb4$counts[1]<501 | updatereg==100)break
    }
}
reg_nb4par=reg_nb4$par
lr=reg_nb4$value

#I hat, V hat, naive variance hat, robust variance hat
I=mvnb_fisher_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb4par,link,dlink,ph)
V=mvnb_scorevar_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb4par,link,dlink,ph)
reg_nb4var_nam=solve(I)/n
reg_nb4var_rbm=solve(I)%*%V%*%solve(I)/n

#three kinds of tests
LR_reg_nb4=c()
reg_nb4llh=c()
walds_reg_nb4=c()
score_reg_nb4=c()
for(q in c(1:(p+1))){
    A=mvnb_A_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb4par,link,dlink,ph,q)
    B=mvnb_B_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb4par,link,dlink,ph,q)
    
    #LR test 
    max_f0=function(beta0){
        beta=rep(0,(p+1))
        beta[-q]=beta0
        ans=-mvnb_loglike_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,ph)
        return(ans)
    }
    opt0=optim(c(0,0,0,0,0),max_f0,method="Nelder-Mead",hessian = T)
    updatereg=0
    if(opt0$counts[1]>=501){
        repeat{
            updatereg=updatereg+1
            temppar=opt0$par
            opt0=optim(temppar,max_f0,method="Nelder-Mead",hessian = T)
            if(opt0$counts[1]<501 | updatereg==100)break
        }
    }
    beta_0_hat=opt0$par
    lr0=opt0$value
    LR_reg_nb4[q]=2*(A/B)*(lr0-lr)
    reg_nb4llh[q]=(A/B)*(-lr)
    
    #walds test
    w1=A%*%matrix(reg_nb4par[q]-0,1,1)
    walds_reg_nb4[q]=n*t(w1)%*%solve(B)%*%w1
    
    #score test
    beta=rep(0,(p+1))
    beta[-q]=beta_0_hat
    sA=mvnb_A_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    sB=mvnb_B_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    s1=mvnb_score_beta_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    score_reg_nb4[q]=1/n*t(s1)%*%solve(sB)%*%s1
}
reg_nb4var_na=diag(reg_nb4var_nam)
reg_nb4var_rb=diag(reg_nb4var_rbm)
reg_nb4p_w=1-pchisq(walds_reg_nb4,1)
reg_nb4p_s=1-pchisq(score_reg_nb4,1)
reg_nb4p_lr=1-pchisq(LR_reg_nb4,1)
CInb4L=c(reg_nb4par-qnorm(0.975)*sqrt(reg_nb4var_rb))
CInb4U=c(reg_nb4par+qnorm(0.975)*sqrt(reg_nb4var_rb))

#result
#robust 1,LR test, 2,score test 1,Wald test, p-value 
reg_nb4p_lr
reg_nb4p_s
reg_nb4p_w
#estimates of parameters
reg_nb4par
#estimates of standard deviation(robust)
as.numeric(sqrt(reg_nb4var_rb))
as.numeric(sqrt(reg_nb4var_na))
#Confident interval of Wald test(robust)
as.numeric(CInb4L)
as.numeric(CInb4U)
#==============================================================================


#===================================gynoid=====================================
d=5

#log likelihood function
max_f=function(beta){
    -mvnb_loglike_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,ph)
}

#optimize
reg_nb5=optim(c(0,0,0,0,0,0),max_f,method="Nelder-Mead",hessian = T)
updatereg=0
if(reg_nb5$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        temppar=reg_nb5$par
        reg_nb5=optim(temppar,max_f,method="Nelder-Mead",hessian = T)
        if(reg_nb5$counts[1]<501 | updatereg==100)break
    }
}
reg_nb5par=reg_nb5$par
lr=reg_nb5$value

#I hat, V hat, naive variance hat, robust variance hat
I=mvnb_fisher_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb5par,link,dlink,ph)
V=mvnb_scorevar_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb5par,link,dlink,ph)
reg_nb5var_nam=solve(I)/n
reg_nb5var_rbm=solve(I)%*%V%*%solve(I)/n

#three kinds of tests
LR_reg_nb5=c()
reg_nb5llh=c()
walds_reg_nb5=c()
score_reg_nb5=c()
for(q in c(1:(p+1))){
    A=mvnb_A_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb5par,link,dlink,ph,q)
    B=mvnb_B_reg_rem(matrix(YYY[,d],n,1),XXX,reg_nb5par,link,dlink,ph,q)
    
    #LR test 
    max_f0=function(beta0){
        beta=rep(0,(p+1))
        beta[-q]=beta0
        ans=-mvnb_loglike_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,ph)
        return(ans)
    }
    opt0=optim(c(0,0,0,0,0),max_f0,method="Nelder-Mead",hessian = T)
    updatereg=0
    if(opt0$counts[1]>=501){
        repeat{
            updatereg=updatereg+1
            temppar=opt0$par
            opt0=optim(temppar,max_f0,method="Nelder-Mead",hessian = T)
            if(opt0$counts[1]<501 | updatereg==100)break
        }
    }
    beta_0_hat=opt0$par
    lr0=opt0$value
    LR_reg_nb5[q]=2*(A/B)*(lr0-lr)
    reg_nb5llh[q]=(A/B)*(-lr)
    
    #walds test
    w1=A%*%matrix(reg_nb5par[q]-0,1,1)
    walds_reg_nb5[q]=n*t(w1)%*%solve(B)%*%w1
    
    #score test
    beta=rep(0,(p+1))
    beta[-q]=beta_0_hat
    sA=mvnb_A_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    sB=mvnb_B_reg_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    s1=mvnb_score_beta_rem(matrix(YYY[,d],n,1),XXX,beta,link,dlink,ph,q)
    score_reg_nb5[q]=1/n*t(s1)%*%solve(sB)%*%s1
}
reg_nb5var_na=diag(reg_nb5var_nam)
reg_nb5var_rb=diag(reg_nb5var_rbm)
reg_nb5p_w=1-pchisq(walds_reg_nb5,1)
reg_nb5p_s=1-pchisq(score_reg_nb5,1)
reg_nb5p_lr=1-pchisq(LR_reg_nb5,1)
CInb5L=c(reg_nb5par-qnorm(0.975)*sqrt(reg_nb5var_rb))
CInb5U=c(reg_nb5par+qnorm(0.975)*sqrt(reg_nb5var_rb))

#result
#robust 1,LR test, 2,score test 1,Wald test, p-value 
reg_nb5p_lr
reg_nb5p_s
reg_nb5p_w
#estimates of parameters
reg_nb5par
#estimates of standard deviation(robust)
as.numeric(sqrt(reg_nb5var_rb))
as.numeric(sqrt(reg_nb5var_na))
#Confident interval of Wald test(robust)
as.numeric(CInb5L)
as.numeric(CInb5U)
#==============================================================================
#------------------------------------------------------------------------------



#---------------------------------5d model MNB----------------------------------
#setting
dm=5
p=5
df=30
n=length(bfpY[,1])
ph=1
YYY=as.matrix(bfpY)
XXX=matrix(0,n*dm,p)
for(j in c(1:dm)){
    XXX[c(dm*(c(1:n))-(dm-j)),]=as.matrix(bfpX)
}

#initial point (from the result of univariate model)
initnb=c(reg_nb1par,reg_nb2par,reg_nb3par,reg_nb4par,reg_nb5par)

#log likelihood function
max_f=function(beta){
    -mvnb_loglike_reg_rem_full_id.p(YYY,XXX,beta,link,ph)
}
loglikevalue=max_f(initnb)
initnb
loglikevalue

#optimize
reg_nb5d=optim(initnb,max_f,method="Nelder-Mead",hessian = T)
reg_nb5d
reg_nb5d$par
reg_nb5d$value
loglikevalue=reg_nb5d$value
temppar=reg_nb5d$par
#To ensure it is maximum, repeat optim() until convergence
set.seed(52210)
updatereg=0
a=0
if(reg_nb5d$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        reg_nb5d=optim(c(temppar+runif(30,-0.00001,0.00001)),max_f,method="Nelder-Mead",hessian = T)
        print(c(reg_nb5d$value,(reg_nb5d$value-loglikevalue)))
        egvalue_nb=eigen(reg_nb5d$hessian)$value
        detect=length(egvalue_nb[egvalue_nb<=0])
        #first-check the hessian matrix is positive-define
        if(detect<=0){
            print(c("positive define"))
            #second-reach convergence tolerance
            if((reg_nb5d$value-loglikevalue)<0){
                print(c("update point"))
                if(abs(reg_nb5d$value-loglikevalue)<0.000001){
                    a=1
                }
                loglikevalue=reg_nb5d$value
                temppar=reg_nb5d$par
            }
        }
        if(a==1)break
        if(reg_nb5d$counts[1]<501 | updatereg==600)break
    }
}
#To ensure it is maximum, repeat optim() until convergence(double check)
if(reg_nb5d$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        reg_nb5d=optim(c(temppar),max_f,method="Nelder-Mead",hessian = T)
        print(c(reg_nb5d$value,(reg_nb5d$value-loglikevalue)))
        egvalue_nb=eigen(reg_nb5d$hessian)$value
        detect=length(egvalue_nb[egvalue_nb<=0])
        if(detect<=0){
            print(c("positive define"))
            if((reg_nb5d$value-loglikevalue)<=0){
                print(c("update point"))
                if(abs(reg_nb5d$value-loglikevalue)<0.000001){
                    a=1
                }
                loglikevalue=reg_nb5d$value
                temppar=reg_nb5d$par
            }
        }
        if(a==1)break
        if(reg_nb5d$counts[1]<501 | updatereg==300)break
    }
}
reg_nb5dpar=temppar
lr=loglikevalue

#I hat, V hat, naive variance hat, robust variance hat
I=mvnb_fisher_reg_rem_full_id.p(YYY,XXX,reg_nb5dpar,link,dlink,ph)
V=mvnb_scorevar_reg_rem_full_id.p(YYY,XXX,reg_nb5dpar,link,dlink,ph)
reg_nb5dvar_nam=solve(I)/n
reg_nb5dvar_rbm=solve(I)%*%V%*%solve(I)/n

##walds test
walds_reg_nb5d=c()
for(q in c(1:df)){
    A=mvnb_A_reg_rem_full_id.p(YYY,XXX,reg_nb5dpar,link,dlink,ph,q)
    B=mvnb_B_reg_rem_full_id.p(YYY,XXX,reg_nb5dpar,link,dlink,ph,q)
    w1=A%*%matrix(reg_nb5dpar[q]-0,1,1)
    walds_reg_nb5d[q]=n*t(w1)%*%solve(B)%*%w1
}
reg_nb5dvar_na=diag(reg_nb5dvar_nam)
reg_nb5dvar_rb=diag(reg_nb5dvar_rbm)
reg_nb5dp_w=1-pchisq(walds_reg_nb5d,1)
CInb5dL=c(reg_nb5dpar-qnorm(0.975)*sqrt(reg_nb5dvar_rb))
CInb5dU=c(reg_nb5dpar+qnorm(0.975)*sqrt(reg_nb5dvar_rb))

#result
#likelihood
lr
#robust Wald test p-value 
1-pchisq(walds_reg_nb5d,1)
#estimates of parameters
reg_nb5dpar
#estimates of variance & standard deviation(robust)
reg_nb5dvar_rb
sqrt(reg_nb5dvar_rb)
#Confident interval of Wald test(robust)
CInb5dL
CInb5dU

table_bpd_robust=matrix(c(reg_nb5dpar,sqrt(reg_nb5dvar_rb),(1-pchisq(walds_reg_nb5d,1)),
                          CInb5dL,CInb5dU,c(CInb5dU-CInb5dL)),df,6)
colnames(table_bpd_robust)=c("par","sd","pvalue","CI_L","CI_U","length")
row.names(table_bpd_robust)=c("beta01","beta11","beta21","beta31","beta41","beta51",
                              "beta02","beta12","beta22","beta32","beta42","beta52",
                              "beta03","beta13","beta23","beta33","beta43","beta53",
                              "beta04","beta14","beta24","beta34","beta44","beta54",
                              "beta05","beta15","beta25","beta35","beta45","beta55")
table_bpd_robust[1:6,]
table_bpd_robust[7:12,]
table_bpd_robust[13:18,]
table_bpd_robust[19:24,]
table_bpd_robust[25:30,]
#------------------------------------------------------------------------------


#---------------------------------2d model MNB---------------------------------
#(arms & trunk) => d=1,t=3 ; (legs & android) => d=2,t=4
#setting
dm=2
p=5
df=12
n=length(bfpY[,1])
ph=1
d=2
t=4
YYY=as.matrix(bfpY)[,c(d,t)]
XXX=matrix(0,n*dm,p)
for(j in c(1:dm)){
    XXX[c(dm*(c(1:n))-(dm-j)),]=as.matrix(bfpX)
}

#initial point (from the result of univariate model)
initnb=c(reg_nb1par,reg_nb2par,reg_nb3par,reg_nb4par,reg_nb5par)
initnb=initnb[c(c(1:(p+1))+(p+1)*(d-1),c(1:(p+1))+(p+1)*(t-1))]

#log likelihood function
max_f=function(beta){
    -mvnb_loglike_reg_rem_full_id.p(YYY,XXX,beta,link,ph)
}
loglikevalue=max_f(initnb)
initnb
#In (legs & android) case, use the initial point like this is better
if(d==2&t==4){
    initnb[2]=0
}
loglikevalue


reg_nb2d=optim(initnb,max_f,method="Nelder-Mead",hessian = T)
reg_nb2d
reg_nb2d$par
reg_nb2d$value
loglikevalue=reg_nb2d$value
temppar=reg_nb2d$par
#To ensure it is maximum, repeat optim() until convergence
set.seed(52210)
updatereg=0
a=0
if(reg_nb2d$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        reg_nb2d=optim(c(temppar+runif(df,-0.00001,0.00001)),max_f,method="Nelder-Mead",hessian = T)
        print(c(reg_nb2d$value,(reg_nb2d$value-loglikevalue)))
        egvalue_nb=eigen(reg_nb2d$hessian)$value
        detect=length(egvalue_nb[egvalue_nb<=0])
        #first-check the hessian matrix is positive-define
        if(detect<=0){
            print(c("positive define"))
            #second-reach convergence tolerance
            if((reg_nb2d$value-loglikevalue)<0){
                print(c("update point"))
                if(abs(reg_nb2d$value-loglikevalue)<0.000001){
                    a=1
                }
                loglikevalue=reg_nb2d$value
                temppar=reg_nb2d$par
            }
        }
        if(a==1)break
        if(reg_nb2d$counts[1]<501 | updatereg==300)break
    }
}
#To ensure it is maximum, repeat optim() until convergence(double check)
if(reg_nb2d$counts[1]>=501){
    repeat{
        updatereg=updatereg+1
        reg_nb2d=optim(c(temppar),max_f,method="Nelder-Mead",hessian = T)
        print(c(reg_nb2d$value,(reg_nb2d$value-loglikevalue)))
        egvalue_nb=eigen(reg_nb2d$hessian)$value
        detect=length(egvalue_nb[egvalue_nb<=0])
        if(detect<=0){
            print(c("positive define"))
            if((reg_nb2d$value-loglikevalue)<=0){
                print(c("update point"))
                if(abs(reg_nb2d$value-loglikevalue)<0.000001){
                    a=1
                }
                loglikevalue=reg_nb2d$value
                temppar=reg_nb2d$par
            }
        }
        if(a==1)break
        if(reg_nb2d$counts[1]<501 | updatereg==300)break
    }
}
reg_nb2dpar=temppar
lr=loglikevalue

#I hat, V hat, naive variance hat, robust variance hat
I=mvnb_fisher_reg_rem_full_id.p(YYY,XXX,reg_nb2dpar,link,dlink,ph)
V=mvnb_scorevar_reg_rem_full_id.p(YYY,XXX,reg_nb2dpar,link,dlink,ph)
reg_nb2dvar_nam=solve(I)/n
reg_nb2dvar_rbm=solve(I)%*%V%*%solve(I)/n

#three kinds of tests
LR_reg_nb2d=c()
walds_reg_nb2d=c()
score_reg_nb2d=c()
for(q in c(1:df)){
    A=mvnb_A_reg_rem_full_id.p(YYY,XXX,reg_nb2dpar,link,dlink,ph,q)
    B=mvnb_B_reg_rem_full_id.p(YYY,XXX,reg_nb2dpar,link,dlink,ph,q)
    
    #LR test 
    max_f0=function(beta0){
        beta=rep(0,(p+1)*dm)
        beta[-q]=beta0
        ans=-mvnb_loglike_reg_rem_full_id.p(YYY,XXX,beta,link,ph)
        return(ans)
    }
    opt0=optim(reg_nb2dpar[-q],max_f0,method="Nelder-Mead",hessian = T)
    updatereg=0
    if(opt0$counts[1]>=501){
        repeat{
            print(c(q,"updating",opt0$value))
            updatereg=updatereg+1
            temppar=opt0$par
            opt0=optim(temppar,max_f0,method="Nelder-Mead",hessian = T)
            if(opt0$counts[1]<501 | updatereg==100)break
        }
    }
    beta_0_hat=opt0$par
    lr0=opt0$value
    LR_reg_nb2d[q]=2*(A/B)*(lr0-lr)
    
    #walds test
    w1=A%*%matrix(reg_nb2dpar[q]-0,1,1)
    walds_reg_nb2d[q]=n*t(w1)%*%solve(B)%*%w1
    
    #score test
    beta=rep(0,(p+1)*dm)
    beta[-q]=beta_0_hat
    sA=mvnb_A_reg_rem_full_id.p(YYY,XXX,beta,link,dlink,ph,q)
    sB=mvnb_B_reg_rem_full_id.p(YYY,XXX,beta,link,dlink,ph,q)
    s1=mvnb_score_beta_rem_full_id.p(YYY,XXX,beta,link,dlink,ph,q)
    score_reg_nb2d[q]=1/n*t(s1)%*%solve(sB)%*%s1
    
    print(c(q,"finish"))
}
reg_nb2dvar_na=diag(reg_nb2dvar_nam)
reg_nb2dvar_rb=diag(reg_nb2dvar_rbm)
reg_nb2dp_w=1-pchisq(walds_reg_nb2d,1)
reg_nb2dp_s=1-pchisq(score_reg_nb2d,1)
reg_nb2dp_lr=1-pchisq(LR_reg_nb2d,1)
CInb5dL=c(reg_nb2dpar-qnorm(0.975)*sqrt(reg_nb2dvar_rb))
CInb5dU=c(reg_nb2dpar+qnorm(0.975)*sqrt(reg_nb2dvar_rb))

#result
#likelihood
lr
#robust 1,LR test, 2,score test 1,Wald test, p-value 
1-pchisq(walds_reg_nb2d,1)
1-pchisq(score_reg_nb2d,1)
1-pchisq(LR_reg_nb2d,1)
#estimates of parameters
reg_nb2dpar
#estimates of variance & standard deviation(robust)
reg_nb2dvar_rb
sqrt(reg_nb2dvar_rb)
#Confident interval of Wald test(robust)
CInb5dL
CInb5dU

table_bpd_robust_2d=matrix(c(reg_nb2dpar,sqrt(reg_nb2dvar_rb),(1-pchisq(walds_reg_nb2d,1)),
                             CInb5dL,CInb5dU,c(CInb5dU-CInb5dL)),df,6)
colnames(table_bpd_robust_2d)=c("par","sd","pvalue","CI_L","CI_U","length")
row.names(table_bpd_robust_2d)=c("beta0d","beta1d","beta2d","beta3d","beta4d","beta5d",
                                 "beta0t","beta1t","beta2t","beta3t","beta4t","beta5t")
table_bpd_robust_2d[1:6,]
table_bpd_robust_2d[7:12,]
#------------------------------------------------------------------------------




