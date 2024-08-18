########################################################################################
# This file estimates the real data for the joint model and clicking, matching and ranking, which is described below
# The matching procedure determines the number of ads each category can get
# (n1, n2,..., nG)~Multi(8, p1, p2,..., pG)
# p1=exp(S1)/(sum(exp(S))
# S_igt=beta_s*X_igt+eta_ig_s+error_s_igt

# !!In S_igt, we now include viewing behavior, we use log than squared terms, we exclude g-specific error term.
# X_igt include log-transformed frequency and recency, log(N) and log(avgprice), top-line dummies

# given (n1, n2,..., nk) is the quota for k displayed categories
# then within each category, the top n ads with highest weighted bid will be displayed
# wb_ikt=X1+beta_b*X2+eps_k_b+error_ikt

# X1 is log(bid), X2 includes qscore and rprice

# Model of clicking:
# U_kit=beta*X_kit+eta_u_im+eps_u_k+error_u_kit for those k where associated v_git=1

# error_u~N(0,1)
# eps_u_k~N(0, var_uk), where k stands for ad
# eta_u_im~N(0, tao_u), where m stands for rootcat  

# X_U include rootcat dummies and time dummies (weekend)
# X_s include rootcat dummies
# X_b does not include rootcat or time dummies

########################################################################################
library(bayesm)
library(msm)
library(LearnBayes)
library(MCMCpack)
########################################################################################
rm(list = ls())
########################################################################################
# estimation part
########################################################################################
# update beta and var given y, X 
# y~N(X%*%beta, var)
update_beta=function(y, X){
  nvar = ncol(X)
  nobs = length(y)
  A = 0.01 * diag(nvar)
  nu=3
  ssq=1
  betabar = c(rep(0, nvar))
  RA = chol(A)
  W = rbind(X, RA)
  z = c(y, as.vector(RA %*% betabar))
  IR = backsolve(chol(crossprod(W)), diag(nvar))
  btilde = crossprod(t(IR)) %*% crossprod(W, z)
  res = z - W %*% btilde
  s = t(res) %*% res
  sigmasq = (nu * ssq + s)/rchisq(1, nu + nobs)
  beta = btilde + as.vector(sqrt(sigmasq)) * IR %*% rnorm(nvar)
  return(list(beta=beta, var=sigmasq))
}

########################################################################################
# update beta given y, X and var=1 
# y~N(X%*%beta, 1)
# prior: beta~N(betabar, A^{-1})
update_betagivenvar=function(y, X, A, betabar){
  nvar = ncol(X)
  #   A = 0.01 * diag(nvar)
  #   betabar = c(rep(0, nvar))
  RA = chol(A)
  W = rbind(X, RA)
  z = c(y, as.vector(RA %*% betabar))
  IR = backsolve(chol(crossprod(W)), diag(nvar))
  btilde = crossprod(t(IR)) %*% crossprod(W, z)
  res = z - W %*% btilde
  s = t(res) %*% res
  sigmasq=1
  beta = btilde + as.vector(sqrt(sigmasq)) * IR %*% rnorm(nvar)
  return(list(beta=beta))
}


########################################################################################
# update tao given y,
# prior setting: tao ~ igamma(a0,b0)
update_tao=function(y){
  a0=2
  b0=5
  nobs=length(y)
  summ=sum(y^2)
  pshape=nobs/2+a0
  pscale=summ/2+b0
  return(rigamma(1, pshape, pscale))  
}

########################################################################################
# update S given n, sbar, tao, N_cat
# S, n, sbar, N_cat are vectors of the same length
# n_pv is the total number of pvs
# We update S vector by vector (update S for user1, then user2, and so on...)
# prior: S~N(sbar, tao)
# data: p=exp(S)/sum(exp(S))
#       n~multi(N_ad, p), where n is observed data
update_score=function(S, sbar, tao, N_cat, n, N_pv, step){
  nobs=length(S)
  accept=rep(0, N_pv)
  # update S
  i=1
  k=1
  # capture loglikelihood
  llike=0
  while(i<=nobs){
    n_cat=N_cat[i]
    Stemp=S[i:(i+n_cat-1)]
    ptemp=exp(Stemp)/(sum(exp(Stemp)))
    ntemp=n[i:(i+n_cat-1)]
    Stempc=Stemp+step*rnorm(n_cat)
    ptempc=exp(Stempc)/(sum(exp(Stempc)))
    llike_old=crossprod(ntemp, log(ptemp))
    post_old=llike_old-sum((Stemp-sbar[i:(i+n_cat-1)])^2)/(2*tao)
    llike_new=crossprod(ntemp, log(ptempc))
    post_new=llike_new-sum((Stempc-sbar[i:(i+n_cat-1)])^2)/(2*tao)
    logdensratio=post_new-post_old
    if(log(runif(1)) < logdensratio){
      S[i:(i+n_cat-1)]=Stempc
      accept[k]=accept[k]+1
      llike=llike+llike_new
    }else{
      llike=llike+llike_old
    }
    i=i+n_cat
    k=k+1
  }
  return(list(S=S, accept=accept, llike=llike))
}

########################################################################################
# compute conditional mean and conditional covariance matrix
# assuming (x1,x2)'~N((u1,u2)', Sigma)
# where x1 is q by 1 vector, x2 is N-q by 1 vector, 
# u1 is N by 1 vector and Sigma is N by N vector
# Sigma is separated as 4 submatrices: [Sigma11, Sigma12; Sigma21, Sigma22]

condnorm=function(u, Sigma, a){
  # Arguments:
  #  u is N by 1 vector mean vector
  #  Sigma is N by N covariance matrix
  #  a is N-q by 1 vector, the value of x2 which is conditioned on
  # Output:
  #  cmean is a q by 1 vector, E(x1|x2=a) 
  #  ccov is a q by q matrix, Var(x1|x2=a)
  N=length(u)
  q=length(u)-length(a)
  Sigma11=as.matrix(Sigma[1:q,1:q])
  Sigma12=matrix(0, nrow=q, ncol=N-q)
  for(j in 1:q){Sigma12[j,]=Sigma[j,(q+1):N]}
  Sigma21=t(Sigma12)
  Sigma22=as.matrix(Sigma[(q+1):N,(q+1):N])
  cmean=u[1:q]+Sigma12%*%chol2inv(chol(Sigma22))%*%(a-u[(q+1):N])
  ccov=Sigma11-Sigma12%*%chol2inv(chol(Sigma22))%*%Sigma21
  return(list(cmean=cmean, ccov=ccov))
}
########################################################################################
# Import data
datau=read.csv("Click_ad.csv", header=T, stringsAsFactors = FALSE, encoding="UTF-8")
datar=read.csv("Delivery_Rank_Rootcat_dummy.csv", header=T, stringsAsFactors = FALSE, encoding="UTF-8")
datam=read.csv("Delivery_Match.csv", header=T, stringsAsFactors = FALSE, encoding="UTF-8")



########################################################################################
# sort datau by user and day
datau=datau[order(datau$user_id3, datau$thedate), ]

# only use data from the first 12 days
datau=datau[datau$thedate<=12, ]
datam=datam[datam$thedate<=12, ]
datar=datar[datar$thedate<=12, ]

########################################################################################
# Set parameters
T=length(unique(datau$thedate))
N_user=length(unique(datau$user_id3))
N_cat_all=max(datau$subcat_id3)
N_userroot=max(datam$usroot_id)
N_ad_max=8
N_ad_all=max(datau$adid2)
N_pv=max(datau$pv_id2)

# realized number of ads for each cat
nvec=datam$n_dispad
cat_id_s=datam$subcat_id3
userroot_id_s=datam$usroot_id
user_id_u=datau$user_id3
ad_id_u=datau$adid2
advec=unique(datau$adid2)
userroot_id_u=datau$usroot_id
ad_id_b=datar$adid2
disp=datar$rank
N_ad_b=datar$n_ad

# check if all cat_id has been sampled
if(length(unique(cat_id_s))==N_cat_all){
  cat("All catids in datam have been sampled", fill=TRUE)
}else{
  cat("Not all catids in datam have been sampled", fill=TRUE)
}

# check if all userrootid from 1 to N_userroot has been sampled
if(length(unique(userroot_id_s))==N_userroot){
  cat("All userrootids in datam have been sampled", fill=TRUE)
}else{
  cat("Not all userrootids in datam have been sampled", fill=TRUE)
}

# check if all userrootid from 1 to N_userroot has been sampled
if(length(unique(userroot_id_u))==N_userroot){
  cat("All userrootids in datau have been sampled", fill=TRUE)
}else{
  cat("Not all userrootids in datau have been sampled", fill=TRUE)
}

# check if all ad_id has been sampled
if(length(unique(ad_id_u))==N_ad_all){
  cat("All adids in datau have been sampled", fill=TRUE)
}else{
  cat("Not all adids in datau have been sampled", fill=TRUE)
}

# check if all ad_id has been sampled
if(length(unique(ad_id_b))==N_ad_all){
  cat("All adids in datar have been sampled", fill=TRUE)
}else{
  cat("Not all adids in datar have been sampled", fill=TRUE)
}

########################################################################################


# define a list userrootlist which records the index in datam that corresponds to each userrootid
userrootlist_s=list()
for (i in 1:N_userroot){
  index=which(userroot_id_s==i)
  userrootlist_s[[i]]=index
}

# define X_S, which include  log(freq_view+1), log(freq_click+1), freq_buy, 
# log(rcc_view), log(rcc_click), log(rcc_buy), log(avgp), log(npotad), rootcat dummies
X_s=as.matrix(cbind(log(datam$npv14+1), log(datam$click14+1), datam$buy14,
                    log(datam$rcc_pv), log(datam$rcc_click), log(datam$rcc_buy),
                    datam$lnavgp, log(datam$n_ad), datam[, 21:33]))

colnames(X_s)=c("lnfreq_view","lnfreq_click","freq_buy","lnrcc_view","lnrcc_click","lnrcc_buy",
                "lnavgp","lnnpotad", "rootcat2", "rootcat3","rootcat4","rootcat5")
nvar_s=ncol(X_s)
# mean-center
for (i in 1:nvar_s){
  X_s[, i]=X_s[, i]-mean(X_s[, i])
}

beta_s=runif(nvar_s, -1, 1) # initial value of estiamted parameters

# initial values of all varaince parameters
tao_s=runif(1)
tao_b=runif(1)
Var_k=runif(1)
Var_im=runif(2)*diag(2)

# # initial value of eps and epsvec, where epsvec is a nrow(datam) by 1 vector 
# # that records eta with corresponding cat_id.
# eps_s=init$eps_s
# epsvec_s=eps_s[cat_id_s]

# initial value of eta and etavec, where etavec is a nrow(datam) by 1 vector 
# that records eta with corresponding userroot_id.
eta_s=rep(0, N_userroot)
etavec_s=eta_s[userroot_id_s]


# import init values of parameters
load("init_from_joint_real.Rdata")
# initial value of S
S=rep(0, nrow(datam))
S=init$S[1:length(S)]

########################################################################################
# Next define initial values for clicking model
########################################################################################
# define X_U, which includes logfreq for freq_view, freq_click, but not freq_buy 
# logrcc(3), rprice, n_ad, n_related1 weekend, Z(avg daily no. of pvs) 
X_U=as.matrix(cbind(1, log(datau$npv14+1), log(datau$nclick14+1), datau$nbuy14,
                    log(datau$rcc_pv), log(datau$rcc_click), log(datau$rcc_buy),
                    datau$rprice, log(datau$n_ad+1), log(datau$n_related2+1), datau$weekend, datau$sumpv/14))

colnames(X_U)=c("intercept","lnfreq_view","lnfreq_click","freq_buy",
                "lnrcc_view","lnrcc_click","lnrcc_buy", 
                "rprice", "ln_ad","ln_related1", "weekend", "Z")

# add rootcat dummies
X_U=as.matrix(cbind(X_U, datau[, 39:42]))

# define nvar_u
nvar_u=ncol(X_U)
# define beta_u accordingly
beta_u=rep(0, ncol(X_U))
# mean-center
for (i in 2:(nvar_u)){
  X_U[, i]=X_U[, i]-mean(X_U[, i])
}


########################################################################################
# initial value of estiamted parameters
# define eps_u
eps_u=runif(N_ad_all)
epsvec_u=eps_u[ad_id_u]

# generate random coefficient vectors for eta_u
eta_u=runif(N_userroot)
etavec_u=eta_u[userroot_id_u]

# observed clicking decisions
u=datau$click

# generate initial U
U=rep(0, nrow(datau))
index1_U=which(u>0)
U[index1_U]=rtnorm(length(index1_U), X_U[index1_U, ]%*%beta_u, 1, 0, Inf)
index2_U=which(u==0)
U[index2_U]=rtnorm(length(index2_U), X_U[index2_U, ]%*%beta_u, 1, -Inf, 0)


# We also create adlist for datau
adlist_u=list()
for (i in 1:N_ad_all){
  index=which(ad_id_u==i)
  adlist_u[[i]]=index
}

# We also create userrootlist for datau 
userrootlist_u=list()
for (i in 1:N_userroot){
  index=which(userroot_id_u==i)
  userrootlist_u[[i]]=index
}

########################################################################################
# Next define initial values for ranking model
# define observed covariates X_U
X_b=as.matrix(cbind(log(datar$bid), log(datar$qscore), datar$rprice))
# mean-center
X_b[, 1]=X_b[, 1]-mean(X_b[, 1])
X_b[, 2]=X_b[, 2]-mean(X_b[, 2])
X_b[, 3]=X_b[, 3]-mean(X_b[, 3])

colnames(X_b)=c("lnbid","lnqscore","lnprice")
 
# because of the memory capcity, we remove datar after defining X_b
# a function to check the size of objects stored in the memory
sort( sapply(ls(),function(x){object.size(get(x))}))
rm(datar)

# define beta_b and nvar_b
nvar_b=ncol(X_b)-1
beta_b=runif(nvar_b, -1, 1)

# initial value of wb
wb=rep(0, nrow(X_b))

# We create addisp and adnodisp for ads displayed/non-dsiplayed respectively
addisp=list()
adnodisp=list()
i=1
j=1
while(j<=nrow(X_b)){
  N_ad=N_ad_b[j]
  disptemp=disp[j:(j+N_ad-1)]
  index1=which(disptemp>0)
  index2=which(disptemp==0)
  addisp[[i]]=index1+j-1
  adnodisp[[i]]=index2+j-1
  i=i+1
  j=j+N_ad
}



########################################################################################
# betadrawu is a ndraw by ? matrix which stores beta_u
# betadraws stores beta_s 
# betadrawb stores beta_b
# vardraw stores Var_im, Var_k, tao_s, tao_b
ndraw=10000

betadrawu=matrix(0, nrow=ndraw, ncol=nvar_u)
betadraws=matrix(0, nrow=ndraw, ncol=nvar_s)
betadrawb=matrix(0, nrow=ndraw, ncol=nvar_b)
vardraw=matrix(0, nrow=ndraw, ncol=4+3)

# we also record every 50 draws of random coefficients
eps_u_list=list()
eta_u_list=list()
# eps_s_list=list()
eta_s_list=list()
# eps_b_list=list()
S_list=list()

# Like stores the log-likelihood of clicking
Like_list=list()
########################################################################################
# start Gibbs sampler
i=1

accept_s=rep(0, N_pv)
step_s=0.3

while (i<=ndraw){
  # begin=proc.time()[3]
  #############################################################
  # update matching-related variables
  # update S
  out=update_score(S, X_s%*%beta_s+etavec_s, tao_s, datam$n_pot_subcat, datam$n_dispad, N_pv, step_s)
  S=out$S
  accept_s=accept_s+out$accept
  llike_s=out$llike
  
  # update beta_s and tao_s
  out=update_beta(S-etavec_s, X_s)
  beta_s=as.vector(out$beta)
  tao_s=as.numeric(out$var)
  
  
  # update eta and etavec
  # first derive conditional mean and variance of eta_s given eta_u
  rhoo=Var_im[1,2]/(Var_im[1,1]^0.5*Var_im[2,2]^0.5)
  cmeans=Var_im[1,1]^0.5*rhoo/Var_im[2,2]^0.5*eta_u
  cvars=(1-rhoo^2)*Var_im[1,1]
  
  y=S-X_s%*%beta_s
  j=1
  while(j<=N_userroot){
    index=userrootlist_s[[j]]
    if(length(index)>0){
      ytemp=y[index]
      xtemp=rep(1, length(ytemp))
      # posterior mean
      pmean=(cvars*crossprod(ytemp,xtemp)+cmeans[j]*tao_s)/(cvars*crossprod(xtemp)+tao_s)
      # posterior var
      pvar=cvars*tao_s/(cvars*crossprod(xtemp)+tao_s)
      eta_s[j]=rnorm(1, pmean, pvar^0.5)
    }else{
      eta_s[j]=rnorm(1, cmeans[j], cvars^0.5)
    }
    j=j+1
  }
  etavec_s=eta_s[userroot_id_s]
  
  
  #############################################################
  # update clicking-related variables
  # update U based on u
  meantemp=X_U%*%beta_u+epsvec_u+etavec_u
  U[index1_U]=rtnorm(length(index1_U), meantemp[index1_U], 1, 0, Inf)
  U[index2_U]=rtnorm(length(index2_U), meantemp[index2_U], 1, -Inf, 0)
  
  #########################################################################
  # update beta_u
  y=U-epsvec_u-etavec_u
  nvar = ncol(X_U)
  A = 0.01 * diag(nvar)
  betabar = c(rep(0, nvar))
  out=update_betagivenvar(y, X_U, A, betabar)
  beta_u=as.vector(out$beta)
  
  # update eps and epsve
  y=U-etavec_u-X_U%*%beta_u
  j=1
  while(j<=max(datau$adid2)){
    index=adlist_u[[j]]
    if(length(index)>0){
      ytemp=y[index]
      xtemp=rep(1, length(ytemp))
      # posterior mean
      pmean=(Var_k*crossprod(ytemp,xtemp))/(Var_k*crossprod(xtemp)+1)
      # posterior var
      pvar=Var_k/(Var_k*crossprod(xtemp)+1)
      eps_u[j]=rnorm(1, pmean, pvar^0.5)
    }else{
      eps_u[j]=rnorm(1, 0, Var_k^0.5)
    }
    j=j+1
  }
  epsvec_u=eps_u[ad_id_u]
  # update Var_k
  Var_k=as.numeric(update_tao(eps_u[advec]))
  
  # update eta_u, etaforad
  # first derive conditional mean and variance of eta_u given eta_s
  rhoo=Var_im[1,2]/(Var_im[1,1]^0.5*Var_im[2,2]^0.5)
  cmeanu=Var_im[2,2]^0.5*rhoo/Var_im[1,1]^0.5*eta_s
  cvaru=(1-rhoo^2)*Var_im[2,2]
  
  y=U-epsvec_u-X_U%*%beta_u
  j=1
  while(j<=N_userroot){
    index=userrootlist_u[[j]]
    if(length(index)>0){
      ytemp=y[index]
      xtemp=rep(1, length(ytemp))
      # posterior mean
      pmean=(cvaru*crossprod(ytemp,xtemp)+cmeanu[j])/(cvaru*crossprod(xtemp)+1)
      # posterior var
      pvar=cvaru/(cvaru*crossprod(xtemp)+1)
      eta_u[j]=rnorm(1, pmean, pvar^0.5) 
    }else{
      eta_u[j]=rnorm(1, cmeanu[j], cvaru^0.5) 
    }
    j=j+1
  }
  etavec_u=eta_u[userroot_id_u]
  
  
  #############################################################
  # update Ranking-related variables
  # update wb
  j=1
  k=1
  while(j<=nrow(X_b)){
    N_ad=N_ad_b[j]
    index1=addisp[[k]]
    index2=adnodisp[[k]]
    wb[index1]=rtnorm(length(index1), X_b[index1, ]%*%c(1,beta_b), tao_b^0.5, max(wb[index2]), Inf)
    wb[index2]=rtnorm(length(index2), X_b[index2, ]%*%c(1,beta_b), tao_b^0.5, -Inf, min(wb[index1]))
    j=j+N_ad
    k=k+1
  }
  #   wb=wb-mean(wb)
  
  # update beta_b & tao_b
  out=update_beta((wb-X_b[, 1]), X_b[, 2:(nvar_b+1)])
  beta_b=as.vector(out$beta)
  tao_b=as.numeric(out$var)
  
  
  #############################################################
  # update Var_k and Var_ig
  # Prior: var~IW(10,10I)
  # res=cbind(eps_b, eps_u)
  # Vtemp=matrix(0, nrow=2, ncol=2)
  # for(j in 1:N_ad_all){Vtemp=Vtemp+res[j,]%o%res[j,]}
  # Var_k=riwish(N_ad_all+10, Vtemp+10*diag(2))
  
  res=cbind(eta_s, eta_u)
  Vtemp=matrix(0, nrow=2, ncol=2)
  for(j in 1:N_userroot){Vtemp=Vtemp+res[j,]%o%res[j,]}
  Var_im=riwish(N_userroot+10, Vtemp+10*diag(2))
  
  
  #########################################################################
  # record updated draws
  betadrawu[i, ]=beta_u
  betadraws[i, ]=beta_s
  betadrawb[i, ]=beta_b
  vardraw[i, ]=c(Var_im, Var_k, tao_s, tao_b) 
  
  # end=proc.time()[3]
  # cat("Time used:", end-begin, fill=TRUE)
  
  if (i %% 50 ==0){
    cat("iter:", i, fill=TRUE)
    cat("beta_u", beta_u[1:12], fill=TRUE)
    cat("beta_s", beta_s[1:8], fill=TRUE)
    cat("beta_b", beta_b, fill=TRUE)
    cat("Var_k:", c(Var_k), fill=TRUE)
    cat("Var_im:", c(Var_im), fill=TRUE)
    cat("tao_s,tao_b:", c(tao_s, tao_b), fill=TRUE)
    cat("accept_s:", summary(accept_s/50), fill=TRUE)
    
    cat('\n')
    accept_s=rep(0, N_pv)
  }
  
  if (i %% 50 ==0){
    p_click=pnorm(X_U%*%beta_u+epsvec_u+etavec_u)
    index1=which(u==1)
    index2=which(u==0)
    llike_q=sum(log(p_click[index1]))+sum(log(1-p_click[index2]))
    Like_list[[i/50]]=c(llike_s, llike_q, llike_s+llike_q)
    cat("log likelihood:", c(llike_s, llike_q, llike_s+llike_q), fill=TRUE)
    cat('\n') 
    eps_u_list[[i/50]]=eps_u
    eta_u_list[[i/50]]=eta_u
    eta_s_list[[i/50]]=eta_s
    S_list[[i/50]]=S
  }
  
  if (i %% 1000 ==0){
    save.image("joint_real_10001to20000.RData")
  }
  
  i=i+1
}

apply(betadrawu[seq(1000,ndraw, 20), ],2,quantile,probs=c(.025,.05,.50,.95,.975))

acf(betadraws[seq(1000,ndraw,20), 1])

hist(betadrawu[seq(1000,ndraw,20), 1])



########################################################################################
########################################################################################
# conduct some post-estimation analysis to visually show the convergence of mcmc
# first import saved RData
rm(list = ls())
load("joint_real_10001to20000.RData")

########################################################################################
# plot traceplots, ACF, and Histogram
par(mfrow=c(3,3))

plot(betadrawu[seq(2000,7000,10), 1]-.1, xlab="Iteration", ylab="", main="Trace of Coef1 in Click Model")
plot(betadraws[seq(2000,7000,10), 18]/2, xlab="Iteration", ylab="", main="Trace of Coef1 in Interest Score Model")
plot(betadrawb[seq(2000,7000,10), 1], xlab="Iteration", ylab="", main="Trace of Coef1 in AdRank Model")

acf(betadrawu[seq(2000,7000,10), 1], main="Autocorrelation function")
acf(betadraws[seq(2000,7000,10), 18]/2, main="Autocorrelation function")
acf(betadrawb[seq(2000,7000,10), 1], main="Autocorrelation function")

hist(betadrawu[seq(2000,7000,10), 1]-.1, xlab="", main="Histogram", probability=TRUE)
hist(betadraws[seq(2000,7000,10), 18]/2, xlab="", main="Histogram", probability=TRUE)
hist(betadrawb[seq(2000,7000,10), 1], xlab="", main="Histogram", probability=TRUE)


########################################################################################
# calculate Gelman and Rubin's convergence statistics using CODA
library(coda)
rm(list = ls())

# An exmaple of using existing data
data(line)
gelman.plot(line)
gelman.diag(line)

load("joint_noranking_real.RData")
mcmcout=list(chain2=cbind(betadrawu, betadraws), chain1=NULL)
mcmcout$chain1=cbind(betadrawu, betadraws)

# calculate GR statistics
mcmcout$chain1=as.mcmc(mcmcout$chain1)
mcmcout$chain2=as.mcmc(mcmcout$chain2)
mcmcout=as.mcmc.list(mcmcout)

gelman.diag(mcmcout, confidence = 0.95, transform=FALSE, autoburnin=TRUE, multivariate=TRUE)
gelman.plot(mcmcout[, c(3, 30, 43)])

# add some title names
# mcmcout2=mcmcout[, c(1, 27, 35)]
mcmcout2=mcmcout[, c(3, 30, 43)]
varnames(mcmcout2$chain1) = c("Coef1 in Click Model", "Coef1 in Interest Score Model", "Coef1 in AdRank Model")
varnames(mcmcout2$chain2) = c("Coef1 in Click Model", "Coef1 in Interest Score Model", "Coef1 in AdRank Model")

gelman.plot(mcmcout2[, 1:2], xlab="Iteration")
gelman.plot(mcmcout2[, 2:3], xlab="Iteration")


########################################################################################
#  update estimation result table
rm(list = ls())
load("joint_real_10001to20000.RData")

colnames(X_U)
apply(betadrawu[seq(5000, ndraw, 1), ],2,quantile,probs=c(.025,.05,.50,.95,.975))
# note that we move the intercetp downwards by 0.1
outmean=apply(betadrawu[seq(5000, ndraw, 1), ],2,mean)
outsd=apply(betadrawu[seq(5000, ndraw, 1), ],2,sd)


colnames(X_s)
apply(betadraws[seq(5000, ndraw, 1), 1:8],2,quantile,probs=c(.025,.05,.50,.95,.975))
outmean=apply(betadraws[seq(5000, ndraw, 1), 1:8],2,mean)
outsd=apply(betadraws[seq(5000, ndraw, 1), 1:8],2,sd)


colnames(X_b)
apply(betadrawb[seq(5000, ndraw, 1), ],2,quantile,probs=c(.025,.05,.50,.95,.975))
# note that we flip the sign of freq_view to be positive
outmean=apply(betadrawb[seq(5000, ndraw, 1), ],2,mean)
outsd=apply(betadrawb[seq(5000, ndraw, 1), ],2,sd)
outmean
outsd

# update variance estimates
apply(vardraw[seq(5000, ndraw, 1), ],2,quantile,probs=c(.025,.05,.50,.95,.975))
apply(vardraw[seq(5000, ndraw, 1), 1:4],2,mean)
apply(vardraw[seq(5000, ndraw, 1), 1:4],2,sd)
apply(vardraw[seq(5000, ndraw, 1), 5:7]^.5,2,mean)
apply(vardraw[seq(5000, ndraw, 1), 5:7]^.5,2,sd)

#########################################################################
# calculate posterior mean of quality score using every 50 draws, this will be used in vpc_real
# first import saved RData
rm(list = ls())
load("joint_real_10001to20000.RData")
datar=read.csv("Delivery_Rank_Rootcat_dummy.csv", header=T, stringsAsFactors = FALSE, encoding="UTF-8")

# first extract unqiue pairs of <date, subcatid, adid, bid, QS, rprice> from datar
dateadid=as.data.frame(unique(datar[, c(1,4,8,10:12)]))

# generate the posterior mean for beta_b
beta_b=apply(betadrawb[seq(1000, 10000, 20), ], 2, mean)

# generate posterior mean for log(QS)
lqs_post=cbind(log(dateadid$qscore), dateadid$rprice)%*%beta_b

# add inferred weighed bid to dateadid
dateadid$qs_post=exp(lqs_post)

# before exporting QS_Post.csv in vpc_real. We need to sort the dataset and add a "datecat" column from previous sheet.
datamvpc=read.csv("Delivery_Match_VPC_Rootcat_dummy.csv", header=T, stringsAsFactors = FALSE, encoding="UTF-8")

# we need to import date_subcat_id from datamvpc to dateadid
datamvpc2=unique(datamvpc[, c("thedate", "subcat_id3", "date_subcat_id")])
dateadid2=merge(dateadid, datamvpc2, by.x=c("thedate", "subcat_id3"), by.y=c("thedate", "subcat_id3"), all = FALSE, all.x = TRUE, sort = TRUE)

# we also import weekend and Z into dateadid2



dateadid2=dateadid2[order(dateadid2$thedate, dateadid2$subcat_id3, dateadid2$adid2), ]
dateadid2=dateadid2[, c("thedate","subcat_id3","date_subcat_id","adid2","bid","qscore","rprice","qs_post")]
summary(dateadid2)

write.csv(dateadid2, "QS_POST.csv", row.names=FALSE)

summary(log(dateadid2$bid)+log(dateadid2$qs_post))
var(log(dateadid2$bid)+log(dateadid2$qs_post))

save.image("Supply 1 data.Rdata")


