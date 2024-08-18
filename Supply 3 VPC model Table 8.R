########################################################################################
# The purpose of this file is regress inferred vpc on advertiser characteristics
# log(vpc)=beta*X+eta_k+eps
# eps~N(0,tao)
# eta_k~N(0, var_k)

library(stargazer)



########################################################################################
library(bayesm)
library(msm)
library(LearnBayes)
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
  nu=10
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
# update beta given y, X and tao 
# y~N(X%*%beta, tao)
# prior: beta~N(betabar, A^{-1})
update_betagivenvar=function(y, X, A, betabar, tao){
  nvar = ncol(X)
  RA = chol(A)
  W = rbind(X, RA)
  z = c(y, as.vector(RA %*% betabar))
  IR = backsolve(chol(crossprod(W)), diag(nvar))
  btilde = crossprod(t(IR)) %*% crossprod(W, z)
  res = z - W %*% btilde
  s = t(res) %*% res
  sigmasq=tao
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
# Import original data inferred from vpc_real
datamraw=read.csv("Delivery_Match_VPC.csv", header=T)

databraw=read.csv("QS_POST.csv", header=T)

vpcraw=read.csv("vpc_ratio=0_full.csv", header=T)

datam=datamraw[datamraw$thedate<=12, ]
datab=databraw[databraw$thedate<=12, ]

# we add "vpc" and "flag" into datab
datab=cbind(datab, vpcraw$vpc, vpcraw$flag)
colnames(datab)[9:10]=c("vpc", "flag")

# we first fix abnormalies in vpc
datab$ratio=datab$vpc/datab$bid
# indexabn=which(datab$flag>1 | datab$ratio>5 | datab$ratio<.2)
indexabn=which(datab$flag>1)

# percentage of abnormal
length(indexabn)/nrow(datab)
avgratio=mean(datab$vpc[-indexabn]/datab$bid[-indexabn])

# fix those vpc with flag>2 using avgratio between vpc and bid
datab$vpc[indexabn]=datab$bid[indexabn]*avgratio
summary(datab$bid/datab$vpc)


########################################################################################
load("Supply 1 data.Rdata")

# transfer lists into matricies
S_list=matrix(unlist(S_list), ncol = 50, byrow = FALSE)
# take the mean of S_list
# may not be useful in subsequent analysis
S_vec=apply(S_list, 1, mean)

# take the mean of beta_u_mat & beta_s_mat
beta_u_vec=apply(beta_u_mat, 2, mean)
beta_s_vec=apply(beta_s_mat, 2, mean)

# We calculate iscore by subtracting the impact of n_ad from S_vec
iscore=S_vec-beta_s_vec[9]*(log(datam$n_ad)-mean(log(datam$n_ad)))
iscore=exp(iscore)
summary(iscore)

datam=cbind(datam, iscore)

########################################################################################
# Next we need to compute the average iscore for each datesubcat, also import corrresponding npotad, cat dummies, weekend
subcatmat=aggregate(x=cbind(datam[, c("iscore", "n_ad", "weekend")], datam[, 22:34]),
                    by=list(date_subcat_id=datam$date_subcat_id), FUN="mean")

summary(subcatmat)

# merge subcatmat with databraw
datab2=merge(datab, subcatmat, by.x="date_subcat_id", by.y="date_subcat_id", all.x=TRUE)
datab2=as.data.frame(datab2)
# update ratio between bid and vpc
datab2$ratio=datab2$bid/datab2$vpc
summary(datab2$ratio)

# export datab2
write.csv(datab2, "VPC_Reg_ratio=0.csv",row.names=FALSE)


########################################################################################
# next we do some plots
########################################################################################
library(bayesm)
library(msm)
library(LearnBayes)
library(ggplot2)
########################################################################################
rm(list = ls())
# import data
datab2 = read.csv("VPC_Reg_ratio=0.csv", header = T)

# we first plot the histogram of vpc and ratio
par(mfrow=c(1,2))

hist(datab2$vpc, xlab="(a) Value-Per-Click", col="6", xlim=c(0,20), probability =TRUE, breaks=200, right=FALSE, main=NULL)
axis(1, at=seq(0, 20, by=1) ) 

hist(datab2$ratio, xlab="(b) Ratio=Bid/Value-Per-Click", col="6", xlim=c(0.2,1.5), probability =TRUE, breaks=200, right=FALSE, main=NULL)
axis(1, at=seq(0.2, 1.5, by=.1) ) 

# next we do scatterplots
index=which(datab2$bid<50 & datab2$vpc<60)
p1=ggplot(datab2[index, ], aes(x=vpc, y=bid))+geom_point(alpha=0.3)+geom_smooth(color="red", fill="blue")+scale_size_area()+
  scale_x_continuous("(a) Value-per-click") + 
  scale_y_continuous("Bid")
print(p1)

png(filename="bid_vpc.png", 
    width = 1200, height = 1200, res=200) 
print(p1)
dev.off()

p2=ggplot(datab2, aes(x=qscore, y=bid))+geom_point(alpha=0.5)+geom_smooth(color="red", fill="blue")+
  scale_x_continuous("(b) Quality Score") + 
  scale_y_continuous("Bid")
print(p2)

png(filename="bid_qs.png", 
    width = 1200, height = 1200, res=200) # perhaps width/height as well
print(p2)
dev.off()



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(p1, p2, cols=2)



########################################################################################
########################################################################################
# Now we are ready for estimation
########################################################################################
# Set parameters
ad_id=datab2$adid2
N_ad=max(ad_id)
lvpc=log(datab2$vpc)

# We create a adlist for datab
# also create a ad_exist as a N_ad by 1 vec, if ad_id==i appears in data, then ad_exist has an element i
# ad_exist2 recoreds all ads that have at least two observations in data
adlist=list()
ad_exist=rep(0, N_ad)
ad_exist2=rep(0, N_ad)

for (i in 1:N_ad){
  index=which(ad_id==i)
  adlist[[i]]=index
  if(length(index)>0){
    ad_exist[i]=i
  }
  if(length(index)>=2){
    ad_exist2[i]=i
  }
}
temp=sapply(adlist, length)
length(which(temp>=1))/N_ad
length(which(temp>=2))/N_ad 

ad_exist=ad_exist[which(ad_exist>0)]
ad_exist2=ad_exist2[which(ad_exist2>0)]


# exclude quality score
X=as.matrix(cbind(1, datab2$rprice, log(datab2$iscore), log(datab2$n_ad), datab2$weekend,
                  datab2[, 15:18]))
colnames(X)=c("intercept","rprice", "lniscore", "lnnpotad", "weekend",
              "rootcat2", "rootcat3","rootcat4","rootcat5")
nvar=ncol(X)

# mean-center
for (i in 2:nvar){
  X[, i]=X[, i]-mean(X[, i])
}

# define initial values for parameters
beta=rep(0, nvar)
tao=0.1
var_k=0.1
eta_k=rnorm(N_ad, 0, var_k^0.5)
etavec=eta_k[ad_id]


ndraw=20000
# betadraw is a ndraw by ? matrix which stores beta, tao in each row
betadraw=matrix(0, nrow=ndraw, ncol=nvar)
vardraw=matrix(0, nrow=ndraw, ncol=2)
# we also record every 100 draws of random coefficients
eta_list=list()
# Like stores the log-likelihood of lvpc
Like=list()
########################################################################################
# start Gibbs sampler

i=1

while (i<=ndraw){
  
  #########################################################################
  # update beta 
  y=lvpc-etavec
  out=update_beta(y, X)
  beta=as.vector(out$beta)
  tao=as.numeric(out$var)
  
  # update eta_k & var_k
  y=lvpc-X%*%beta
  j=1
  while(j<=N_ad){
    index=adlist[[j]]
    if(length(index)>=2){
      ytemp=y[index]
      xtemp=rep(1, length(ytemp))
      # posterior mean
      pmean=(var_k*crossprod(ytemp,xtemp))/(var_k*crossprod(xtemp)+tao)
      # posterior var
      pvar=var_k*tao/(var_k*crossprod(xtemp)+tao)
      eta_k[j]=rnorm(1, pmean, pvar^0.5)  
    }else{
      eta_k[j]=0
    }
    j=j+1
  }
  etavec=eta_k[ad_id]
  var_k=as.numeric(update_tao(eta_k[ad_exist2]))
  
  
  # record updated draws
  betadraw[i, ]=beta
  vardraw[i, ]=c(tao, var_k)
  
  if (i %% 100 ==0){
    cat("iter:", i, fill=TRUE)
    cat("beta:", beta, fill=TRUE)
    cat("tao, var_k:", c(tao, var_k), fill=TRUE)
  
    eta_list[[i/100]]=eta_k
    Like[[i/100]]=-sum((lvpc-X%*%beta-etavec)^2/(2*tao))-length(lvpc)/2*log(tao)
    cat("loglikelihood:", Like[[i/100]], fill=TRUE)
    cat('\n')
  }
  i=i+1
}

apply(betadraw[seq(10000,ndraw,20),],2,quantile,probs=c(.025,.05,.50,.95,.975))

acf(betadraw[seq(10000,ndraw,20), 1])

hist(betadraw[seq(5000,ndraw,20), 1])


apply(sqrt(vardraw[seq(10000,ndraw,20),]),2,quantile,probs=c(.025,.05,.50,.95,.975))




table8result1 = cbind( apply(betadraw[seq(10000,ndraw,20),],2, mean) ,    apply(betadraw[seq(10000,ndraw,20),],2, sd) )
table8result2 = cbind( apply( sqrt(vardraw[seq(10000,ndraw,20),]),2, mean) ,    apply( sqrt(vardraw[seq(10000,ndraw,20),]) ,2, sd) )
table8result = rbind(table8result1[1:5, ] , table8result2)
rownames(table8result) = c("intercept","rprice", "lniscore", "lnnpotad", "weekend" , "err1 sdvkd" , "err2 sdvk"  ) 

