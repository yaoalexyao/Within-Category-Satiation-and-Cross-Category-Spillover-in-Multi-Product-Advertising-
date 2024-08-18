########################################################################################
# This file aims to infer advertiser's (daily) value-per-click based on the FOC of their 
# profit-maximization problems.
# We use the estimation results from Supply 1 code to simulate the profit function
# Please note: 
# we reconstruct the profit function by including an integration of the stochastic term in weighted bid
#    in other words, the functions of ECPC, EPROFIT, PROBPOS are all different from before





########################################################################################
library(bayesm)
library(msm)
library(LearnBayes)
library(MCMCpack)
library(numDeriv)
library(fastGHQuad)
library(doParallel)
library(apollo)
########################################################################################
rm(list = ls())
########################################################################################
# function of drawing a N:N (highest) order statistics from F(x)~N(a, var)
OSTOPDRAW=function(N, a, var, u){
  # Arguments
  # N, a, var are explained
  # u: a draw from Unif(0, 1)
  
  # first compute u^(1/N)
  x=exp(log(u)/N)
  out=qnorm(x, mean=a, sd=var^0.5)
  return(out)
}

########################################################################################
# function of drawing the (N-M)th order statistics from N conditional on the (N-M+1)th equals x
# M can be regarded as the position
# where F~N(a, var)
OSDRAWBELOW=function(x, N, M, a, var, u){
  y=exp(log(u)/(N-M))
  y=y*pnorm(x, a, var^0.5)
  out=qnorm(y, a, var^0.5)
  return(out)
}

########################################################################################
# function of drawing the (N-M+2)th order statistics from N conditional on the (N-M+1)th equals x
# M can be regarded as the position
# where F~N(a, var)
OSDRAWABOVE=function(x, N, M, a, var, u){
  y=1-exp(log(1-u)/(M-1))
  z=pnorm(x, a, var^0.5)
  y=(1-z)*y+z
  out=qnorm(y, a, var^0.5)
  return(out)
}

########################################################################################
# function of calculation expected CPC given the no. of ads for display in each cat, 
# the position of focal ad in its cat and the distribution of others' weighted bid
ECPC=function(n_disp, n_pot, i, pos, b, QS, a, var, u){
  # Arguments
  # n_disp: a N by 1 vector, each element records the number of ads for display for each potential cat,
  #         where N is N_cat_pot, and sum(n_disp)=8 because the ad quota is 8.
  # n_pot: a N by 1 vector, each element records the number of competing ads for that category.
  # i: the ith category in n_disp is for the focal advertiser
  # pos: pos is the rank of focal advertiser in its category. pos>=1 and <=n_disp[i]
  # b: bid submitted
  # QS: quality score
  # a, var: the distribution of others weighted bid follow ln(wb)~N(a, var)
  # u: a 8 by 1 vector drawn from unif(8); u is used to create weighted bid for all other 7 advertisers.
  
  # first check if the cat of focal advertiser has positive quota
  if(n_disp[i]==0){
    CPC=0
  }else{
    # we start the draw for the cat of focal advertiser
    lwb=log(b)+log(QS)
    # we draw the lwb ranked below focal advertiser given the position
    lwb_b=OSDRAWBELOW(lwb, n_pot[i], pos, a, var, u[1])
    
    # next we compute CPC for the focal advertiser
    CPC=exp(lwb_b)/QS
  }
  return(CPC)
}

########################################################################################
# function of calculating the probability of getting a position for the focal advertiser in its category
PROBPOS=function(eps, pos, n_pot, b, QS, a, var, tao){
  # Arguments
  # eps: an additive error term in log(weighted bid)
  # pos: position
  # b: bid
  # QS: quality score
  # the distribution of others' ln(wb)~N(a, var)
  # the stochastic error in ranking model is error~N(0, tao)
  # next we calculate the probability of getting a fixed position when bidding b
  lwb=log(b)+log(QS)+eps*tao^0.5
  
  out=rep(0, length(eps))
  for (i in 1:length(eps)){
    # the CDF for lwb is
    cdf=pnorm(lwb[i], mean=a, sd=(var+tao)^0.5)
    out[i]=choose(n_pot-1, pos-1)*(1-cdf)^(pos-1)*(cdf)^(n_pot-pos)
  }
  return(out)
}

########################################################################################
# function of calculating profit from clicks given the position and relevant parameters for that impression
PROFITPOS1=function(eps, n_disp, n_pot, i, pos, b, QS, v, ratio, a, var, tao, u){
  # Arguments
  # eps: an additive error term in log(weighted bid)
  # n_disp: a N by 1 vector, each element records the number of ads for display for each potential cat,
  #         where N is N_cat_pot, and sum(n_disp)=8 because the ad quota is 8.
  # n_pot: a N by 1 vector, each element records the number of competing ads for that category.
  # i: the ith category in n_disp is for the focal advertiser
  # pos: pos is the rank of focal advertiser in its category. pos>=1 and <=n_disp[i]
  # b: bid submitted
  # QS: quality score
  # v: value-per-click
  # ratio: vpm/vpc, vpm is value per thousand impressions
  # a, var, tao: the distribution of others weighted bid follow ln(wb)~N(a, var+tao)
  # u: a 8 by 1 vector drawn from unif(8); u is used to create weighted bid for all other 7 advertisers.
  
  # profit=prob of being displayed * (v-cpc)
  CPC=ECPC(n_disp, n_pot, i, pos, b, QS, a, var, u)
  prob=PROBPOS(eps, pos, n_pot[i], b, QS, a, var, tao)
  out=prob*(v-CPC)  
  return(out)
}


########################################################################################
# function of calculating profit from impressions given the position and relevant parameters for that impression
# this function is no longer useful
PROFITPOS2=function(eps, n_disp, n_pot, i, pos, b, QS, v, ratio, a, var, tao, u){
  # Arguments
  # eps: an addtive error term in log(weighted bid)
  # n_disp: a N by 1 vector, each element records the number of ads for display for each potential cat,
  #         where N is N_cat_pot, and sum(n_disp)=8 because the ad quota is 8.
  # n_pot: a N by 1 vector, each element records the number of competing ads for that category.
  # i: the ith category in n_disp is for the focal advertiser
  # pos: pos is the rank of focal advertiser in its category. pos>=1 and <=n_disp[i]
  # b: bid submitted
  # QS: quality score
  # v: value-per-click
  # ratio: vpm/vpc, vpm is value per thousand impressions
  # a, var, tao: the distribution of others weighted bid follow ln(wb)~N(a, var+tao)
  # u: a 8 by 1 vector drawn from unif(8); u is used to create weighted bid for all other 7 advertisers.
  
  # profit=prob of being displayed * v * ratio/1000
  prob=PROBPOS(eps, pos, n_pot[i], b, QS, a, var, tao)
  out=prob*v*ratio/1000 
  return(out)
}

########################################################################################
# function of integrating out eps in function "PROFITPOS1"
EPROFITPOS1=function(n_disp, n_pot, i, pos, b, QS, v, ratio, a, var, tao, u, rule){
  out=ghQuad(PROFITPOS1, rule, n_disp, n_pot, i, pos, b, QS, v, ratio, a, var, tao, u)
  return(out)
}

# function of integrating out eps in function "PROFITPOS2"
EPROFITPOS2=function(n_disp, n_pot, i, pos, b, QS, v, ratio, a, var, tao, u, rule){
  out=ghQuad(PROFITPOS2, rule, n_disp, n_pot, i, pos, b, QS, v, ratio, a, var, tao, u)
  return(out)
}

rule10 <- gaussHermiteData(10)
rule20 <- gaussHermiteData(20)
rule30 <- gaussHermiteData(30)
rule50 <- gaussHermiteData(50)
rule100 <- gaussHermiteData(100)

########################################################################################
# function of assigning ad quota to each category based on log(interest score) S
ADASSIGN=function(S, n_total){
  # S: log of interest score for a set of potential categories
  # n_total: 8 in our case
  
  # first compute the probability p vector
  p=exp(S)/sum(exp(S))
  # record the floor of p*n_total
  out1=floor(p*n_total)
  out2=p*n_total-out1
  # then we assign the rest of quota (n_total-sum(out1)) based on the ranking of out2
  n_remain=n_total-sum(out1)
  index=rank(-out2)
  for(i in 1:n_remain){
    temp=which(index==i)
    out1[temp]=out1[temp]+1
  }
  return(out1) 
}

########################################################################################
# function of calculating expected daily profit for the focal advertiser given multiple draws of S and pclick
# we change the EPROFIT from previous version by changing the input u
# u is a matrix
# !!the new parameter ratio stands for ratio between value per thousand impressions and one click
EPROFIT=function(n_ad_disp, n_pot_cat, n_pot_ad, cat_index, pconsider, utility, X_std, beta_u,
                 b, QS, v, ratio, a, var, tao, u, rule){
  # Arguments
  # n_ad_disp: a N by 1 vector which records the assigned ad quota for each potential category for each consumer during a day,
  #            N=sum(n_pot_cat)
  # n_pot_cat: a n_user by 1 vector records the number of potential categories for each consumer. its length is shorter than S
  # n_pot_ad: a vector records the number of competing ads in each potential category.it has the same length as S
  # cat_index: a n_user by 1 vector records the position of the focal cat among all potential cat for each consumer. it has the same length as n_pot_cat
  # pconsider: a n_user by 1 vector records the user's probability of considering clicking for a particular impression
  # utility: a n_user by 1 vector records the consumers' utility of clicking on the focal ad (except for the effect of n_ad, n_related)
  # X_std: a n_user by 2 matrix, each row records the standardized n_ad, n_related
  # beta_u: a n_user by 2 matrix, each row records the user's coefficient for n_ad, n_related
  
  # b: bid
  # QS: quality score
  # v: value-per-click
  # ratio: vpm/vpc
  # a, var: the log(weighted bid) of others follow N(a, var)
  # tao: variance of error term in ranking model
  # u: a N by 8 matrix, where each element is drawn from unif(0,1)
  
  # define n_user
  n_user=length(n_pot_cat)
  
  # define profit as a n_user by 1 vector
  profit=rep(0, n_user)
  i=1
  j=1
  while(i<=n_user){
    # pre-calculated assigned quota for each potential cat
    n_disp=n_ad_disp[j:(j+n_pot_cat[i]-1)]
    npotadtemp=n_pot_ad[j:(j+n_pot_cat[i]-1)]
    indextemp=cat_index[i]
    # We need to make sure that n_pot_ad>=n_disp, otherwise we cannot simulate CPC
    npotadtemp=mapply(max, npotadtemp, n_disp)
    # define the vector of random draws
    uvec=u[i, ]
    
    n_ad=n_disp[indextemp]
    if(n_ad==0){ # if no position for the cat of focal ad
      profit[i]=0
    }else{
      # create a n_disp[indextemp] vector recording the EProfit of position 1,2... 
      profitemp1=rep(0, n_ad) # save profit from clicks
      profitemp2=rep(0, n_ad) # save profit from impressions
      for(pos in 1:n_ad){
        profitemp1[pos]=EPROFITPOS1(n_disp, npotadtemp, indextemp, pos, b, QS, v, ratio, a, var, tao, uvec, rule)
        profitemp2[pos]=EPROFITPOS2(n_disp, npotadtemp, indextemp, pos, b, QS, v, ratio, a, var, tao, uvec, rule)
      }
      # compute expected click
      pclick=pnorm(utility[i]+crossprod(X_std[i, ], beta_u))*pconsider[i]
                     
      profit[i]=pclick*sum(profitemp1)+sum(profitemp2)
    }
    j=j+n_pot_cat[i]
    i=i+1
  }
  out=sum(profit)
  
  return(out) 
}

########################################################################################
# function of calculating expected daily profit for the focal advertiser given multiple draws of S and pclick
# we change the EPROFIT from previous version by changing the input u
# u is a matrix
# compared to EPROFIT, the changes include
#  ratio=0 (no longer capture value from impressions)
#  expected clicks matrix provided, denoted by pclick
EPROFIT_V2=function(n_ad_disp, n_pot_cat, n_pot_ad, cat_index, pclick,
                 b, QS, v, a, var, tao, u, rule){
  # Arguments
  # n_ad_disp: a N by 1 vector which records the assigned ad quota for each potential category for each consumer during a day,
  #            N=sum(n_pot_cat)
  # n_pot_cat: a n_user by 1 vector records the number of potential categories for each consumer. its length is shorter than S
  # n_pot_ad: a vector records the number of competing ads in each potential category.it has the same length as S
  # cat_index: a n_user by 1 vector records the position of the focal cat among all potential cat for each consumer. it has the same length as n_pot_cat
  # pclick: a n_user by 1 vector records the user's expected number of clicks on the focal ad
  
  # b: bid
  # QS: quality score
  # v: value-per-click
  # a, var: the log(weighted bid) of others follow N(a, var)
  # tao: variance of error term in ranking model
  # u: a N by 8 matrix, where each element is drawn from unif(0,1)
  
  # define n_user
  n_user=length(n_pot_cat)
  
  # define profit as a n_user by 1 vector
  profit=rep(0, n_user)
  i=1
  j=1
  while(i<=n_user){
    # pre-calculated assigned quota for each potential cat
    n_disp=n_ad_disp[j:(j+n_pot_cat[i]-1)]
    npotadtemp=n_pot_ad[j:(j+n_pot_cat[i]-1)]
    indextemp=cat_index[i]
    # We need to make sure that n_pot_ad>=n_disp, otherwise we cannot simulate CPC
    npotadtemp=mapply(max, npotadtemp, n_disp)
    # define the vector of random draws
    uvec=u[i, ]
    
    n_ad=n_disp[indextemp]
    if(n_ad==0){ # if no position for the cat of focal ad
      profit[i]=0
    }else{
      # create a n_disp[indextemp] vector recording the EProfit of position 1,2... 
      profitemp1=rep(0, n_ad) # save profit from clicks
      for(pos in 1:n_ad){
        profitemp1[pos]=EPROFITPOS1(n_disp, npotadtemp, indextemp, pos, b, QS, v, ratio=0, a, var, tao, uvec, rule)
      }
      profit[i]=pclick[i]*sum(profitemp1)
    }
    j=j+n_pot_cat[i]
    i=i+1
  }
  out=sum(profit)
  
  return(out) 
}


########################################################################################
# function of calculating FOC of the expected profit function w.r.t bid
FOC=function(v, b, n_ad_disp, n_pot_cat, n_pot_ad, cat_index, pconsider, utility, X_std, beta_u, 
             QS, ratio, a, var, tao, u, rule){
  out=grad(EPROFIT, b, method="simple", v=v, n_ad_disp=n_ad_disp, n_pot_cat=n_pot_cat, n_pot_ad=n_pot_ad, cat_index=cat_index, 
           pconsider=pconsider, utility=utility, X_std=X_std, beta_u=beta_u,
           QS=QS, ratio=ratio, a=a, var=var, tao=tao, u=u, rule=rule)
  return(out)
}

########################################################################################
# function of calculating FOC of the expected profit function w.r.t bid
FOC_V2=function(v, b, n_ad_disp, n_pot_cat, n_pot_ad, cat_index, pclick,
             QS, a, var, tao, u, rule){
  out=grad(EPROFIT_V2, b, method="simple", v=v, n_ad_disp=n_ad_disp, n_pot_cat=n_pot_cat, n_pot_ad=n_pot_ad, cat_index=cat_index, 
           pclick=pclick, QS=QS, a=a, var=var, tao=tao, u=u, rule=rule)
  return(out)
}



########################################################################################
# Import data
datam=read.csv("Delivery_Match_VPC.csv", header=T)

datab=read.csv("QS_POST.csv", header=T)

load("Supply 1 data.Rdata")

# transfer lists into matricies
S_list=matrix(unlist(S_list), ncol = 50, byrow = FALSE)
# take the mean of S_list
S_vec=apply(S_list, 1, mean)

# take the mean of beta_u_mat
beta_u_vec=apply(beta_u_mat, 2, mean)

eps_u_list=matrix(unlist(eps_u_list), ncol = 50, byrow = FALSE)
eps_u_vec=apply(eps_u_list, 1, mean)

eta_u_list=matrix(unlist(eta_u_list), ncol = 50, byrow = FALSE)
eta_u_vec=apply(eta_u_list, 1, mean)


########################################################################################
# Set parameters
T=12
# note we only use first 12 days in model estimation
datam=datam[datam$thedate<=T, ]
datab=datab[datab$thedate<=T, ]

N_user=length(unique(datam$user_id3))
N_ad=8 # no of ads for display
N_cat_all=max(datam$subcat_id3)
N_userroot=max(datam$usroot_id) # here the userroot refers to <user, rootcat>
N_ad_all=max(datab$adid2)
N_datecat_all=max(datab$date_subcat_id)
N_pv=max(datam$pv_id2) # no of pvs

# some useful index vectors
pv_id=datam$pv_id2
cat_id_s=datam$subcat_id3
userroot_id_s=datam$usroot_id
date_cat_id_s=datam$date_subcat_id
ad_id_b=datab$adid2
date_cat_id_b=datab$date_subcat_id

# define the distribution of log(weighted bid)~N(a, var) and the tao, the variance of error term in ranking model
a=mean(log(datab$bid)+log(datab$qs_post))
var=var(log(datab$bid)+log(datab$qs_post))
tao=1.208


########################################################################################
# create a matrix of random draws from uniform distribution for each <date, adid> in datab
set.seed(123)
ranmat=matrix(runif(N_pv*N_ad), N_pv, N_ad)
rm(.Random.seed)


########################################################################################
# transfer S_vec into a X by 1 vector, where each element is the corresponding n_ad
# generate a vec Sim_vec, N_pv by 1, each element is the corresponding n_related2 for an impression
n_ad_vec=rep(0, nrow(datam))
Sim_vec=rep(0, N_pv)
i=1
j=1
set.seed(123)

while(i<=nrow(datam)){
  n_pot=datam$n_pot_subcat[i]
  S_temp=S_vec[i:(i+n_pot-1)]
  # transfer S_temp into n_ad
  n_ad_temp=ADASSIGN(S_temp, 8)
  n_ad_vec[i:(i+n_pot-1)]=n_ad_temp
  
  # next compute similarity
  # find those categories where the user has visited
  RCC_temp=datam$rcc_click[i:(i+n_pot-1)]
  
  # create a n_pot by 1 vector where a cat that has been clicked
  click=rep(0, n_pot)
  click[RCC_temp<=2]=1
  
  Sim_vec[j]=sum(n_ad_temp*click)
  
  i=i+n_pot
  j=j+1
  if(j%%100==0){
    cat(j, fill=TRUE)
  }
}
rm(.Random.seed)



# ########################################################################################
# # we use demand model code to predict expected number of clicks given datam and predicted n_ad_vec under S_vec
# ########################################################################################
# 

database1 = read.csv("cleaned data.csv",header=TRUE)


i=1

while(i<=nrow(database1)){
  pv_id_temp = database1$pv_id2[i]
  # find all index where datam$pv_id == pv_id_temp
  index_temp = which(datam$pv_id2 == pv_id_temp)
  if(length(index_temp)==0){
    i=i+1
  }else{
    # create a 8 by 1 vector recording the subcat_id3 with positive n_ad for this pv
    index_temp2 = which(n_ad_vec[index_temp]>0)
    n_ad_temp = n_ad_vec[index_temp][index_temp2]

    # update subcat_id
    database1[i, 16:23] = rep(datam$subcat_id3[index_temp[index_temp2]], times = n_ad_temp)

    # update npv, click, buy
    index_temp3 = rep(index_temp[index_temp2], times = n_ad_temp)
    database1[i, 32:39] = log(datam$npv14[index_temp3]+1)
    database1[i, 41:48] = log(datam$click14[index_temp3]+1)
    database1[i, 50:57] = datam$buy14[index_temp3]

    # update rcc_pv, rcc_click, rcc_buy
    database1[i, 59:66] = log(datam$rcc_pv[index_temp3]+1)
    database1[i, 68:75] = log(datam$rcc_click[index_temp3]+1)
    database1[i, 77:84] = log(datam$rcc_buy[index_temp3]+1)

    # update availnew1 to availnew9
    temp = database1[i, 17:23] - database1[i, 16:22]
    temp = c(as.numeric(database1[i, 17]), as.vector(temp))
    temp2 = rep(0, 8)
    temp2[temp>0] =1
    database1[i, 86:93] = temp2

    if(i%%500==0){
      cat(i, fill=TRUE)
    }
    i=i+1
  }
}

database = database1



### Initialise code
apollo_initialise()

apollo_control = list(
  modelName  ="Demand eMDC model with latent class",
  modelDescr ="Demand eMDC model with latent class",
  indivID    ="user_id3",
  nCores     = 7,
  workInLogs=F
)





#####This function calculates the ad click probability using subcategory information, predicted clicks from MDCEV and choice probability
add_restrict_imputevalue = function(x1) {
  x = unlist(x1[1:8]) 
  y = x1[9:16] 
  z = x1[17:24] 
  
  for (i in 1:length(unique(x))) {
    temp = max(y[x == unique(x)[i]]) # the highest clicks in the subcategory i 
    y[x == unique(x)[i]] = min(temp, sum( x == unique(x)[i] ) )
    z[x == unique(x)[i] ] = z[x == unique(x)[i] ] / sum(z[x == unique(x)[i] ]) * y[x == unique(x)[i]]
  }
  z[z>1] = 1
  return(z)
}


apollo_beta1 = c(
  gamma_v1      = 1.301, 
  gamma_v2      = 0.600, 
  gamma_v3      = .501, 
  gamma_v4      = .477, 
  gamma_v5      = .664, 
  delta_v1      = -3.967,
  delta_v2         = -3.849,
  delta_v3       = -3.897,
  delta_v4     = -3.357,
  delta_v5 = -3.348,
  delta_outside      = 0,
  
  
  

  bweekend = 0.088, bnpv = 0.026,  bnclick = 0.109,  bnbuy = 	-0.080  , nrccpv=0.010, nrccclick=0.348, nrccbuy=0.025, 
  badstockpsai = 0.637 ,  bresidual = 0.139 , bvariety = 0.043,
  # Compl/subst
  d12= 0.005, d13= 0.014, 
  d14= 0.022, d15= 0.003, 
  d23= 0.015, d24= 0.017, 
  d25= 0.014, d34= 0.023, 
  d35= 0.042, d45= -0.002
  
  
  
)

apollo_beta2 = c(
  
  gamma_v1      = 0.903	, 
  gamma_v2      = 0.557, 
  gamma_v3      = 0.489, 
  gamma_v4      = 0.397, 
  gamma_v5      = 0.197, 
  delta_v1      = -2.536,
  delta_v2         = -2.465,
  delta_v3       = -2.394,
  delta_v4     = -2.529,
  delta_v5 = -2.412,
  delta_outside      = 0,
  
  bweekend = -0.029	, bnpv = 0.024,  bnclick = 0.062,  bnbuy = -0.727 , nrccpv=0.037, nrccclick=0.196, nrccbuy=	-0.476, 
  badstockpsai = 0.305	 ,  bresidual = 0.137 , bvariety = 0.026	,
  d12= 0.009, d13= 0.001, 
  d14= 0.006, d15= -0.006, 
  d23= 0.005, d24= 0.008, 
  d25= -0.007, d34= 0.023, 
  d35= 0.013, d45= -0.008
)

segmentco = c(6.001,-4.084, -2.799)

membership0 = read.csv("membership.csv",header=TRUE)
membership = membership0[,2]

database = database0

apollo_fixed = c( "delta_outside"    ) 



t1<-Sys.time()
t1


#######################Segment 1  
apollo_beta = apollo_beta1

apollo_probabilities1=function(apollo_beta = apollo_beta1, apollo_inputs, functionality="estimate"){
  
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### Define individual alternatives
  alternatives = c("v1", 
                   "v2", 
                   "v3", 
                   "v4", 
                   "v5", 
                   "v6", 
                   "v7", 
                   "v8" 
  )
  
  ### Define availabilities
  avail = list(v1  = availnew1,   
               v2     = availnew2,
               v3   = availnew3,
               v4 = availnew4,
               v5 = availnew5,
               v6   = availnew6,
               v7 = availnew7,
               v8 = availnew8 
  )
  
  ### Define continuous consumption for individual alternatives
  continuousChoice = list(v1  =clicknew1 * availnew1,
                          v2     =clicknew2 * availnew2,
                          v3   =clicknew3 * availnew3,
                          v4 =clicknew4 * availnew4,
                          v5 =clicknew5 * availnew5,
                          v6   =clicknew6 * availnew6,
                          v7 =clicknew7 * availnew7,
                          v8 =clicknew8 * availnew8#,
  )
  
  
  ### Define alpha parameters
  alpha = list(v1  = 1e-3 , 
               v2     = 1e-3 , 
               v3   = 1e-3 , 
               v4 = 1e-3 , 
               v5 = 1e-3 ,
               v6   = 1e-3 ,
               v7  = 1e-3 , 
               v8 = 1e-3 #, 
  )
  
  
  ### Define costs for individual alternatives
  cost = list(v1      = 1, 
              v2         = 1,
              v3       = 1,
              v4     = 1,
              v5 = 1,
              v6       = 1, 
              v7      = 1,
              v8     = 1 #,
  )
  
  ### Define budget
  budget = clicktotal
  
  emdc_settings <- list(continuousChoice = continuousChoice, 
                        avail            = avail,
                        budget           = budget,
                        sigma            = 0.99, 
                        cost             = cost)
  
  
  ### ### Compute class-specific utilities
  V = list()
  
  V[["v1"    ]] =   delta_v1*(rootcat1==1) +  delta_v2*(rootcat1==2) +  delta_v3*(rootcat1==3) +  
    delta_v4*(rootcat1==4) +  delta_v5*(rootcat1==5) + 
    bweekend*weekend + bnpv*npvnew1 + bnclick*nclicknew1 + bnbuy*nbuynew1 + 
    nrccpv*rccpvnew1 + nrccclick*rccclicknew1 + nrccbuy*rccbuynew1 + badstockpsai*ad11 + bresidual*res1 + bvariety*Lag_Distinct
  
  V[["v2"    ]] = delta_v1*(rootcat2==1) +  delta_v2*(rootcat2==2) +  delta_v3*(rootcat2==3) +  
    delta_v4*(rootcat2==4) +  delta_v5*(rootcat2==5) + 
    bweekend*weekend + bnpv*npvnew2 + bnclick*nclicknew2 + bnbuy*nbuynew2 + 
    nrccpv*rccpvnew2 + nrccclick*rccclicknew2 + nrccbuy*rccbuynew2 + badstockpsai*ad21 + bresidual*res2 + bvariety*Lag_Distinct
  
  V[["v3"  ]] = delta_v1*(rootcat3==1) +  delta_v2*(rootcat3==2) +  delta_v3*(rootcat3==3) +  
    delta_v4*(rootcat3==4) +  delta_v5*(rootcat3==5) + 
    bweekend*weekend + bnpv*npvnew3 + bnclick*nclicknew3 + bnbuy*nbuynew3 + 
    nrccpv*rccpvnew3 + nrccclick*rccclicknew3 + nrccbuy*rccbuynew3 + badstockpsai*ad31 + bresidual*res3 + bvariety*Lag_Distinct
  
  V[["v4"]] = delta_v1*(rootcat4==1) +  delta_v2*(rootcat4==2) +  delta_v3*(rootcat4==3) +  
    delta_v4*(rootcat4==4) +  delta_v5*(rootcat4==5) + 
    bweekend*weekend + bnpv*npvnew4 + bnclick*nclicknew4 + bnbuy*nbuynew4 + 
    nrccpv*rccpvnew4 + nrccclick*rccclicknew4 + nrccbuy*rccbuynew4 + badstockpsai*ad41 + bresidual*res4 + bvariety*Lag_Distinct
  
  V[["v5"]] = delta_v1*(rootcat5==1) +  delta_v2*(rootcat5==2) +  delta_v3*(rootcat5==3) +  
    delta_v4*(rootcat5==4) +  delta_v5*(rootcat5==5) + 
    bweekend*weekend + bnpv*npvnew5 + bnclick*nclicknew5 + bnbuy*nbuynew5 + 
    nrccpv*rccpvnew5 + nrccclick*rccclicknew5 + nrccbuy*rccbuynew5 + badstockpsai*ad51 + bresidual*res5 + bvariety*Lag_Distinct
  
  V[["v6"  ]] = delta_v1*(rootcat6==1) +  delta_v2*(rootcat6==2) +  delta_v3*(rootcat6==3) +  
    delta_v4*(rootcat6==4) +  delta_v5*(rootcat6==5) + 
    bweekend*weekend + bnpv*npvnew6 + bnclick*nclicknew6 + bnbuy*nbuynew6 + 
    nrccpv*rccpvnew6 + nrccclick*rccclicknew6 + nrccbuy*rccbuynew6 + badstockpsai*ad61 + bresidual*res6 + bvariety*Lag_Distinct
  
  V[["v7"]] = delta_v1*(rootcat7==1) +  delta_v2*(rootcat7==2) +  delta_v3*(rootcat7==3) +  
    delta_v4*(rootcat7==4) +  delta_v5*(rootcat7==5) + 
    bweekend*weekend + bnpv*npvnew7 + bnclick*nclicknew7 + bnbuy*nbuynew7 + 
    nrccpv*rccpvnew7 + nrccclick*rccclicknew7 + nrccbuy*rccbuynew7 + badstockpsai*ad71 + bresidual*res7 + bvariety*Lag_Distinct
  
  V[["v8"]] = delta_v1*(rootcat8==1) +  delta_v2*(rootcat8==2) +  delta_v3*(rootcat8==3) +  
    delta_v4*(rootcat8==4) +  delta_v5*(rootcat8==5) + 
    bweekend*weekend + bnpv*npvnew8 + bnclick*nclicknew8 + bnbuy*nbuynew8 + 
    nrccpv*rccpvnew8 + nrccclick*rccclicknew8 + nrccbuy*rccbuynew8 + badstockpsai*ad81 + bresidual*res8 + bvariety*Lag_Distinct

  
  
  ### Define gamma parameters
  gamma = list(v1      = exp( gamma_v1*(rootcat1==1) +  gamma_v2*(rootcat1==2) +  gamma_v3*(rootcat1==3) +  
                                gamma_v4*(rootcat1==4) +  gamma_v5*(rootcat1==5) ),
               
               v2      = exp( gamma_v1*(rootcat2==1) +  gamma_v2*(rootcat2==2) +  gamma_v3*(rootcat2==3) +  
                                gamma_v4*(rootcat2==4) +  gamma_v5*(rootcat2==5) ),
               
               
               v3      = exp( gamma_v1*(rootcat3==1) +  gamma_v2*(rootcat3==2) +  gamma_v3*(rootcat3==3) +  
                                gamma_v4*(rootcat3==4) +  gamma_v5*(rootcat3==5) ),
               
               v4      = exp( gamma_v1*(rootcat4==1) +  gamma_v2*(rootcat4==2) +  gamma_v3*(rootcat4==3) +  
                                gamma_v4*(rootcat4==4) +  gamma_v5*(rootcat4==5) ), 
               
               v5      = exp( gamma_v1*(rootcat5==1) +  gamma_v2*(rootcat5==2) +  gamma_v3*(rootcat5==3) +  
                                gamma_v4*(rootcat5==4) +  gamma_v5*(rootcat5==5) ), 
               
               v6      = exp( gamma_v1*(rootcat6==1) +  gamma_v2*(rootcat6==2) +  gamma_v3*(rootcat6==3) +  
                                gamma_v4*(rootcat6==4) +  gamma_v5*(rootcat6==5) ),
               
               v7      = exp( gamma_v1*(rootcat7==1) +  gamma_v2*(rootcat7==2) +  gamma_v3*(rootcat7==3) +  
                                gamma_v4*(rootcat7==4) +  gamma_v5*(rootcat7==5) ), 
               
               v8      = exp( gamma_v1*(rootcat8==1) +  gamma_v2*(rootcat8==2) +  gamma_v3*(rootcat8==3) +  
                                gamma_v4*(rootcat8==4) +  gamma_v5*(rootcat8==5) )

  )
  
  delta = list( list(0,0,0,0,0,0,0 ,0 ) , 
                list(   d12 * ( (rootcat2 == 1) * (rootcat1 == 2) + (rootcat2 == 2) * (rootcat1 == 1) )    +
                          d13 * ( (rootcat2 == 1) * (rootcat1 == 3) + (rootcat2 == 3) * (rootcat1 == 1) )  +
                          d14 * ( (rootcat2 == 1) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 1) )  +
                          d15 * ( (rootcat2 == 1) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 1) )  +
                          d23 * ( (rootcat2 == 2) * (rootcat1 == 3) + (rootcat2 == 3) * (rootcat1 == 2) )  +
                          d24 * ( (rootcat2 == 2) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 2) )  +
                          d25 * ( (rootcat2 == 2) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 2) )  +
                          d34 * ( (rootcat2 == 3) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 3) )  +
                          d35 * ( (rootcat2 == 3) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 3) )  +
                          d45 * ( (rootcat2 == 4) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 4) )  ,
                        0,0,0,0,0,0,0 ) ,
                
                list(   d12 * ( (rootcat3 == 1) * (rootcat1 == 2) + (rootcat3 == 2) * (rootcat1 == 1) )    +
                          d13 * ( (rootcat3 == 1) * (rootcat1 == 3) + (rootcat3 == 3) * (rootcat1 == 1) )  +
                          d14 * ( (rootcat3 == 1) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 1) )  +
                          d15 * ( (rootcat3 == 1) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 1) )  +
                          d23 * ( (rootcat3 == 2) * (rootcat1 == 3) + (rootcat3 == 3) * (rootcat1 == 2) )  +
                          d24 * ( (rootcat3 == 2) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 2) )  +
                          d25 * ( (rootcat3 == 2) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 2) )  +
                          d34 * ( (rootcat3 == 3) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 3) )  +
                          d35 * ( (rootcat3 == 3) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 3) )  +
                          d45 * ( (rootcat3 == 4) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 4) )  ,
                        
                        d12 * ( (rootcat3 == 1) * (rootcat2 == 2) + (rootcat3 == 2) * (rootcat2 == 1) )    +
                          d13 * ( (rootcat3 == 1) * (rootcat2 == 3) + (rootcat3 == 3) * (rootcat2 == 1) )  +
                          d14 * ( (rootcat3 == 1) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 1) )  +
                          d15 * ( (rootcat3 == 1) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 1) )  +
                          d23 * ( (rootcat3 == 2) * (rootcat2 == 3) + (rootcat3 == 3) * (rootcat2 == 2) )  +
                          d24 * ( (rootcat3 == 2) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 2) )  +
                          d25 * ( (rootcat3 == 2) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 2) )  +
                          d34 * ( (rootcat3 == 3) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 3) )  +
                          d35 * ( (rootcat3 == 3) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 3) )  +
                          d45 * ( (rootcat3 == 4) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 4) )  ,
                        0,0,0,0,0,0 ) , 
                
                list(     d12 * ( (rootcat4 == 1) * (rootcat1 == 2) + (rootcat4 == 2) * (rootcat1 == 1) )    +
                            d13 * ( (rootcat4 == 1) * (rootcat1 == 3) + (rootcat4 == 3) * (rootcat1 == 1) )  +
                            d14 * ( (rootcat4 == 1) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 1) )  +
                            d15 * ( (rootcat4 == 1) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 1) )  +
                            d23 * ( (rootcat4 == 2) * (rootcat1 == 3) + (rootcat4 == 3) * (rootcat1 == 2) )  +
                            d24 * ( (rootcat4 == 2) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 2) )  +
                            d25 * ( (rootcat4 == 2) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 2) )  +
                            d34 * ( (rootcat4 == 3) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 3) )  +
                            d35 * ( (rootcat4 == 3) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 3) )  +
                            d45 * ( (rootcat4 == 4) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 4) )  ,
                          
                          d12 * ( (rootcat4 == 1) * (rootcat2 == 2) + (rootcat4 == 2) * (rootcat2 == 1) )    +
                            d13 * ( (rootcat4 == 1) * (rootcat2 == 3) + (rootcat4 == 3) * (rootcat2 == 1) )  +
                            d14 * ( (rootcat4 == 1) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 1) )  +
                            d15 * ( (rootcat4 == 1) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 1) )  +
                            d23 * ( (rootcat4 == 2) * (rootcat2 == 3) + (rootcat4 == 3) * (rootcat2 == 2) )  +
                            d24 * ( (rootcat4 == 2) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 2) )  +
                            d25 * ( (rootcat4 == 2) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 2) )  +
                            d34 * ( (rootcat4 == 3) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 3) )  +
                            d35 * ( (rootcat4 == 3) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 3) )  +
                            d45 * ( (rootcat4 == 4) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 4) )  ,
                          
                          d12 * ( (rootcat4 == 1) * (rootcat3 == 2) + (rootcat4 == 2) * (rootcat3 == 1) )    +
                            d13 * ( (rootcat4 == 1) * (rootcat3 == 3) + (rootcat4 == 3) * (rootcat3 == 1) )  +
                            d14 * ( (rootcat4 == 1) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 1) )  +
                            d15 * ( (rootcat4 == 1) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 1) )  +
                            d23 * ( (rootcat4 == 2) * (rootcat3 == 3) + (rootcat4 == 3) * (rootcat3 == 2) )  +
                            d24 * ( (rootcat4 == 2) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 2) )  +
                            d25 * ( (rootcat4 == 2) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 2) )  +
                            d34 * ( (rootcat4 == 3) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 3) )  +
                            d35 * ( (rootcat4 == 3) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 3) )  +
                            d45 * ( (rootcat4 == 4) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 4) )  ,
                          0,0,0,0,0 ) ,
                
                list(       d12 * ( (rootcat5 == 1) * (rootcat1 == 2) + (rootcat5 == 2) * (rootcat1 == 1) )    +
                              d13 * ( (rootcat5 == 1) * (rootcat1 == 3) + (rootcat5 == 3) * (rootcat1 == 1) )  +
                              d14 * ( (rootcat5 == 1) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 1) )  +
                              d15 * ( (rootcat5 == 1) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 1) )  +
                              d23 * ( (rootcat5 == 2) * (rootcat1 == 3) + (rootcat5 == 3) * (rootcat1 == 2) )  +
                              d24 * ( (rootcat5 == 2) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 2) )  +
                              d25 * ( (rootcat5 == 2) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 2) )  +
                              d34 * ( (rootcat5 == 3) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 3) )  +
                              d35 * ( (rootcat5 == 3) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 3) )  +
                              d45 * ( (rootcat5 == 4) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 4) )  ,
                            
                            d12 * ( (rootcat5 == 1) * (rootcat2 == 2) + (rootcat5 == 2) * (rootcat2 == 1) )    +
                              d13 * ( (rootcat5 == 1) * (rootcat2 == 3) + (rootcat5 == 3) * (rootcat2 == 1) )  +
                              d14 * ( (rootcat5 == 1) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 1) )  +
                              d15 * ( (rootcat5 == 1) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 1) )  +
                              d23 * ( (rootcat5 == 2) * (rootcat2 == 3) + (rootcat5 == 3) * (rootcat2 == 2) )  +
                              d24 * ( (rootcat5 == 2) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 2) )  +
                              d25 * ( (rootcat5 == 2) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 2) )  +
                              d34 * ( (rootcat5 == 3) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 3) )  +
                              d35 * ( (rootcat5 == 3) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 3) )  +
                              d45 * ( (rootcat5 == 4) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 4) )  ,
                            
                            d12 * ( (rootcat5 == 1) * (rootcat3 == 2) + (rootcat5 == 2) * (rootcat3 == 1) )    +
                              d13 * ( (rootcat5 == 1) * (rootcat3 == 3) + (rootcat5 == 3) * (rootcat3 == 1) )  +
                              d14 * ( (rootcat5 == 1) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 1) )  +
                              d15 * ( (rootcat5 == 1) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 1) )  +
                              d23 * ( (rootcat5 == 2) * (rootcat3 == 3) + (rootcat5 == 3) * (rootcat3 == 2) )  +
                              d24 * ( (rootcat5 == 2) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 2) )  +
                              d25 * ( (rootcat5 == 2) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 2) )  +
                              d34 * ( (rootcat5 == 3) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 3) )  +
                              d35 * ( (rootcat5 == 3) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 3) )  +
                              d45 * ( (rootcat5 == 4) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 4) )  ,
                            
                            d12 * ( (rootcat5 == 1) * (rootcat4 == 2) + (rootcat5 == 2) * (rootcat4 == 1) )    +
                              d13 * ( (rootcat5 == 1) * (rootcat4 == 3) + (rootcat5 == 3) * (rootcat4 == 1) )  +
                              d14 * ( (rootcat5 == 1) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 1) )  +
                              d15 * ( (rootcat5 == 1) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 1) )  +
                              d23 * ( (rootcat5 == 2) * (rootcat4 == 3) + (rootcat5 == 3) * (rootcat4 == 2) )  +
                              d24 * ( (rootcat5 == 2) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 2) )  +
                              d25 * ( (rootcat5 == 2) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 2) )  +
                              d34 * ( (rootcat5 == 3) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 3) )  +
                              d35 * ( (rootcat5 == 3) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 3) )  +
                              d45 * ( (rootcat5 == 4) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 4) )  ,
                            
                            0,0,0,0) ,
                
                list(       d12 * ( (rootcat6 == 1) * (rootcat1 == 2) + (rootcat6 == 2) * (rootcat1 == 1) )    +
                              d13 * ( (rootcat6 == 1) * (rootcat1 == 3) + (rootcat6 == 3) * (rootcat1 == 1) )  +
                              d14 * ( (rootcat6 == 1) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 1) )  +
                              d15 * ( (rootcat6 == 1) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 1) )  +
                              d23 * ( (rootcat6 == 2) * (rootcat1 == 3) + (rootcat6 == 3) * (rootcat1 == 2) )  +
                              d24 * ( (rootcat6 == 2) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 2) )  +
                              d25 * ( (rootcat6 == 2) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 2) )  +
                              d34 * ( (rootcat6 == 3) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 3) )  +
                              d35 * ( (rootcat6 == 3) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 3) )  +
                              d45 * ( (rootcat6 == 4) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 4) )  ,
                            
                            d12 * ( (rootcat6 == 1) * (rootcat2 == 2) + (rootcat6 == 2) * (rootcat2 == 1) )    +
                              d13 * ( (rootcat6 == 1) * (rootcat2 == 3) + (rootcat6 == 3) * (rootcat2 == 1) )  +
                              d14 * ( (rootcat6 == 1) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 1) )  +
                              d15 * ( (rootcat6 == 1) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 1) )  +
                              d23 * ( (rootcat6 == 2) * (rootcat2 == 3) + (rootcat6 == 3) * (rootcat2 == 2) )  +
                              d24 * ( (rootcat6 == 2) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 2) )  +
                              d25 * ( (rootcat6 == 2) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 2) )  +
                              d34 * ( (rootcat6 == 3) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 3) )  +
                              d35 * ( (rootcat6 == 3) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 3) )  +
                              d45 * ( (rootcat6 == 4) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 4) )  ,
                            
                            d12 * ( (rootcat6 == 1) * (rootcat3 == 2) + (rootcat6 == 2) * (rootcat3 == 1) )    +
                              d13 * ( (rootcat6 == 1) * (rootcat3 == 3) + (rootcat6 == 3) * (rootcat3 == 1) )  +
                              d14 * ( (rootcat6 == 1) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 1) )  +
                              d15 * ( (rootcat6 == 1) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 1) )  +
                              d23 * ( (rootcat6 == 2) * (rootcat3 == 3) + (rootcat6 == 3) * (rootcat3 == 2) )  +
                              d24 * ( (rootcat6 == 2) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 2) )  +
                              d25 * ( (rootcat6 == 2) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 2) )  +
                              d34 * ( (rootcat6 == 3) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 3) )  +
                              d35 * ( (rootcat6 == 3) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 3) )  +
                              d45 * ( (rootcat6 == 4) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 4) )  ,
                            
                            d12 * ( (rootcat6 == 1) * (rootcat4 == 2) + (rootcat6 == 2) * (rootcat4 == 1) )    +
                              d13 * ( (rootcat6 == 1) * (rootcat4 == 3) + (rootcat6 == 3) * (rootcat4 == 1) )  +
                              d14 * ( (rootcat6 == 1) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 1) )  +
                              d15 * ( (rootcat6 == 1) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 1) )  +
                              d23 * ( (rootcat6 == 2) * (rootcat4 == 3) + (rootcat6 == 3) * (rootcat4 == 2) )  +
                              d24 * ( (rootcat6 == 2) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 2) )  +
                              d25 * ( (rootcat6 == 2) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 2) )  +
                              d34 * ( (rootcat6 == 3) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 3) )  +
                              d35 * ( (rootcat6 == 3) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 3) )  +
                              d45 * ( (rootcat6 == 4) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 4) )  ,
                            
                            d12 * ( (rootcat6 == 1) * (rootcat5 == 2) + (rootcat6 == 2) * (rootcat5 == 1) )    +
                              d13 * ( (rootcat6 == 1) * (rootcat5 == 3) + (rootcat6 == 3) * (rootcat5 == 1) )  +
                              d14 * ( (rootcat6 == 1) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 1) )  +
                              d15 * ( (rootcat6 == 1) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 1) )  +
                              d23 * ( (rootcat6 == 2) * (rootcat5 == 3) + (rootcat6 == 3) * (rootcat5 == 2) )  +
                              d24 * ( (rootcat6 == 2) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 2) )  +
                              d25 * ( (rootcat6 == 2) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 2) )  +
                              d34 * ( (rootcat6 == 3) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 3) )  +
                              d35 * ( (rootcat6 == 3) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 3) )  +
                              d45 * ( (rootcat6 == 4) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 4) )  ,
                            
                            0,0,0) , 
                
                list(       d12 * ( (rootcat7 == 1) * (rootcat1 == 2) + (rootcat7 == 2) * (rootcat1 == 1) )    +
                              d13 * ( (rootcat7 == 1) * (rootcat1 == 3) + (rootcat7 == 3) * (rootcat1 == 1) )  +
                              d14 * ( (rootcat7 == 1) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 1) )  +
                              d15 * ( (rootcat7 == 1) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 1) )  +
                              d23 * ( (rootcat7 == 2) * (rootcat1 == 3) + (rootcat7 == 3) * (rootcat1 == 2) )  +
                              d24 * ( (rootcat7 == 2) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 2) )  +
                              d25 * ( (rootcat7 == 2) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 2) )  +
                              d34 * ( (rootcat7 == 3) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 3) )  +
                              d35 * ( (rootcat7 == 3) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 3) )  +
                              d45 * ( (rootcat7 == 4) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 4) )  ,
                            
                            d12 * ( (rootcat7 == 1) * (rootcat2 == 2) + (rootcat7 == 2) * (rootcat2 == 1) )    +
                              d13 * ( (rootcat7 == 1) * (rootcat2 == 3) + (rootcat7 == 3) * (rootcat2 == 1) )  +
                              d14 * ( (rootcat7 == 1) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 1) )  +
                              d15 * ( (rootcat7 == 1) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 1) )  +
                              d23 * ( (rootcat7 == 2) * (rootcat2 == 3) + (rootcat7 == 3) * (rootcat2 == 2) )  +
                              d24 * ( (rootcat7 == 2) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 2) )  +
                              d25 * ( (rootcat7 == 2) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 2) )  +
                              d34 * ( (rootcat7 == 3) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 3) )  +
                              d35 * ( (rootcat7 == 3) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 3) )  +
                              d45 * ( (rootcat7 == 4) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 4) )  ,
                            
                            d12 * ( (rootcat7 == 1) * (rootcat3 == 2) + (rootcat7 == 2) * (rootcat3 == 1) )    +
                              d13 * ( (rootcat7 == 1) * (rootcat3 == 3) + (rootcat7 == 3) * (rootcat3 == 1) )  +
                              d14 * ( (rootcat7 == 1) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 1) )  +
                              d15 * ( (rootcat7 == 1) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 1) )  +
                              d23 * ( (rootcat7 == 2) * (rootcat3 == 3) + (rootcat7 == 3) * (rootcat3 == 2) )  +
                              d24 * ( (rootcat7 == 2) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 2) )  +
                              d25 * ( (rootcat7 == 2) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 2) )  +
                              d34 * ( (rootcat7 == 3) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 3) )  +
                              d35 * ( (rootcat7 == 3) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 3) )  +
                              d45 * ( (rootcat7 == 4) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 4) )  ,
                            
                            d12 * ( (rootcat7 == 1) * (rootcat4 == 2) + (rootcat7 == 2) * (rootcat4 == 1) )    +
                              d13 * ( (rootcat7 == 1) * (rootcat4 == 3) + (rootcat7 == 3) * (rootcat4 == 1) )  +
                              d14 * ( (rootcat7 == 1) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 1) )  +
                              d15 * ( (rootcat7 == 1) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 1) )  +
                              d23 * ( (rootcat7 == 2) * (rootcat4 == 3) + (rootcat7 == 3) * (rootcat4 == 2) )  +
                              d24 * ( (rootcat7 == 2) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 2) )  +
                              d25 * ( (rootcat7 == 2) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 2) )  +
                              d34 * ( (rootcat7 == 3) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 3) )  +
                              d35 * ( (rootcat7 == 3) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 3) )  +
                              d45 * ( (rootcat7 == 4) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 4) )  ,
                            
                            d12 * ( (rootcat7 == 1) * (rootcat5 == 2) + (rootcat7 == 2) * (rootcat5 == 1) )    +
                              d13 * ( (rootcat7 == 1) * (rootcat5 == 3) + (rootcat7 == 3) * (rootcat5 == 1) )  +
                              d14 * ( (rootcat7 == 1) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 1) )  +
                              d15 * ( (rootcat7 == 1) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 1) )  +
                              d23 * ( (rootcat7 == 2) * (rootcat5 == 3) + (rootcat7 == 3) * (rootcat5 == 2) )  +
                              d24 * ( (rootcat7 == 2) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 2) )  +
                              d25 * ( (rootcat7 == 2) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 2) )  +
                              d34 * ( (rootcat7 == 3) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 3) )  +
                              d35 * ( (rootcat7 == 3) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 3) )  +
                              d45 * ( (rootcat7 == 4) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 4) )  ,
                            
                            d12 * ( (rootcat7 == 1) * (rootcat6 == 2) + (rootcat7 == 2) * (rootcat6 == 1) )    +
                              d13 * ( (rootcat7 == 1) * (rootcat6 == 3) + (rootcat7 == 3) * (rootcat6 == 1) )  +
                              d14 * ( (rootcat7 == 1) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 1) )  +
                              d15 * ( (rootcat7 == 1) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 1) )  +
                              d23 * ( (rootcat7 == 2) * (rootcat6 == 3) + (rootcat7 == 3) * (rootcat6 == 2) )  +
                              d24 * ( (rootcat7 == 2) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 2) )  +
                              d25 * ( (rootcat7 == 2) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 2) )  +
                              d34 * ( (rootcat7 == 3) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 3) )  +
                              d35 * ( (rootcat7 == 3) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 3) )  +
                              d45 * ( (rootcat7 == 4) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 4) )  ,
                            
                            0,0 ) , 
                
                list(       d12 * ( (rootcat8 == 1) * (rootcat1 == 2) + (rootcat8 == 2) * (rootcat1 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat1 == 3) + (rootcat8 == 3) * (rootcat1 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat1 == 3) + (rootcat8 == 3) * (rootcat1 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 4) )  ,
                            
                            d12 * ( (rootcat8 == 1) * (rootcat2 == 2) + (rootcat8 == 2) * (rootcat2 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat2 == 3) + (rootcat8 == 3) * (rootcat2 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat2 == 3) + (rootcat8 == 3) * (rootcat2 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 4) )  ,
                            
                            d12 * ( (rootcat8 == 1) * (rootcat3 == 2) + (rootcat8 == 2) * (rootcat3 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat3 == 3) + (rootcat8 == 3) * (rootcat3 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat3 == 3) + (rootcat8 == 3) * (rootcat3 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 4) )  ,
                            
                            d12 * ( (rootcat8 == 1) * (rootcat4 == 2) + (rootcat8 == 2) * (rootcat4 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat4 == 3) + (rootcat8 == 3) * (rootcat4 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat4 == 3) + (rootcat8 == 3) * (rootcat4 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 4) )  ,
                            
                            d12 * ( (rootcat8 == 1) * (rootcat5 == 2) + (rootcat8 == 2) * (rootcat5 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat5 == 3) + (rootcat8 == 3) * (rootcat5 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat5 == 3) + (rootcat8 == 3) * (rootcat5 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 4) )  ,
                            
                            d12 * ( (rootcat8 == 1) * (rootcat6 == 2) + (rootcat8 == 2) * (rootcat6 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat6 == 3) + (rootcat8 == 3) * (rootcat6 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat6 == 3) + (rootcat8 == 3) * (rootcat6 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 4) )  ,
                            
                            d12 * ( (rootcat8 == 1) * (rootcat7 == 2) + (rootcat8 == 2) * (rootcat7 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat7 == 3) + (rootcat8 == 3) * (rootcat7 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat7 == 3) + (rootcat8 == 3) * (rootcat7 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 4) )  ,
                            0 ) 
  )
  
  
  
  emdc_settings$utilityOutside = delta_outside
  emdc_settings$utilities = V
  emdc_settings$gamma = gamma
  emdc_settings$delta = delta
  
  
  # 
  P[["model"]] = apollo_emdc(emdc_settings, functionality)
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs ,functionality)
  
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

apollo_prediction = function (model, apollo_probabilities = apollo_probabilities1, apollo_inputs, prediction_settings = list(), modelComponent = NA) {
  prediction_settings = list()
  modelComponent = NA
  
  if (is.null(prediction_settings$modelComponent)) {
    if (exists("modelComponent")) 
      prediction_settings$modelComponent = modelComponent
    else prediction_settings$modelComponent = NA
  }
  if (is.null(prediction_settings$runs))   prediction_settings$runs = 1
  if (is.null(prediction_settings$silent))   prediction_settings$silent = FALSE
  silent = prediction_settings$silent
  if (!is.null(apollo_inputs$silent) && apollo_inputs$silent) silent = TRUE
  
  if (is.null(prediction_settings$nRep)) prediction_settings$nRep <- 100L
  
  if (is.null(prediction_settings$summary)) prediction_settings$summary <- TRUE
  
  modelComponent = prediction_settings$modelComponent
  runs = prediction_settings$runs
  apollo_inputs$nRep <- prediction_settings$nRep
 
  apollo_compareInputs(apollo_inputs)
  apollo_randCoeff = apollo_inputs[["apollo_randCoeff"]]
  apollo_lcPars = apollo_inputs[["apollo_lcPars"]]
  apollo_checkArguments(apollo_probabilities, apollo_randCoeff, 
                        apollo_lcPars)
  if (!silent)     apollo_print("Running predictions from model using parameter estimates...")
  

  predictions = apollo_probabilities(apollo_beta, apollo_inputs, 
                                     functionality = "prediction")
  

  predictionamount = predictions$model[,1:9]
  
  return(predictionamount)
}

apollo_inputs = apollo_validateInputs(apollo_beta = apollo_beta1)

predictionprob1 = apollo_prediction(model, apollo_probabilities = apollo_probabilities1, apollo_inputs, prediction_settings=list(runs=1), modelComponent = NA)




#######################Segment 2

apollo_beta = apollo_beta2

apollo_probabilities2=function(apollo_beta = apollo_beta2, apollo_inputs, functionality="estimate"){
  
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### Define individual alternatives
  alternatives = c("v1", 
                   "v2", 
                   "v3", 
                   "v4", 
                   "v5", 
                   "v6", 
                   "v7", 
                   "v8" #, 
                   # "outside"
  )
  
  ### Define availabilities
  avail = list(v1  = availnew1,    ######Need to check whether it works here
               v2     = availnew2,
               v3   = availnew3,
               v4 = availnew4,
               v5 = availnew5,
               v6   = availnew6,
               v7 = availnew7,
               v8 = availnew8 #,
               #outside  = availnew9
  )
  
  ### Define continuous consumption for individual alternatives
  continuousChoice = list(v1  =clicknew1 * availnew1,
                          v2     =clicknew2 * availnew2,
                          v3   =clicknew3 * availnew3,
                          v4 =clicknew4 * availnew4,
                          v5 =clicknew5 * availnew5,
                          v6   =clicknew6 * availnew6,
                          v7 =clicknew7 * availnew7,
                          v8 =clicknew8 * availnew8#,
                          #outside  =outsideclick
  )
  
  
  
  ### Define alpha parameters
  alpha = list(v1  = 1e-3 , 
               v2     = 1e-3 , 
               v3   = 1e-3 , 
               v4 = 1e-3 , 
               v5 = 1e-3 ,
               v6   = 1e-3 ,
               v7  = 1e-3 , 
               v8 = 1e-3 #, 
               #outside  = 1e-3 
  )
  
  
  ### Define costs for individual alternatives
  cost = list(v1      = 1, 
              v2         = 1,
              v3       = 1,
              v4     = 1,
              v5 = 1,
              v6       = 1, 
              v7      = 1,
              v8     = 1 #,
              #outside       = 1
  )
  
  ### Define budget
  budget = clicktotal

  emdc_settings <- list(continuousChoice = continuousChoice, 
                        avail            = avail,
                        budget           = budget,
                        sigma            = 0.99, 
                        cost             = cost)
 
  ### ### Compute class-specific utilities
  V = list()
  
  V[["v1"    ]] =   delta_v1*(rootcat1==1) +  delta_v2*(rootcat1==2) +  delta_v3*(rootcat1==3) +  
    delta_v4*(rootcat1==4) +  delta_v5*(rootcat1==5) + 
    bweekend*weekend + bnpv*npvnew1 + bnclick*nclicknew1 + bnbuy*nbuynew1 + 
    nrccpv*rccpvnew1 + nrccclick*rccclicknew1 + nrccbuy*rccbuynew1 + badstockpsai*ad11 + bresidual*res1 + bvariety*Lag_Distinct
  
  V[["v2"    ]] = delta_v1*(rootcat2==1) +  delta_v2*(rootcat2==2) +  delta_v3*(rootcat2==3) +  
    delta_v4*(rootcat2==4) +  delta_v5*(rootcat2==5) + 
    bweekend*weekend + bnpv*npvnew2 + bnclick*nclicknew2 + bnbuy*nbuynew2 + 
    nrccpv*rccpvnew2 + nrccclick*rccclicknew2 + nrccbuy*rccbuynew2 + badstockpsai*ad21 + bresidual*res2 + bvariety*Lag_Distinct
  
  V[["v3"  ]] = delta_v1*(rootcat3==1) +  delta_v2*(rootcat3==2) +  delta_v3*(rootcat3==3) +  
    delta_v4*(rootcat3==4) +  delta_v5*(rootcat3==5) + 
    bweekend*weekend + bnpv*npvnew3 + bnclick*nclicknew3 + bnbuy*nbuynew3 + 
    nrccpv*rccpvnew3 + nrccclick*rccclicknew3 + nrccbuy*rccbuynew3 + badstockpsai*ad31 + bresidual*res3 + bvariety*Lag_Distinct
  
  V[["v4"]] = delta_v1*(rootcat4==1) +  delta_v2*(rootcat4==2) +  delta_v3*(rootcat4==3) +  
    delta_v4*(rootcat4==4) +  delta_v5*(rootcat4==5) + 
    bweekend*weekend + bnpv*npvnew4 + bnclick*nclicknew4 + bnbuy*nbuynew4 + 
    nrccpv*rccpvnew4 + nrccclick*rccclicknew4 + nrccbuy*rccbuynew4 + badstockpsai*ad41 + bresidual*res4 + bvariety*Lag_Distinct
  
  V[["v5"]] = delta_v1*(rootcat5==1) +  delta_v2*(rootcat5==2) +  delta_v3*(rootcat5==3) +  
    delta_v4*(rootcat5==4) +  delta_v5*(rootcat5==5) + 
    bweekend*weekend + bnpv*npvnew5 + bnclick*nclicknew5 + bnbuy*nbuynew5 + 
    nrccpv*rccpvnew5 + nrccclick*rccclicknew5 + nrccbuy*rccbuynew5 + badstockpsai*ad51 + bresidual*res5 + bvariety*Lag_Distinct
  
  V[["v6"  ]] = delta_v1*(rootcat6==1) +  delta_v2*(rootcat6==2) +  delta_v3*(rootcat6==3) +  
    delta_v4*(rootcat6==4) +  delta_v5*(rootcat6==5) + 
    bweekend*weekend + bnpv*npvnew6 + bnclick*nclicknew6 + bnbuy*nbuynew6 + 
    nrccpv*rccpvnew6 + nrccclick*rccclicknew6 + nrccbuy*rccbuynew6 + badstockpsai*ad61 + bresidual*res6 + bvariety*Lag_Distinct
  
  V[["v7"]] = delta_v1*(rootcat7==1) +  delta_v2*(rootcat7==2) +  delta_v3*(rootcat7==3) +  
    delta_v4*(rootcat7==4) +  delta_v5*(rootcat7==5) + 
    bweekend*weekend + bnpv*npvnew7 + bnclick*nclicknew7 + bnbuy*nbuynew7 + 
    nrccpv*rccpvnew7 + nrccclick*rccclicknew7 + nrccbuy*rccbuynew7 + badstockpsai*ad71 + bresidual*res7 + bvariety*Lag_Distinct
  
  V[["v8"]] = delta_v1*(rootcat8==1) +  delta_v2*(rootcat8==2) +  delta_v3*(rootcat8==3) +  
    delta_v4*(rootcat8==4) +  delta_v5*(rootcat8==5) + 
    bweekend*weekend + bnpv*npvnew8 + bnclick*nclicknew8 + bnbuy*nbuynew8 + 
    nrccpv*rccpvnew8 + nrccclick*rccclicknew8 + nrccbuy*rccbuynew8 + badstockpsai*ad81 + bresidual*res8 + bvariety*Lag_Distinct

  ### Define gamma parameters
  gamma = list(v1      = exp( gamma_v1*(rootcat1==1) +  gamma_v2*(rootcat1==2) +  gamma_v3*(rootcat1==3) +  
                                gamma_v4*(rootcat1==4) +  gamma_v5*(rootcat1==5) ), 
               
               v2      = exp( gamma_v1*(rootcat2==1) +  gamma_v2*(rootcat2==2) +  gamma_v3*(rootcat2==3) +  
                                gamma_v4*(rootcat2==4) +  gamma_v5*(rootcat2==5) ), 
               
               
               v3      = exp( gamma_v1*(rootcat3==1) +  gamma_v2*(rootcat3==2) +  gamma_v3*(rootcat3==3) +  
                                gamma_v4*(rootcat3==4) +  gamma_v5*(rootcat3==5) ), 
               
               v4      = exp( gamma_v1*(rootcat4==1) +  gamma_v2*(rootcat4==2) +  gamma_v3*(rootcat4==3) +  
                                gamma_v4*(rootcat4==4) +  gamma_v5*(rootcat4==5) ), 
               
               v5      = exp( gamma_v1*(rootcat5==1) +  gamma_v2*(rootcat5==2) +  gamma_v3*(rootcat5==3) +  
                                gamma_v4*(rootcat5==4) +  gamma_v5*(rootcat5==5) ), 
               
               v6      = exp( gamma_v1*(rootcat6==1) +  gamma_v2*(rootcat6==2) +  gamma_v3*(rootcat6==3) +  
                                gamma_v4*(rootcat6==4) +  gamma_v5*(rootcat6==5) ),
               
               v7      = exp( gamma_v1*(rootcat7==1) +  gamma_v2*(rootcat7==2) +  gamma_v3*(rootcat7==3) +  
                                gamma_v4*(rootcat7==4) +  gamma_v5*(rootcat7==5) ), 
               
               v8      = exp( gamma_v1*(rootcat8==1) +  gamma_v2*(rootcat8==2) +  gamma_v3*(rootcat8==3) +  
                                gamma_v4*(rootcat8==4) +  gamma_v5*(rootcat8==5) )
  )
  
  delta = list( list(0,0,0,0,0,0,0 ,0 ) , 
                list(   d12 * ( (rootcat2 == 1) * (rootcat1 == 2) + (rootcat2 == 2) * (rootcat1 == 1) )    +
                          d13 * ( (rootcat2 == 1) * (rootcat1 == 3) + (rootcat2 == 3) * (rootcat1 == 1) )  +
                          d14 * ( (rootcat2 == 1) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 1) )  +
                          d15 * ( (rootcat2 == 1) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 1) )  +
                          d23 * ( (rootcat2 == 2) * (rootcat1 == 3) + (rootcat2 == 3) * (rootcat1 == 2) )  +
                          d24 * ( (rootcat2 == 2) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 2) )  +
                          d25 * ( (rootcat2 == 2) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 2) )  +
                          d34 * ( (rootcat2 == 3) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 3) )  +
                          d35 * ( (rootcat2 == 3) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 3) )  +
                          d45 * ( (rootcat2 == 4) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 4) )  ,
                        0,0,0,0,0,0,0 ) ,
                
                list(   d12 * ( (rootcat3 == 1) * (rootcat1 == 2) + (rootcat3 == 2) * (rootcat1 == 1) )    +
                          d13 * ( (rootcat3 == 1) * (rootcat1 == 3) + (rootcat3 == 3) * (rootcat1 == 1) )  +
                          d14 * ( (rootcat3 == 1) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 1) )  +
                          d15 * ( (rootcat3 == 1) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 1) )  +
                          d23 * ( (rootcat3 == 2) * (rootcat1 == 3) + (rootcat3 == 3) * (rootcat1 == 2) )  +
                          d24 * ( (rootcat3 == 2) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 2) )  +
                          d25 * ( (rootcat3 == 2) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 2) )  +
                          d34 * ( (rootcat3 == 3) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 3) )  +
                          d35 * ( (rootcat3 == 3) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 3) )  +
                          d45 * ( (rootcat3 == 4) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 4) )  ,
                        
                        d12 * ( (rootcat3 == 1) * (rootcat2 == 2) + (rootcat3 == 2) * (rootcat2 == 1) )    +
                          d13 * ( (rootcat3 == 1) * (rootcat2 == 3) + (rootcat3 == 3) * (rootcat2 == 1) )  +
                          d14 * ( (rootcat3 == 1) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 1) )  +
                          d15 * ( (rootcat3 == 1) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 1) )  +
                          d23 * ( (rootcat3 == 2) * (rootcat2 == 3) + (rootcat3 == 3) * (rootcat2 == 2) )  +
                          d24 * ( (rootcat3 == 2) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 2) )  +
                          d25 * ( (rootcat3 == 2) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 2) )  +
                          d34 * ( (rootcat3 == 3) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 3) )  +
                          d35 * ( (rootcat3 == 3) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 3) )  +
                          d45 * ( (rootcat3 == 4) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 4) )  ,
                        0,0,0,0,0,0 ) , 
                
                list(     d12 * ( (rootcat4 == 1) * (rootcat1 == 2) + (rootcat4 == 2) * (rootcat1 == 1) )    +
                            d13 * ( (rootcat4 == 1) * (rootcat1 == 3) + (rootcat4 == 3) * (rootcat1 == 1) )  +
                            d14 * ( (rootcat4 == 1) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 1) )  +
                            d15 * ( (rootcat4 == 1) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 1) )  +
                            d23 * ( (rootcat4 == 2) * (rootcat1 == 3) + (rootcat4 == 3) * (rootcat1 == 2) )  +
                            d24 * ( (rootcat4 == 2) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 2) )  +
                            d25 * ( (rootcat4 == 2) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 2) )  +
                            d34 * ( (rootcat4 == 3) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 3) )  +
                            d35 * ( (rootcat4 == 3) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 3) )  +
                            d45 * ( (rootcat4 == 4) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 4) )  ,
                          
                          d12 * ( (rootcat4 == 1) * (rootcat2 == 2) + (rootcat4 == 2) * (rootcat2 == 1) )    +
                            d13 * ( (rootcat4 == 1) * (rootcat2 == 3) + (rootcat4 == 3) * (rootcat2 == 1) )  +
                            d14 * ( (rootcat4 == 1) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 1) )  +
                            d15 * ( (rootcat4 == 1) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 1) )  +
                            d23 * ( (rootcat4 == 2) * (rootcat2 == 3) + (rootcat4 == 3) * (rootcat2 == 2) )  +
                            d24 * ( (rootcat4 == 2) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 2) )  +
                            d25 * ( (rootcat4 == 2) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 2) )  +
                            d34 * ( (rootcat4 == 3) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 3) )  +
                            d35 * ( (rootcat4 == 3) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 3) )  +
                            d45 * ( (rootcat4 == 4) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 4) )  ,
                          
                          d12 * ( (rootcat4 == 1) * (rootcat3 == 2) + (rootcat4 == 2) * (rootcat3 == 1) )    +
                            d13 * ( (rootcat4 == 1) * (rootcat3 == 3) + (rootcat4 == 3) * (rootcat3 == 1) )  +
                            d14 * ( (rootcat4 == 1) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 1) )  +
                            d15 * ( (rootcat4 == 1) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 1) )  +
                            d23 * ( (rootcat4 == 2) * (rootcat3 == 3) + (rootcat4 == 3) * (rootcat3 == 2) )  +
                            d24 * ( (rootcat4 == 2) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 2) )  +
                            d25 * ( (rootcat4 == 2) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 2) )  +
                            d34 * ( (rootcat4 == 3) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 3) )  +
                            d35 * ( (rootcat4 == 3) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 3) )  +
                            d45 * ( (rootcat4 == 4) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 4) )  ,
                          0,0,0,0,0 ) ,
                
                list(       d12 * ( (rootcat5 == 1) * (rootcat1 == 2) + (rootcat5 == 2) * (rootcat1 == 1) )    +
                              d13 * ( (rootcat5 == 1) * (rootcat1 == 3) + (rootcat5 == 3) * (rootcat1 == 1) )  +
                              d14 * ( (rootcat5 == 1) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 1) )  +
                              d15 * ( (rootcat5 == 1) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 1) )  +
                              d23 * ( (rootcat5 == 2) * (rootcat1 == 3) + (rootcat5 == 3) * (rootcat1 == 2) )  +
                              d24 * ( (rootcat5 == 2) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 2) )  +
                              d25 * ( (rootcat5 == 2) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 2) )  +
                              d34 * ( (rootcat5 == 3) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 3) )  +
                              d35 * ( (rootcat5 == 3) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 3) )  +
                              d45 * ( (rootcat5 == 4) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 4) )  ,
                            
                            d12 * ( (rootcat5 == 1) * (rootcat2 == 2) + (rootcat5 == 2) * (rootcat2 == 1) )    +
                              d13 * ( (rootcat5 == 1) * (rootcat2 == 3) + (rootcat5 == 3) * (rootcat2 == 1) )  +
                              d14 * ( (rootcat5 == 1) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 1) )  +
                              d15 * ( (rootcat5 == 1) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 1) )  +
                              d23 * ( (rootcat5 == 2) * (rootcat2 == 3) + (rootcat5 == 3) * (rootcat2 == 2) )  +
                              d24 * ( (rootcat5 == 2) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 2) )  +
                              d25 * ( (rootcat5 == 2) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 2) )  +
                              d34 * ( (rootcat5 == 3) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 3) )  +
                              d35 * ( (rootcat5 == 3) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 3) )  +
                              d45 * ( (rootcat5 == 4) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 4) )  ,
                            
                            d12 * ( (rootcat5 == 1) * (rootcat3 == 2) + (rootcat5 == 2) * (rootcat3 == 1) )    +
                              d13 * ( (rootcat5 == 1) * (rootcat3 == 3) + (rootcat5 == 3) * (rootcat3 == 1) )  +
                              d14 * ( (rootcat5 == 1) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 1) )  +
                              d15 * ( (rootcat5 == 1) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 1) )  +
                              d23 * ( (rootcat5 == 2) * (rootcat3 == 3) + (rootcat5 == 3) * (rootcat3 == 2) )  +
                              d24 * ( (rootcat5 == 2) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 2) )  +
                              d25 * ( (rootcat5 == 2) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 2) )  +
                              d34 * ( (rootcat5 == 3) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 3) )  +
                              d35 * ( (rootcat5 == 3) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 3) )  +
                              d45 * ( (rootcat5 == 4) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 4) )  ,
                            
                            d12 * ( (rootcat5 == 1) * (rootcat4 == 2) + (rootcat5 == 2) * (rootcat4 == 1) )    +
                              d13 * ( (rootcat5 == 1) * (rootcat4 == 3) + (rootcat5 == 3) * (rootcat4 == 1) )  +
                              d14 * ( (rootcat5 == 1) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 1) )  +
                              d15 * ( (rootcat5 == 1) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 1) )  +
                              d23 * ( (rootcat5 == 2) * (rootcat4 == 3) + (rootcat5 == 3) * (rootcat4 == 2) )  +
                              d24 * ( (rootcat5 == 2) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 2) )  +
                              d25 * ( (rootcat5 == 2) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 2) )  +
                              d34 * ( (rootcat5 == 3) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 3) )  +
                              d35 * ( (rootcat5 == 3) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 3) )  +
                              d45 * ( (rootcat5 == 4) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 4) )  ,
                            
                            0,0,0,0) ,
                
                list(       d12 * ( (rootcat6 == 1) * (rootcat1 == 2) + (rootcat6 == 2) * (rootcat1 == 1) )    +
                              d13 * ( (rootcat6 == 1) * (rootcat1 == 3) + (rootcat6 == 3) * (rootcat1 == 1) )  +
                              d14 * ( (rootcat6 == 1) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 1) )  +
                              d15 * ( (rootcat6 == 1) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 1) )  +
                              d23 * ( (rootcat6 == 2) * (rootcat1 == 3) + (rootcat6 == 3) * (rootcat1 == 2) )  +
                              d24 * ( (rootcat6 == 2) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 2) )  +
                              d25 * ( (rootcat6 == 2) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 2) )  +
                              d34 * ( (rootcat6 == 3) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 3) )  +
                              d35 * ( (rootcat6 == 3) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 3) )  +
                              d45 * ( (rootcat6 == 4) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 4) )  ,
                            
                            d12 * ( (rootcat6 == 1) * (rootcat2 == 2) + (rootcat6 == 2) * (rootcat2 == 1) )    +
                              d13 * ( (rootcat6 == 1) * (rootcat2 == 3) + (rootcat6 == 3) * (rootcat2 == 1) )  +
                              d14 * ( (rootcat6 == 1) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 1) )  +
                              d15 * ( (rootcat6 == 1) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 1) )  +
                              d23 * ( (rootcat6 == 2) * (rootcat2 == 3) + (rootcat6 == 3) * (rootcat2 == 2) )  +
                              d24 * ( (rootcat6 == 2) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 2) )  +
                              d25 * ( (rootcat6 == 2) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 2) )  +
                              d34 * ( (rootcat6 == 3) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 3) )  +
                              d35 * ( (rootcat6 == 3) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 3) )  +
                              d45 * ( (rootcat6 == 4) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 4) )  ,
                            
                            d12 * ( (rootcat6 == 1) * (rootcat3 == 2) + (rootcat6 == 2) * (rootcat3 == 1) )    +
                              d13 * ( (rootcat6 == 1) * (rootcat3 == 3) + (rootcat6 == 3) * (rootcat3 == 1) )  +
                              d14 * ( (rootcat6 == 1) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 1) )  +
                              d15 * ( (rootcat6 == 1) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 1) )  +
                              d23 * ( (rootcat6 == 2) * (rootcat3 == 3) + (rootcat6 == 3) * (rootcat3 == 2) )  +
                              d24 * ( (rootcat6 == 2) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 2) )  +
                              d25 * ( (rootcat6 == 2) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 2) )  +
                              d34 * ( (rootcat6 == 3) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 3) )  +
                              d35 * ( (rootcat6 == 3) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 3) )  +
                              d45 * ( (rootcat6 == 4) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 4) )  ,
                            
                            d12 * ( (rootcat6 == 1) * (rootcat4 == 2) + (rootcat6 == 2) * (rootcat4 == 1) )    +
                              d13 * ( (rootcat6 == 1) * (rootcat4 == 3) + (rootcat6 == 3) * (rootcat4 == 1) )  +
                              d14 * ( (rootcat6 == 1) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 1) )  +
                              d15 * ( (rootcat6 == 1) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 1) )  +
                              d23 * ( (rootcat6 == 2) * (rootcat4 == 3) + (rootcat6 == 3) * (rootcat4 == 2) )  +
                              d24 * ( (rootcat6 == 2) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 2) )  +
                              d25 * ( (rootcat6 == 2) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 2) )  +
                              d34 * ( (rootcat6 == 3) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 3) )  +
                              d35 * ( (rootcat6 == 3) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 3) )  +
                              d45 * ( (rootcat6 == 4) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 4) )  ,
                            
                            d12 * ( (rootcat6 == 1) * (rootcat5 == 2) + (rootcat6 == 2) * (rootcat5 == 1) )    +
                              d13 * ( (rootcat6 == 1) * (rootcat5 == 3) + (rootcat6 == 3) * (rootcat5 == 1) )  +
                              d14 * ( (rootcat6 == 1) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 1) )  +
                              d15 * ( (rootcat6 == 1) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 1) )  +
                              d23 * ( (rootcat6 == 2) * (rootcat5 == 3) + (rootcat6 == 3) * (rootcat5 == 2) )  +
                              d24 * ( (rootcat6 == 2) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 2) )  +
                              d25 * ( (rootcat6 == 2) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 2) )  +
                              d34 * ( (rootcat6 == 3) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 3) )  +
                              d35 * ( (rootcat6 == 3) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 3) )  +
                              d45 * ( (rootcat6 == 4) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 4) )  ,
                            
                            0,0,0) , 
                
                list(       d12 * ( (rootcat7 == 1) * (rootcat1 == 2) + (rootcat7 == 2) * (rootcat1 == 1) )    +
                              d13 * ( (rootcat7 == 1) * (rootcat1 == 3) + (rootcat7 == 3) * (rootcat1 == 1) )  +
                              d14 * ( (rootcat7 == 1) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 1) )  +
                              d15 * ( (rootcat7 == 1) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 1) )  +
                              d23 * ( (rootcat7 == 2) * (rootcat1 == 3) + (rootcat7 == 3) * (rootcat1 == 2) )  +
                              d24 * ( (rootcat7 == 2) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 2) )  +
                              d25 * ( (rootcat7 == 2) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 2) )  +
                              d34 * ( (rootcat7 == 3) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 3) )  +
                              d35 * ( (rootcat7 == 3) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 3) )  +
                              d45 * ( (rootcat7 == 4) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 4) )  ,
                            
                            d12 * ( (rootcat7 == 1) * (rootcat2 == 2) + (rootcat7 == 2) * (rootcat2 == 1) )    +
                              d13 * ( (rootcat7 == 1) * (rootcat2 == 3) + (rootcat7 == 3) * (rootcat2 == 1) )  +
                              d14 * ( (rootcat7 == 1) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 1) )  +
                              d15 * ( (rootcat7 == 1) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 1) )  +
                              d23 * ( (rootcat7 == 2) * (rootcat2 == 3) + (rootcat7 == 3) * (rootcat2 == 2) )  +
                              d24 * ( (rootcat7 == 2) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 2) )  +
                              d25 * ( (rootcat7 == 2) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 2) )  +
                              d34 * ( (rootcat7 == 3) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 3) )  +
                              d35 * ( (rootcat7 == 3) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 3) )  +
                              d45 * ( (rootcat7 == 4) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 4) )  ,
                            
                            d12 * ( (rootcat7 == 1) * (rootcat3 == 2) + (rootcat7 == 2) * (rootcat3 == 1) )    +
                              d13 * ( (rootcat7 == 1) * (rootcat3 == 3) + (rootcat7 == 3) * (rootcat3 == 1) )  +
                              d14 * ( (rootcat7 == 1) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 1) )  +
                              d15 * ( (rootcat7 == 1) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 1) )  +
                              d23 * ( (rootcat7 == 2) * (rootcat3 == 3) + (rootcat7 == 3) * (rootcat3 == 2) )  +
                              d24 * ( (rootcat7 == 2) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 2) )  +
                              d25 * ( (rootcat7 == 2) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 2) )  +
                              d34 * ( (rootcat7 == 3) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 3) )  +
                              d35 * ( (rootcat7 == 3) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 3) )  +
                              d45 * ( (rootcat7 == 4) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 4) )  ,
                            
                            d12 * ( (rootcat7 == 1) * (rootcat4 == 2) + (rootcat7 == 2) * (rootcat4 == 1) )    +
                              d13 * ( (rootcat7 == 1) * (rootcat4 == 3) + (rootcat7 == 3) * (rootcat4 == 1) )  +
                              d14 * ( (rootcat7 == 1) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 1) )  +
                              d15 * ( (rootcat7 == 1) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 1) )  +
                              d23 * ( (rootcat7 == 2) * (rootcat4 == 3) + (rootcat7 == 3) * (rootcat4 == 2) )  +
                              d24 * ( (rootcat7 == 2) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 2) )  +
                              d25 * ( (rootcat7 == 2) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 2) )  +
                              d34 * ( (rootcat7 == 3) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 3) )  +
                              d35 * ( (rootcat7 == 3) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 3) )  +
                              d45 * ( (rootcat7 == 4) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 4) )  ,
                            
                            d12 * ( (rootcat7 == 1) * (rootcat5 == 2) + (rootcat7 == 2) * (rootcat5 == 1) )    +
                              d13 * ( (rootcat7 == 1) * (rootcat5 == 3) + (rootcat7 == 3) * (rootcat5 == 1) )  +
                              d14 * ( (rootcat7 == 1) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 1) )  +
                              d15 * ( (rootcat7 == 1) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 1) )  +
                              d23 * ( (rootcat7 == 2) * (rootcat5 == 3) + (rootcat7 == 3) * (rootcat5 == 2) )  +
                              d24 * ( (rootcat7 == 2) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 2) )  +
                              d25 * ( (rootcat7 == 2) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 2) )  +
                              d34 * ( (rootcat7 == 3) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 3) )  +
                              d35 * ( (rootcat7 == 3) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 3) )  +
                              d45 * ( (rootcat7 == 4) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 4) )  ,
                            
                            d12 * ( (rootcat7 == 1) * (rootcat6 == 2) + (rootcat7 == 2) * (rootcat6 == 1) )    +
                              d13 * ( (rootcat7 == 1) * (rootcat6 == 3) + (rootcat7 == 3) * (rootcat6 == 1) )  +
                              d14 * ( (rootcat7 == 1) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 1) )  +
                              d15 * ( (rootcat7 == 1) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 1) )  +
                              d23 * ( (rootcat7 == 2) * (rootcat6 == 3) + (rootcat7 == 3) * (rootcat6 == 2) )  +
                              d24 * ( (rootcat7 == 2) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 2) )  +
                              d25 * ( (rootcat7 == 2) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 2) )  +
                              d34 * ( (rootcat7 == 3) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 3) )  +
                              d35 * ( (rootcat7 == 3) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 3) )  +
                              d45 * ( (rootcat7 == 4) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 4) )  ,
                            
                            0,0 ) , 
                
                list(       d12 * ( (rootcat8 == 1) * (rootcat1 == 2) + (rootcat8 == 2) * (rootcat1 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat1 == 3) + (rootcat8 == 3) * (rootcat1 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat1 == 3) + (rootcat8 == 3) * (rootcat1 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 4) )  ,
                            
                            d12 * ( (rootcat8 == 1) * (rootcat2 == 2) + (rootcat8 == 2) * (rootcat2 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat2 == 3) + (rootcat8 == 3) * (rootcat2 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat2 == 3) + (rootcat8 == 3) * (rootcat2 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 4) )  ,
                            
                            d12 * ( (rootcat8 == 1) * (rootcat3 == 2) + (rootcat8 == 2) * (rootcat3 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat3 == 3) + (rootcat8 == 3) * (rootcat3 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat3 == 3) + (rootcat8 == 3) * (rootcat3 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 4) )  ,
                            
                            d12 * ( (rootcat8 == 1) * (rootcat4 == 2) + (rootcat8 == 2) * (rootcat4 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat4 == 3) + (rootcat8 == 3) * (rootcat4 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat4 == 3) + (rootcat8 == 3) * (rootcat4 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 4) )  ,
                            
                            d12 * ( (rootcat8 == 1) * (rootcat5 == 2) + (rootcat8 == 2) * (rootcat5 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat5 == 3) + (rootcat8 == 3) * (rootcat5 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat5 == 3) + (rootcat8 == 3) * (rootcat5 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 4) )  ,
                            
                            d12 * ( (rootcat8 == 1) * (rootcat6 == 2) + (rootcat8 == 2) * (rootcat6 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat6 == 3) + (rootcat8 == 3) * (rootcat6 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat6 == 3) + (rootcat8 == 3) * (rootcat6 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 4) )  ,
                            
                            d12 * ( (rootcat8 == 1) * (rootcat7 == 2) + (rootcat8 == 2) * (rootcat7 == 1) )    +
                              d13 * ( (rootcat8 == 1) * (rootcat7 == 3) + (rootcat8 == 3) * (rootcat7 == 1) )  +
                              d14 * ( (rootcat8 == 1) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 1) )  +
                              d15 * ( (rootcat8 == 1) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 1) )  +
                              d23 * ( (rootcat8 == 2) * (rootcat7 == 3) + (rootcat8 == 3) * (rootcat7 == 2) )  +
                              d24 * ( (rootcat8 == 2) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 2) )  +
                              d25 * ( (rootcat8 == 2) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 2) )  +
                              d34 * ( (rootcat8 == 3) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 3) )  +
                              d35 * ( (rootcat8 == 3) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 3) )  +
                              d45 * ( (rootcat8 == 4) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 4) )  ,
                            0 ) 
  )
  
  
  
  emdc_settings$utilityOutside = delta_outside
  emdc_settings$utilities = V
  emdc_settings$gamma = gamma
  emdc_settings$delta = delta
 
  # 
  P[["model"]] = apollo_emdc(emdc_settings, functionality)
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs ,functionality)
  
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

apollo_prediction = function (model, apollo_probabilities = apollo_probabilities2, apollo_inputs, prediction_settings = list(), modelComponent = NA) {
  prediction_settings = list()
  modelComponent = NA
  
  if (is.null(prediction_settings$modelComponent)) {
    if (exists("modelComponent")) 
      prediction_settings$modelComponent = modelComponent
    else prediction_settings$modelComponent = NA
  }
  if (is.null(prediction_settings$runs))   prediction_settings$runs = 1
  if (is.null(prediction_settings$silent))   prediction_settings$silent = FALSE
  silent = prediction_settings$silent
  if (!is.null(apollo_inputs$silent) && apollo_inputs$silent) silent = TRUE
  
  if (is.null(prediction_settings$nRep)) prediction_settings$nRep <- 100L
  
  if (is.null(prediction_settings$summary)) prediction_settings$summary <- TRUE
  
  modelComponent = prediction_settings$modelComponent
  runs = prediction_settings$runs
  apollo_inputs$nRep <- prediction_settings$nRep

  apollo_compareInputs(apollo_inputs)
  apollo_randCoeff = apollo_inputs[["apollo_randCoeff"]]
  apollo_lcPars = apollo_inputs[["apollo_lcPars"]]
  apollo_checkArguments(apollo_probabilities, apollo_randCoeff, 
                        apollo_lcPars)
  if (!silent)     apollo_print("Running predictions from model using parameter estimates...")
  
  
  ###################HERE IT REQUIRES APOLLO_PROBABILITIES
  predictions = apollo_probabilities(apollo_beta, apollo_inputs, 
                                     functionality = "prediction")
  
  predictionamount = predictions$model[,1:9] 
  
  return(predictionamount)
}

apollo_inputs = apollo_validateInputs(apollo_beta = apollo_beta2)

predictionprob2 = apollo_prediction(model, apollo_probabilities = apollo_probabilities2, apollo_inputs, prediction_settings=list(runs=1), modelComponent = NA)



t2<-Sys.time()

timeprint = t2-t1
timeprint


finalpredictionprob = predictionprob1 * membership + predictionprob2 * (1 - membership ) 

pricema = database0[,103:110]

adclickma = exp(-0.130*pricema)/(1+exp(-0.130*pricema))

pred_click = t( apply (cbind(database0[,15:22],finalpredictionprob[,2:9], adclickma), 1, add_restrict_imputevalue) )
summary(as.numeric(as.matrix(pred_click)))



########################################################################################
# combine pred_click with database1
colnames(pred_click) = c("pclick1", "pclick2", "pclick3", "pclick4", "pclick5", "pclick6", "pclick7", "pclick8")
database1 = cbind(database1, pred_click)
summary(database1)



# update datam by adding two columns: predicted n_ad and clicks
datam$pred_n_ad = n_ad_vec
datam$pred_clicks = 0

# update pred_clicks in datam
i=1

while(i<=nrow(database1)){
  pv_id_temp = database1$pv_id2[i]
  subcat_temp = database1[i, 16:23]
  subcat_temp2 = subcat_temp[database1[i, 86:93]>0]
  click_temp = database1[i, 204:211] 
  click_temp2 = click_temp[database1[i, 86:93]>0]

  # start looking for match
  index = which(datam$pv_id2==pv_id_temp & datam$subcat_id3 %in% subcat_temp2 == TRUE)

  if(length(index)>0){
    datam$pred_clicks[index] = click_temp2
  }

  if(i%%500==0){
    cat(i, fill=TRUE)
  }

  i=i+1
}



# define rprice, X_U
rprice=datab$rprice-(-0.24)

X_U=as.matrix(cbind(1, log(datam$npv14+1), log(datam$click14+1), datam$buy14, 
                    log(datam$rcc_pv), log(datam$rcc_click), log(datam$rcc_buy), 
                    0, 0, 0, datam$weekend, datam$sumpv/14))

colnames(X_U)=c("intercept","lnfreq_view","lnfreq_click","freq_buy",
                "lnrcc_view","lnrcc_click","lnrcc_buy", 
                "rprice", "ln_ad","ln_related", "weekend", "Z")

# add rootcat dummies
X_U=as.matrix(cbind(X_U, datam[, 22:25]))

X_U_mean=c(1.000,4.170,1.523,0.003,0.565,0.525,3.561,0,0,0,0.318,5.344,0.048,0.011,
           0.029,0.039,0.016,0.063,0.044,0.007,0.017,0.004,0.115,0.033,0.015)

ln_ad_mean=1.233
ln_related_mean=1.702

# mean-center
for (i in 2:ncol(X_U)){
  X_U[, i]=X_U[, i]-X_U_mean[i]
}


########################################################################################
# create a list, which records the following information for each <date, cat>
# index_allcats: the corresponding rows of pvs in datam where the set of potential cats includes <date, cat>
#                this represents all pvs in which the focal cat is potentially interested by the user
# index_focalcat: the corresponding rows in datam where date_cat_id==<date, cat>
# pv_id: all pv_ids associated with the focal <date, cat>
# user_id: all user_ids associated with each pv asscociated with <date, cat>
# n_pot_ad: number of competing ads for each cat for those pvs related to <date, cat>. same length as index
# n_pot_cat: number of cats for each pv. same length as pv_id
# cat_index: the position of the focal cat among all cats in each pv. same length as pv_id
# utility_cat: a n_user by 1 vector which record the category-specific utility in click equation
# pred_clicks (newly added): a n_user by 1 vector which records the category-specific predicted clicks

date_cat_list=list()

i=1

while(i<=N_datecat_all){
  index1=which(date_cat_id_s==i)
  pvidtemp=pv_id[index1]
  index2=which((pv_id %in% pvidtemp)==TRUE)
  useridtemp=datam$user_id3[index1]
  npotadtemp=datam$n_ad[index2]
  npotcattemp=datam$n_pot_subcat[index1]
  # To compute cat_index, we first calculate the index for each pv_id to first appear
  pvid1st=match(pvidtemp, pv_id)
  cat_index=index1-pvid1st+1
  
  # next compute utility
  if(length(index1)>1){
    Xtemp=X_U[index1, ]
  }else{
    Xtemp=matrix(X_U[index1, ], 1)
  }
  
  utility=Xtemp%*%beta_u_vec
  
  # pred_clicks
  pred_clicks_temp=datam$pred_clicks[index1]
  
  date_cat_list[[i]]=list(index_allcats=index2, index_focalcat=index1, pv_id=pvidtemp, user_id=useridtemp, 
                          n_pot_ad=npotadtemp, n_pot_cat=npotcattemp, cat_index=cat_index, utility_cat=utility,
                          pred_clicks=pred_clicks_temp)
  i=i+1
  if(i%%100==0){
    cat(i, fill=TRUE)
  }
}


########################################################################################
# a test of directly solving the FOC condition of profit maximization (previous method)
i=12

# we use CPM/CPC of APAC countries from Facebook benchmark report
ratio=0

begin=proc.time()[3]

datecatid=date_cat_id_b[i]
adid=ad_id_b[i]

index_focalcat=date_cat_list[[datecatid]]$index_focalcat
n=length(index_focalcat)
user_id=date_cat_list[[datecatid]]$user_id

# next we calculate click utility in addition to the effect of n_ad and similarity
utility_cat=date_cat_list[[datecatid]]$utility_cat

# compute ad-specific utility generated by rprice
utility_ad=beta_u_vec[8]*rprice[i]

# utility_eta is a n by ndraw_S matrix, each element is the unobserved error for the corresponding <user, catid>
utility_eta=eta_u_vec[userroot_id_s[index_focalcat]]

# similarly, we define utility_eps, which records unobserved error for each ad
utility_eps=eps_u_vec[adid]

# sum across the four terms above to obtain the sum of clicking utility (excluding the effect of n_ad and simlilarity)
utility=utility_cat+utility_ad+utility_eta+utility_eps

# beta_u is a 2 by 1 vector, records the coef. for ln_ad, ln_related
beta_u=beta_u_vec[9:10]

# define n_ad_disp
index_allcats=date_cat_list[[datecatid]]$index_allcats
n_ad_disp=n_ad_vec[index_allcats]

# define n_pot_cat, n_pot_ad, cat_index
n_pot_cat=date_cat_list[[datecatid]]$n_pot_cat
n_pot_ad=date_cat_list[[datecatid]]$n_pot_ad
# # restrict n_pot_ad to be below 50
# n_pot_ad=mapply(min, n_pot_ad, 50)
cat_index=date_cat_list[[datecatid]]$cat_index
pv_id_temp=date_cat_list[[datecatid]]$pv_id

# define X_std, which is n by 2 matrix, each row includes mean-centered ln_ad, ln_related
ln_ad_mean=1.233
ln_related_mean=1.702
X_std=cbind((log(n_ad_vec[index_focalcat]+1)-ln_ad_mean), (log(Sim_vec[pv_id_temp]+1)-ln_related_mean))

# define pconsider
# pconsider=datam$P_consider_vec[index_focalcat]
pconsider=rep(1, length(index_focalcat))

# define bid, qscore and u
b=datab$bid[i]
QS=datab$qs_post[i]
utemp=ranmat[pv_id_temp, ]

v=3

out=EPROFIT(n_ad_disp, n_pot_cat, n_pot_ad, cat_index, pconsider, utility, X_std, beta_u,
            b, QS, v, ratio, a, var, tao, utemp, rule30)

# infer value-per-click from observed bid
out=uniroot(FOC, c(0.1*b, 10*b), tol=1e-6, maxiter=1e3, b, n_ad_disp, n_pot_cat, n_pot_ad, cat_index, 
            pconsider, utility, X_std, beta_u, QS, ratio, a, var, tao, utemp, rule30)
v=out$root
v
b

end=proc.time()[3] 
cat("Time used:", end-begin, fill=TRUE)

# check SOC
hessian(EPROFIT,b,method="Richardson", v=v, n_ad_disp=n_ad_disp, n_pot_cat=n_pot_cat, n_pot_ad=n_pot_ad, cat_index=cat_index, 
        pconsider=pconsider, utility=utility, X_std=X_std, beta_u=beta_u,
        QS=QS, ratio=ratio, a=a, var=var, tao=tao, u=utemp, rule=rule30)



########################################################################################
# a test of directly solving the FOC condition of profit maximization (new method with new ad-click model)
i=12

# we use CPM/CPC of APAC countries from Facebook benchmark report
ratio=0

begin=proc.time()[3]

datecatid=date_cat_id_b[i]
adid=ad_id_b[i]

index_focalcat=date_cat_list[[datecatid]]$index_focalcat
user_id=date_cat_list[[datecatid]]$user_id

# next we calculate expected clicks
pclick=date_cat_list[[datecatid]]$pred_clicks

# define n_ad_disp
index_allcats=date_cat_list[[datecatid]]$index_allcats
n_ad_disp=n_ad_vec[index_allcats]

# define n_pot_cat, n_pot_ad, cat_index
n_pot_cat=date_cat_list[[datecatid]]$n_pot_cat
n_pot_ad=date_cat_list[[datecatid]]$n_pot_ad
# # restrict n_pot_ad to be below 50
# n_pot_ad=mapply(min, n_pot_ad, 50)
cat_index=date_cat_list[[datecatid]]$cat_index
pv_id_temp=date_cat_list[[datecatid]]$pv_id

# define bid, qscore and u
b=datab$bid[i]
QS=datab$qs_post[i]
utemp=ranmat[pv_id_temp, ]

v=3

out=EPROFIT_V2(n_ad_disp, n_pot_cat, n_pot_ad, cat_index, pclick,
            b, QS, v, a, var, tao, utemp, rule30)

# infer value-per-click from observed bid
out=uniroot(FOC_V2, c(0.1*b, 10*b), tol=1e-6, maxiter=1e3, b, n_ad_disp, n_pot_cat, n_pot_ad, cat_index, 
            pclick, QS, a, var, tao, utemp, rule30)
v=out$root
v
b

end=proc.time()[3] 
cat("Time used:", end-begin, fill=TRUE)

# check SOC
hessian(EPROFIT_V2,b,method="Richardson", v=v, n_ad_disp=n_ad_disp, n_pot_cat=n_pot_cat, n_pot_ad=n_pot_ad, cat_index=cat_index, 
        pclick=pclick, QS=QS, a=a, var=var, tao=tao, u=utemp, rule=rule30)


########################################################################################
# define a function which solve the vpc for each row in datab (based on the new model)
VPCSOLVE_V2=function(i, date_cat_list, date_cat_id_b, ad_id_b,
                     pclick, n_ad_vec, ranmat, datab, a, var, tao, rule){
  # output is a vector which records i, date_subcat_id, adid2, bid, vpc, flag
  # flag=1 means success
  # flag=2 means no impression for associated subcat
  # flag=3 means initial values are out of bounds
  datecatid=date_cat_id_b[i]
  adid=ad_id_b[i]
  
  index_focalcat=date_cat_list[[datecatid]]$index_focalcat
  user_id=date_cat_list[[datecatid]]$user_id
  
  # next we calculate expected clicks
  pclick=date_cat_list[[datecatid]]$pred_clicks
  
  # define n_ad_disp
  index_allcats=date_cat_list[[datecatid]]$index_allcats
  n_ad_disp=n_ad_vec[index_allcats]
  
  # define n_pot_cat, n_pot_ad, cat_index
  n_pot_cat=date_cat_list[[datecatid]]$n_pot_cat
  n_pot_ad=date_cat_list[[datecatid]]$n_pot_ad
  # # restrict n_pot_ad to be below 50
  # n_pot_ad=mapply(min, n_pot_ad, 50)
  cat_index=date_cat_list[[datecatid]]$cat_index
  pv_id_temp=date_cat_list[[datecatid]]$pv_id
  
  # define bid, qscore and u
  b=datab$bid[i]
  QS=datab$qs_post[i]
  utemp=ranmat[pv_id_temp, ]
  
  # meanwhile, we check if the cat of focal ad does not get any exposure
  n_ad_check=max(n_ad_vec[index_focalcat])
  
  # if n_ad_check==0, then we deal with this case later
  if(n_ad_check==0){
    v=0
    flag=2
  }else{
    # infer value-per-click from observed bid
    upperbound=b*4
    lowerbound=b*0.25
    temp1=FOC_V2(lowerbound, b, n_ad_disp, n_pot_cat, n_pot_ad, cat_index, pclick, 
              QS, a, var, tao, utemp, rule)
    temp2=FOC_V2(upperbound, b, n_ad_disp, n_pot_cat, n_pot_ad, cat_index, pclick, 
              QS, a, var, tao, utemp, rule)
    while(temp1*temp2>0 & upperbound<100){
      lowerbound=lowerbound*0.5
      upperbound=upperbound*2
      temp1=FOC_V2(lowerbound, b, n_ad_disp, n_pot_cat, n_pot_ad, cat_index, pclick, 
                QS, a, var, tao, utemp, rule)
      temp2=FOC_V2(upperbound, b, n_ad_disp, n_pot_cat, n_pot_ad, cat_index, pclick, 
                QS, a, var, tao, utemp, rule)
    }
    if(temp1*temp2<0){
      out=uniroot(FOC_V2, c(lowerbound, upperbound), tol=1e-3, maxiter=1e3, b, n_ad_disp, n_pot_cat, n_pot_ad, cat_index, 
                  pclick, QS, a, var, tao, utemp, rule)
      v=out$root
      flag=1
    }else{
      v=0
      flag=3
    }  
  }
  return(c(i, datecatid, adid, b, v, flag))
}




########################################################################################
########################################################################################
# Let's first infer vpc for the first N rows in datab
N=10000

cl = makeCluster(7, outfile="")
registerDoParallel(cl)

begin=proc.time()[3]


out1=foreach(i=1:N, .packages=c('numDeriv', 'fastGHQuad'), .verbose = FALSE) %dopar% {
  VPCSOLVE_V2(i, date_cat_list, date_cat_id_b, ad_id_b, 
              pclick, n_ad_vec, ranmat, datab, a, var, tao, rule30)
}
stopCluster(cl)
end=proc.time()[3] 
cat("Time used:", end-begin, fill=TRUE)

# save out1
write.csv(t(sapply(out1, rbind)), "vpc_ratio=0_1to10000.csv",
          row.names=FALSE)


########################################################################################
# do more from 10001 to end
N=10000

cl = makeCluster(7, outfile="")
registerDoParallel(cl)

begin=proc.time()[3]

out2=foreach(i=(1+N):nrow(datab), .packages=c('numDeriv', 'fastGHQuad'), .verbose = FALSE) %dopar% {
  VPCSOLVE_V2(i, date_cat_list, date_cat_id_b, ad_id_b, 
              pclick, n_ad_vec, ranmat, datab, a, var, tao, rule30)
}
stopCluster(cl)
end=proc.time()[3] 
cat("Time used:", end-begin, fill=TRUE)

# save out2
write.csv(t(sapply(out2, rbind)), "vpc_ratio=0_10001toEnd.csv",
          row.names=FALSE)







  
  



