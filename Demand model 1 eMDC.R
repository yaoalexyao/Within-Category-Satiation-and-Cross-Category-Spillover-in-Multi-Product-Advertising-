########################################################################################
# This file does estimate the eMDC model


library(stargazer)



# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #

### Clear memory and initialise
rm(list = ls())
library(apollo)
apollo_initialise()



apollo_control = list(
  modelName  ="Demand eMDC model with latent class",   
  modelDescr ="Demand eMDC model with latent class",  
  indivID    ="user_id3",
  nCores     = 7,
  workInLogs = F  
)



# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

### Load data from within the Apollo package


database = read.csv("cleaned data.csv",header=TRUE)


# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #


rdn = rnorm(100)


### Parameters starting values 
apollo_beta = c(
  gamma_v1_a      = abs(rdn[1]), 
  gamma_v2_a      = abs(rdn[2]), 
  gamma_v3_a      = abs(rdn[3]), 
  gamma_v4_a      = abs(rdn[4]), 
  gamma_v5_a      = abs(rdn[5]), 

  delta_v1_a      = -abs(rdn[6]), 
  delta_v2_a         = -abs(rdn[7]), 
  delta_v3_a       = -abs(rdn[8]), 
  delta_v4_a     = -abs(rdn[9]), 
  delta_v5_a = -abs(rdn[10]), 

  
  gamma_v1_b      = abs(rdn[11]), 
  gamma_v2_b      = abs(rdn[12]), 
  gamma_v3_b      = abs(rdn[13]), 
  gamma_v4_b      = abs(rdn[14]), 
  gamma_v5_b      = abs(rdn[15]), 

  delta_v1_b      = -abs(rdn[16]), 
  delta_v2_b         = -abs(rdn[17]), 
  delta_v3_b       = -abs(rdn[18]), 
  delta_v4_b     = -abs(rdn[19]), 
  delta_v5_b = -abs(rdn[20]), 

  bweekend_a = rdn[21], 
  bnpv_a = rdn[22], 
  bnclick_a = rdn[23], 
  bnbuy_a = rdn[24], 
  nrccpv_a=rdn[25], 
  nrccclick_a=rdn[26], 
  nrccbuy_a=rdn[27], 
  badstockpsai_a =rdn[28], 
  bresidual_a = rdn[29], 
  bvariety_a = rdn[30], 

  d12_a= rdn[31],
  d13_a= rdn[32],
  d14_a= rdn[33], 
  d15_a= rdn[34],
  d23_a= rdn[35],
  d24_a= rdn[36],
  d25_a= rdn[37],
  d34_a= rdn[38],
  d35_a= rdn[39],
  d45_a= rdn[40],
  
  bweekend_b = rdn[41], 
  bnpv_b = rdn[42], 
  bnclick_b = rdn[43], 
  bnbuy_b = rdn[44], 
  nrccpv_b=rdn[45], 
  nrccclick_b=rdn[46], 
  nrccbuy_b=rdn[47], 
  badstockpsai_b =rdn[48], 
  bresidual_b = rdn[49], 
  bvariety_b = rdn[50], 
  
  d12_b= rdn[51],
  d13_b= rdn[52],
  d14_b= rdn[53], 
  d15_b= rdn[54],
  d23_b= rdn[55],
  d24_b= rdn[56],
  d25_b= rdn[57],
  d34_b= rdn[58],
  d35_b= rdn[59],
  d45_b= rdn[60],
  
  #parameter for segment assignment
  delta_a         = rdn[61],
  gamma_N_pot_subcat_full_a = rdn[62],
  gamma_Z_pvperday_a  = rdn[63]
  
)


apollo_fixed = c( )


apollo_lcPars=function(apollo_beta, apollo_inputs){
  lcpars = list()
  lcpars[["gamma_v1"]] = list(exp(gamma_v1_a), exp(gamma_v1_b) )
  lcpars[["gamma_v2"]] = list(exp(gamma_v2_a), exp(gamma_v2_b)   )
  lcpars[["gamma_v3"]] = list(exp(gamma_v3_a), exp(gamma_v3_b)   )
  lcpars[["gamma_v4"]] = list(exp(gamma_v4_a), exp(gamma_v4_b)   )
  lcpars[["gamma_v5"]] = list(exp(gamma_v5_a), exp(gamma_v5_b)   )
  
  
  lcpars[["delta_v1"]] = list(delta_v1_a, delta_v1_b)
  lcpars[["delta_v2"]] = list(delta_v2_a, delta_v2_b)
  lcpars[["delta_v3"]] = list(delta_v3_a, delta_v3_b)
  lcpars[["delta_v4"]] = list(delta_v4_a, delta_v4_b)
  lcpars[["delta_v5"]] = list(delta_v5_a, delta_v5_b)
  
  lcpars[["bweekend"]] = list(bweekend_a ,  bweekend_b )
  lcpars[["bnpv"]] = list( bnpv_a ,   bnpv_b )
  lcpars[["bnclick"]] = list( bnclick_a ,   bnclick_b )
  lcpars[["bnbuy"]] = list(bnbuy_a ,  bnbuy_b )
  lcpars[["nrccpv"]] = list(nrccpv_a ,  nrccpv_b )
  lcpars[["nrccclick"]] = list(nrccclick_a ,  nrccclick_b )
  lcpars[["nrccbuy"]] = list(nrccbuy_a ,  nrccbuy_b )
  lcpars[["badstockpsai"]] = list(badstockpsai_a ,  badstockpsai_b )
  lcpars[["bresidual"]] = list(bresidual_a, bresidual_b)
  lcpars[["bvariety"]] = list(bvariety_a, bvariety_b)
  
  lcpars[["d12"]] = list(d12_a, d12_b)
  lcpars[["d13"]] = list(d13_a, d13_b)
  lcpars[["d14"]] = list(d14_a, d14_b)
  lcpars[["d15"]] = list(d15_a, d15_b)
  lcpars[["d23"]] = list(d23_a, d23_b)
  lcpars[["d24"]] = list(d24_a, d24_b)
  lcpars[["d25"]] = list(d25_a, d25_b)
  lcpars[["d34"]] = list(d34_a, d34_b)
  lcpars[["d35"]] = list(d35_a, d35_b)
  lcpars[["d45"]] = list(d45_a, d45_b)
  
  #This part is the probability of each class. 
  V=list()
  V[["class_a"]] = delta_a + gamma_N_pot_subcat_full_a*N_pot_subcat_full + gamma_Z_pvperday_a*Z_pvperday
  V[["class_b"]] = 0
  
  #This part is the probability of each class.   
  mnl_settings = list(
    alternatives = c(class_a=1, class_b=2), 
    avail        = 1, 
    choiceVar    = NA, 
    V            = V
  )
  
  lcpars[["pi_values"]] = apollo_mnl(mnl_settings, functionality="raw") #this part returns the probability. It use "raw" to ensure that the probabilities are returned for all alternatives.
  lcpars[["pi_values"]] = apollo_firstRow(lcpars[["pi_values"]], apollo_inputs)
  
  return(lcpars)
}

# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #
apollo_probabilities=function(apollo_beta ){
  
  apollo_inputs = apollo_inputs
  functionality="estimate"
  
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
               v8 = availnew8 #,
  )
  
  ## Define continuous consumption for individual alternatives
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
                        utilityOutside   = 0, 
                        budget           = budget,
                        sigma            = 0.99, 
                        cost             = cost)
  
  # 
  ### Loop over classes
  s=1
  while(s<=2){
    
    ### ### Compute class-specific utilities
    V = list()
    
    V[["v1"    ]] =   delta_v1[[s]]*(rootcat1==1) +  delta_v2[[s]]*(rootcat1==2) +  delta_v3[[s]]*(rootcat1==3) +  
      delta_v4[[s]]*(rootcat1==4) +  delta_v5[[s]]*(rootcat1==5) + 
      bweekend[[s]]*weekend + bnpv[[s]]*npvnew1 + bnclick[[s]]*nclicknew1 + bnbuy[[s]]*nbuynew1 + 
      nrccpv[[s]]*rccpvnew1 + nrccclick[[s]]*rccclicknew1 + nrccbuy[[s]]*rccbuynew1 + badstockpsai[[s]]*ad11 + bresidual[[s]]*res1 + bvariety[[s]]*Lag_Distinct
    
    V[["v2"    ]] = delta_v1[[s]]*(rootcat2==1) +  delta_v2[[s]]*(rootcat2==2) +  delta_v3[[s]]*(rootcat2==3) +  
      delta_v4[[s]]*(rootcat2==4) +  delta_v5[[s]]*(rootcat2==5) + 
      bweekend[[s]]*weekend + bnpv[[s]]*npvnew2 + bnclick[[s]]*nclicknew2 + bnbuy[[s]]*nbuynew2 + 
      nrccpv[[s]]*rccpvnew2 + nrccclick[[s]]*rccclicknew2 + nrccbuy[[s]]*rccbuynew2 + badstockpsai[[s]]*ad21 + bresidual[[s]]*res2 + bvariety[[s]]*Lag_Distinct
    
    V[["v3"  ]] = delta_v1[[s]]*(rootcat3==1) +  delta_v2[[s]]*(rootcat3==2) +  delta_v3[[s]]*(rootcat3==3) +  
      delta_v4[[s]]*(rootcat3==4) +  delta_v5[[s]]*(rootcat3==5) + 
      bweekend[[s]]*weekend + bnpv[[s]]*npvnew3 + bnclick[[s]]*nclicknew3 + bnbuy[[s]]*nbuynew3 + 
      nrccpv[[s]]*rccpvnew3 + nrccclick[[s]]*rccclicknew3 + nrccbuy[[s]]*rccbuynew3 + badstockpsai[[s]]*ad31 + bresidual[[s]]*res3 + bvariety[[s]]*Lag_Distinct
    
    V[["v4"]] = delta_v1[[s]]*(rootcat4==1) +  delta_v2[[s]]*(rootcat4==2) +  delta_v3[[s]]*(rootcat4==3) +  
      delta_v4[[s]]*(rootcat4==4) +  delta_v5[[s]]*(rootcat4==5) + 
      bweekend[[s]]*weekend + bnpv[[s]]*npvnew4 + bnclick[[s]]*nclicknew4 + bnbuy[[s]]*nbuynew4 + 
      nrccpv[[s]]*rccpvnew4 + nrccclick[[s]]*rccclicknew4 + nrccbuy[[s]]*rccbuynew4 + badstockpsai[[s]]*ad41 + bresidual[[s]]*res4 + bvariety[[s]]*Lag_Distinct
    
    V[["v5"]] = delta_v1[[s]]*(rootcat5==1) +  delta_v2[[s]]*(rootcat5==2) +  delta_v3[[s]]*(rootcat5==3) +  
      delta_v4[[s]]*(rootcat5==4) +  delta_v5[[s]]*(rootcat5==5) + 
      bweekend[[s]]*weekend + bnpv[[s]]*npvnew5 + bnclick[[s]]*nclicknew5 + bnbuy[[s]]*nbuynew5 + 
      nrccpv[[s]]*rccpvnew5 + nrccclick[[s]]*rccclicknew5 + nrccbuy[[s]]*rccbuynew5 + badstockpsai[[s]]*ad51 + bresidual[[s]]*res5 + bvariety[[s]]*Lag_Distinct
    
    V[["v6"  ]] = delta_v1[[s]]*(rootcat6==1) +  delta_v2[[s]]*(rootcat6==2) +  delta_v3[[s]]*(rootcat6==3) +  
      delta_v4[[s]]*(rootcat6==4) +  delta_v5[[s]]*(rootcat6==5) + 
      bweekend[[s]]*weekend + bnpv[[s]]*npvnew6 + bnclick[[s]]*nclicknew6 + bnbuy[[s]]*nbuynew6 + 
      nrccpv[[s]]*rccpvnew6 + nrccclick[[s]]*rccclicknew6 + nrccbuy[[s]]*rccbuynew6 + badstockpsai[[s]]*ad61 + bresidual[[s]]*res6 + bvariety[[s]]*Lag_Distinct
    
    V[["v7"]] = delta_v1[[s]]*(rootcat7==1) +  delta_v2[[s]]*(rootcat7==2) +  delta_v3[[s]]*(rootcat7==3) +  
      delta_v4[[s]]*(rootcat7==4) +  delta_v5[[s]]*(rootcat7==5) + 
      bweekend[[s]]*weekend + bnpv[[s]]*npvnew7 + bnclick[[s]]*nclicknew7 + bnbuy[[s]]*nbuynew7 + 
      nrccpv[[s]]*rccpvnew7 + nrccclick[[s]]*rccclicknew7 + nrccbuy[[s]]*rccbuynew7 + badstockpsai[[s]]*ad71 + bresidual[[s]]*res7 + bvariety[[s]]*Lag_Distinct
    
    V[["v8"]] = delta_v1[[s]]*(rootcat8==1) +  delta_v2[[s]]*(rootcat8==2) +  delta_v3[[s]]*(rootcat8==3) +  
      delta_v4[[s]]*(rootcat8==4) +  delta_v5[[s]]*(rootcat8==5) + 
      bweekend[[s]]*weekend + bnpv[[s]]*npvnew8 + bnclick[[s]]*nclicknew8 + bnbuy[[s]]*nbuynew8 + 
      nrccpv[[s]]*rccpvnew8 + nrccclick[[s]]*rccclicknew8 + nrccbuy[[s]]*rccbuynew8 + badstockpsai[[s]]*ad81 + bresidual[[s]]*res8 + bvariety[[s]]*Lag_Distinct

    
    
    ### Define gamma parameters
    gamma = list(v1      = gamma_v1[[s]]*(rootcat1==1) +  gamma_v2[[s]]*(rootcat1==2) +  gamma_v3[[s]]*(rootcat1==3) +  
                   gamma_v4[[s]]*(rootcat1==4) +  gamma_v5[[s]]*(rootcat1==5) , 
                 
                 v2      = gamma_v1[[s]]*(rootcat2==1) +  gamma_v2[[s]]*(rootcat2==2) +  gamma_v3[[s]]*(rootcat2==3) +  
                   gamma_v4[[s]]*(rootcat2==4) +  gamma_v5[[s]]*(rootcat2==5) ,
                 
                 
                 v3      = gamma_v1[[s]]*(rootcat3==1) +  gamma_v2[[s]]*(rootcat3==2) +  gamma_v3[[s]]*(rootcat3==3) +  
                   gamma_v4[[s]]*(rootcat3==4) +  gamma_v5[[s]]*(rootcat3==5) ,
                 
                 v4      = gamma_v1[[s]]*(rootcat4==1) +  gamma_v2[[s]]*(rootcat4==2) +  gamma_v3[[s]]*(rootcat4==3) +  
                   gamma_v4[[s]]*(rootcat4==4) +  gamma_v5[[s]]*(rootcat4==5) , 
                 
                 v5      = gamma_v1[[s]]*(rootcat5==1) +  gamma_v2[[s]]*(rootcat5==2) +  gamma_v3[[s]]*(rootcat5==3) +  
                   gamma_v4[[s]]*(rootcat5==4) +  gamma_v5[[s]]*(rootcat5==5) , 
                 
                 v6      = gamma_v1[[s]]*(rootcat6==1) +  gamma_v2[[s]]*(rootcat6==2) +  gamma_v3[[s]]*(rootcat6==3) +  
                   gamma_v4[[s]]*(rootcat6==4) +  gamma_v5[[s]]*(rootcat6==5) , 
                 
                 v7      = gamma_v1[[s]]*(rootcat7==1) +  gamma_v2[[s]]*(rootcat7==2) +  gamma_v3[[s]]*(rootcat7==3) +  
                   gamma_v4[[s]]*(rootcat7==4) +  gamma_v5[[s]]*(rootcat7==5) , 
                 
                 v8      = gamma_v1[[s]]*(rootcat8==1) +  gamma_v2[[s]]*(rootcat8==2) +  gamma_v3[[s]]*(rootcat8==3) +  
                   gamma_v4[[s]]*(rootcat8==4) +  gamma_v5[[s]]*(rootcat8==5) # , 

    )
    
    delta = list( list(0,0,0,0,0,0,0 ,0 ) , 
                  list(   d12[[s]] * ( (rootcat2 == 1) * (rootcat1 == 2) + (rootcat2 == 2) * (rootcat1 == 1) )    +
                            d13[[s]] * ( (rootcat2 == 1) * (rootcat1 == 3) + (rootcat2 == 3) * (rootcat1 == 1) )  +
                            d14[[s]] * ( (rootcat2 == 1) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 1) )  +
                            d15[[s]] * ( (rootcat2 == 1) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 1) )  +
                            d23[[s]] * ( (rootcat2 == 2) * (rootcat1 == 3) + (rootcat2 == 3) * (rootcat1 == 2) )  +
                            d24[[s]] * ( (rootcat2 == 2) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 2) )  +
                            d25[[s]] * ( (rootcat2 == 2) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 2) )  +
                            d34[[s]] * ( (rootcat2 == 3) * (rootcat1 == 4) + (rootcat2 == 4) * (rootcat1 == 3) )  +
                            d35[[s]] * ( (rootcat2 == 3) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 3) )  +
                            d45[[s]] * ( (rootcat2 == 4) * (rootcat1 == 5) + (rootcat2 == 5) * (rootcat1 == 4) )  ,
                          0,0,0,0,0,0,0 ) ,
                  
                  list(   d12[[s]] * ( (rootcat3 == 1) * (rootcat1 == 2) + (rootcat3 == 2) * (rootcat1 == 1) )    +
                            d13[[s]] * ( (rootcat3 == 1) * (rootcat1 == 3) + (rootcat3 == 3) * (rootcat1 == 1) )  +
                            d14[[s]] * ( (rootcat3 == 1) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 1) )  +
                            d15[[s]] * ( (rootcat3 == 1) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 1) )  +
                            d23[[s]] * ( (rootcat3 == 2) * (rootcat1 == 3) + (rootcat3 == 3) * (rootcat1 == 2) )  +
                            d24[[s]] * ( (rootcat3 == 2) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 2) )  +
                            d25[[s]] * ( (rootcat3 == 2) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 2) )  +
                            d34[[s]] * ( (rootcat3 == 3) * (rootcat1 == 4) + (rootcat3 == 4) * (rootcat1 == 3) )  +
                            d35[[s]] * ( (rootcat3 == 3) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 3) )  +
                            d45[[s]] * ( (rootcat3 == 4) * (rootcat1 == 5) + (rootcat3 == 5) * (rootcat1 == 4) )  ,
                          
                          d12[[s]] * ( (rootcat3 == 1) * (rootcat2 == 2) + (rootcat3 == 2) * (rootcat2 == 1) )    +
                            d13[[s]] * ( (rootcat3 == 1) * (rootcat2 == 3) + (rootcat3 == 3) * (rootcat2 == 1) )  +
                            d14[[s]] * ( (rootcat3 == 1) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 1) )  +
                            d15[[s]] * ( (rootcat3 == 1) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 1) )  +
                            d23[[s]] * ( (rootcat3 == 2) * (rootcat2 == 3) + (rootcat3 == 3) * (rootcat2 == 2) )  +
                            d24[[s]] * ( (rootcat3 == 2) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 2) )  +
                            d25[[s]] * ( (rootcat3 == 2) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 2) )  +
                            d34[[s]] * ( (rootcat3 == 3) * (rootcat2 == 4) + (rootcat3 == 4) * (rootcat2 == 3) )  +
                            d35[[s]] * ( (rootcat3 == 3) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 3) )  +
                            d45[[s]] * ( (rootcat3 == 4) * (rootcat2 == 5) + (rootcat3 == 5) * (rootcat2 == 4) )  ,
                          0,0,0,0,0,0 ) , 
                  
                  list(     d12[[s]] * ( (rootcat4 == 1) * (rootcat1 == 2) + (rootcat4 == 2) * (rootcat1 == 1) )    +
                              d13[[s]] * ( (rootcat4 == 1) * (rootcat1 == 3) + (rootcat4 == 3) * (rootcat1 == 1) )  +
                              d14[[s]] * ( (rootcat4 == 1) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 1) )  +
                              d15[[s]] * ( (rootcat4 == 1) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 1) )  +
                              d23[[s]] * ( (rootcat4 == 2) * (rootcat1 == 3) + (rootcat4 == 3) * (rootcat1 == 2) )  +
                              d24[[s]] * ( (rootcat4 == 2) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 2) )  +
                              d25[[s]] * ( (rootcat4 == 2) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 2) )  +
                              d34[[s]] * ( (rootcat4 == 3) * (rootcat1 == 4) + (rootcat4 == 4) * (rootcat1 == 3) )  +
                              d35[[s]] * ( (rootcat4 == 3) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 3) )  +
                              d45[[s]] * ( (rootcat4 == 4) * (rootcat1 == 5) + (rootcat4 == 5) * (rootcat1 == 4) )  ,
                            
                            d12[[s]] * ( (rootcat4 == 1) * (rootcat2 == 2) + (rootcat4 == 2) * (rootcat2 == 1) )    +
                              d13[[s]] * ( (rootcat4 == 1) * (rootcat2 == 3) + (rootcat4 == 3) * (rootcat2 == 1) )  +
                              d14[[s]] * ( (rootcat4 == 1) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 1) )  +
                              d15[[s]] * ( (rootcat4 == 1) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 1) )  +
                              d23[[s]] * ( (rootcat4 == 2) * (rootcat2 == 3) + (rootcat4 == 3) * (rootcat2 == 2) )  +
                              d24[[s]] * ( (rootcat4 == 2) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 2) )  +
                              d25[[s]] * ( (rootcat4 == 2) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 2) )  +
                              d34[[s]] * ( (rootcat4 == 3) * (rootcat2 == 4) + (rootcat4 == 4) * (rootcat2 == 3) )  +
                              d35[[s]] * ( (rootcat4 == 3) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 3) )  +
                              d45[[s]] * ( (rootcat4 == 4) * (rootcat2 == 5) + (rootcat4 == 5) * (rootcat2 == 4) )  ,
                            
                            d12[[s]] * ( (rootcat4 == 1) * (rootcat3 == 2) + (rootcat4 == 2) * (rootcat3 == 1) )    +
                              d13[[s]] * ( (rootcat4 == 1) * (rootcat3 == 3) + (rootcat4 == 3) * (rootcat3 == 1) )  +
                              d14[[s]] * ( (rootcat4 == 1) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 1) )  +
                              d15[[s]] * ( (rootcat4 == 1) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 1) )  +
                              d23[[s]] * ( (rootcat4 == 2) * (rootcat3 == 3) + (rootcat4 == 3) * (rootcat3 == 2) )  +
                              d24[[s]] * ( (rootcat4 == 2) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 2) )  +
                              d25[[s]] * ( (rootcat4 == 2) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 2) )  +
                              d34[[s]] * ( (rootcat4 == 3) * (rootcat3 == 4) + (rootcat4 == 4) * (rootcat3 == 3) )  +
                              d35[[s]] * ( (rootcat4 == 3) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 3) )  +
                              d45[[s]] * ( (rootcat4 == 4) * (rootcat3 == 5) + (rootcat4 == 5) * (rootcat3 == 4) )  ,
                            0,0,0,0,0 ) ,
                  
                  list(       d12[[s]] * ( (rootcat5 == 1) * (rootcat1 == 2) + (rootcat5 == 2) * (rootcat1 == 1) )    +
                                d13[[s]] * ( (rootcat5 == 1) * (rootcat1 == 3) + (rootcat5 == 3) * (rootcat1 == 1) )  +
                                d14[[s]] * ( (rootcat5 == 1) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 1) )  +
                                d15[[s]] * ( (rootcat5 == 1) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 1) )  +
                                d23[[s]] * ( (rootcat5 == 2) * (rootcat1 == 3) + (rootcat5 == 3) * (rootcat1 == 2) )  +
                                d24[[s]] * ( (rootcat5 == 2) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 2) )  +
                                d25[[s]] * ( (rootcat5 == 2) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 2) )  +
                                d34[[s]] * ( (rootcat5 == 3) * (rootcat1 == 4) + (rootcat5 == 4) * (rootcat1 == 3) )  +
                                d35[[s]] * ( (rootcat5 == 3) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 3) )  +
                                d45[[s]] * ( (rootcat5 == 4) * (rootcat1 == 5) + (rootcat5 == 5) * (rootcat1 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat5 == 1) * (rootcat2 == 2) + (rootcat5 == 2) * (rootcat2 == 1) )    +
                                d13[[s]] * ( (rootcat5 == 1) * (rootcat2 == 3) + (rootcat5 == 3) * (rootcat2 == 1) )  +
                                d14[[s]] * ( (rootcat5 == 1) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 1) )  +
                                d15[[s]] * ( (rootcat5 == 1) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 1) )  +
                                d23[[s]] * ( (rootcat5 == 2) * (rootcat2 == 3) + (rootcat5 == 3) * (rootcat2 == 2) )  +
                                d24[[s]] * ( (rootcat5 == 2) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 2) )  +
                                d25[[s]] * ( (rootcat5 == 2) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 2) )  +
                                d34[[s]] * ( (rootcat5 == 3) * (rootcat2 == 4) + (rootcat5 == 4) * (rootcat2 == 3) )  +
                                d35[[s]] * ( (rootcat5 == 3) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 3) )  +
                                d45[[s]] * ( (rootcat5 == 4) * (rootcat2 == 5) + (rootcat5 == 5) * (rootcat2 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat5 == 1) * (rootcat3 == 2) + (rootcat5 == 2) * (rootcat3 == 1) )    +
                                d13[[s]] * ( (rootcat5 == 1) * (rootcat3 == 3) + (rootcat5 == 3) * (rootcat3 == 1) )  +
                                d14[[s]] * ( (rootcat5 == 1) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 1) )  +
                                d15[[s]] * ( (rootcat5 == 1) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 1) )  +
                                d23[[s]] * ( (rootcat5 == 2) * (rootcat3 == 3) + (rootcat5 == 3) * (rootcat3 == 2) )  +
                                d24[[s]] * ( (rootcat5 == 2) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 2) )  +
                                d25[[s]] * ( (rootcat5 == 2) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 2) )  +
                                d34[[s]] * ( (rootcat5 == 3) * (rootcat3 == 4) + (rootcat5 == 4) * (rootcat3 == 3) )  +
                                d35[[s]] * ( (rootcat5 == 3) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 3) )  +
                                d45[[s]] * ( (rootcat5 == 4) * (rootcat3 == 5) + (rootcat5 == 5) * (rootcat3 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat5 == 1) * (rootcat4 == 2) + (rootcat5 == 2) * (rootcat4 == 1) )    +
                                d13[[s]] * ( (rootcat5 == 1) * (rootcat4 == 3) + (rootcat5 == 3) * (rootcat4 == 1) )  +
                                d14[[s]] * ( (rootcat5 == 1) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 1) )  +
                                d15[[s]] * ( (rootcat5 == 1) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 1) )  +
                                d23[[s]] * ( (rootcat5 == 2) * (rootcat4 == 3) + (rootcat5 == 3) * (rootcat4 == 2) )  +
                                d24[[s]] * ( (rootcat5 == 2) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 2) )  +
                                d25[[s]] * ( (rootcat5 == 2) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 2) )  +
                                d34[[s]] * ( (rootcat5 == 3) * (rootcat4 == 4) + (rootcat5 == 4) * (rootcat4 == 3) )  +
                                d35[[s]] * ( (rootcat5 == 3) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 3) )  +
                                d45[[s]] * ( (rootcat5 == 4) * (rootcat4 == 5) + (rootcat5 == 5) * (rootcat4 == 4) )  ,
                              
                              0,0,0,0) ,
                  
                  list(       d12[[s]] * ( (rootcat6 == 1) * (rootcat1 == 2) + (rootcat6 == 2) * (rootcat1 == 1) )    +
                                d13[[s]] * ( (rootcat6 == 1) * (rootcat1 == 3) + (rootcat6 == 3) * (rootcat1 == 1) )  +
                                d14[[s]] * ( (rootcat6 == 1) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 1) )  +
                                d15[[s]] * ( (rootcat6 == 1) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 1) )  +
                                d23[[s]] * ( (rootcat6 == 2) * (rootcat1 == 3) + (rootcat6 == 3) * (rootcat1 == 2) )  +
                                d24[[s]] * ( (rootcat6 == 2) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 2) )  +
                                d25[[s]] * ( (rootcat6 == 2) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 2) )  +
                                d34[[s]] * ( (rootcat6 == 3) * (rootcat1 == 4) + (rootcat6 == 4) * (rootcat1 == 3) )  +
                                d35[[s]] * ( (rootcat6 == 3) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 3) )  +
                                d45[[s]] * ( (rootcat6 == 4) * (rootcat1 == 5) + (rootcat6 == 5) * (rootcat1 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat6 == 1) * (rootcat2 == 2) + (rootcat6 == 2) * (rootcat2 == 1) )    +
                                d13[[s]] * ( (rootcat6 == 1) * (rootcat2 == 3) + (rootcat6 == 3) * (rootcat2 == 1) )  +
                                d14[[s]] * ( (rootcat6 == 1) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 1) )  +
                                d15[[s]] * ( (rootcat6 == 1) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 1) )  +
                                d23[[s]] * ( (rootcat6 == 2) * (rootcat2 == 3) + (rootcat6 == 3) * (rootcat2 == 2) )  +
                                d24[[s]] * ( (rootcat6 == 2) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 2) )  +
                                d25[[s]] * ( (rootcat6 == 2) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 2) )  +
                                d34[[s]] * ( (rootcat6 == 3) * (rootcat2 == 4) + (rootcat6 == 4) * (rootcat2 == 3) )  +
                                d35[[s]] * ( (rootcat6 == 3) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 3) )  +
                                d45[[s]] * ( (rootcat6 == 4) * (rootcat2 == 5) + (rootcat6 == 5) * (rootcat2 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat6 == 1) * (rootcat3 == 2) + (rootcat6 == 2) * (rootcat3 == 1) )    +
                                d13[[s]] * ( (rootcat6 == 1) * (rootcat3 == 3) + (rootcat6 == 3) * (rootcat3 == 1) )  +
                                d14[[s]] * ( (rootcat6 == 1) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 1) )  +
                                d15[[s]] * ( (rootcat6 == 1) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 1) )  +
                                d23[[s]] * ( (rootcat6 == 2) * (rootcat3 == 3) + (rootcat6 == 3) * (rootcat3 == 2) )  +
                                d24[[s]] * ( (rootcat6 == 2) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 2) )  +
                                d25[[s]] * ( (rootcat6 == 2) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 2) )  +
                                d34[[s]] * ( (rootcat6 == 3) * (rootcat3 == 4) + (rootcat6 == 4) * (rootcat3 == 3) )  +
                                d35[[s]] * ( (rootcat6 == 3) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 3) )  +
                                d45[[s]] * ( (rootcat6 == 4) * (rootcat3 == 5) + (rootcat6 == 5) * (rootcat3 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat6 == 1) * (rootcat4 == 2) + (rootcat6 == 2) * (rootcat4 == 1) )    +
                                d13[[s]] * ( (rootcat6 == 1) * (rootcat4 == 3) + (rootcat6 == 3) * (rootcat4 == 1) )  +
                                d14[[s]] * ( (rootcat6 == 1) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 1) )  +
                                d15[[s]] * ( (rootcat6 == 1) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 1) )  +
                                d23[[s]] * ( (rootcat6 == 2) * (rootcat4 == 3) + (rootcat6 == 3) * (rootcat4 == 2) )  +
                                d24[[s]] * ( (rootcat6 == 2) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 2) )  +
                                d25[[s]] * ( (rootcat6 == 2) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 2) )  +
                                d34[[s]] * ( (rootcat6 == 3) * (rootcat4 == 4) + (rootcat6 == 4) * (rootcat4 == 3) )  +
                                d35[[s]] * ( (rootcat6 == 3) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 3) )  +
                                d45[[s]] * ( (rootcat6 == 4) * (rootcat4 == 5) + (rootcat6 == 5) * (rootcat4 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat6 == 1) * (rootcat5 == 2) + (rootcat6 == 2) * (rootcat5 == 1) )    +
                                d13[[s]] * ( (rootcat6 == 1) * (rootcat5 == 3) + (rootcat6 == 3) * (rootcat5 == 1) )  +
                                d14[[s]] * ( (rootcat6 == 1) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 1) )  +
                                d15[[s]] * ( (rootcat6 == 1) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 1) )  +
                                d23[[s]] * ( (rootcat6 == 2) * (rootcat5 == 3) + (rootcat6 == 3) * (rootcat5 == 2) )  +
                                d24[[s]] * ( (rootcat6 == 2) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 2) )  +
                                d25[[s]] * ( (rootcat6 == 2) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 2) )  +
                                d34[[s]] * ( (rootcat6 == 3) * (rootcat5 == 4) + (rootcat6 == 4) * (rootcat5 == 3) )  +
                                d35[[s]] * ( (rootcat6 == 3) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 3) )  +
                                d45[[s]] * ( (rootcat6 == 4) * (rootcat5 == 5) + (rootcat6 == 5) * (rootcat5 == 4) )  ,
                              
                              0,0,0) , 
                  
                  list(       d12[[s]] * ( (rootcat7 == 1) * (rootcat1 == 2) + (rootcat7 == 2) * (rootcat1 == 1) )    +
                                d13[[s]] * ( (rootcat7 == 1) * (rootcat1 == 3) + (rootcat7 == 3) * (rootcat1 == 1) )  +
                                d14[[s]] * ( (rootcat7 == 1) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 1) )  +
                                d15[[s]] * ( (rootcat7 == 1) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 1) )  +
                                d23[[s]] * ( (rootcat7 == 2) * (rootcat1 == 3) + (rootcat7 == 3) * (rootcat1 == 2) )  +
                                d24[[s]] * ( (rootcat7 == 2) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 2) )  +
                                d25[[s]] * ( (rootcat7 == 2) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 2) )  +
                                d34[[s]] * ( (rootcat7 == 3) * (rootcat1 == 4) + (rootcat7 == 4) * (rootcat1 == 3) )  +
                                d35[[s]] * ( (rootcat7 == 3) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 3) )  +
                                d45[[s]] * ( (rootcat7 == 4) * (rootcat1 == 5) + (rootcat7 == 5) * (rootcat1 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat7 == 1) * (rootcat2 == 2) + (rootcat7 == 2) * (rootcat2 == 1) )    +
                                d13[[s]] * ( (rootcat7 == 1) * (rootcat2 == 3) + (rootcat7 == 3) * (rootcat2 == 1) )  +
                                d14[[s]] * ( (rootcat7 == 1) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 1) )  +
                                d15[[s]] * ( (rootcat7 == 1) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 1) )  +
                                d23[[s]] * ( (rootcat7 == 2) * (rootcat2 == 3) + (rootcat7 == 3) * (rootcat2 == 2) )  +
                                d24[[s]] * ( (rootcat7 == 2) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 2) )  +
                                d25[[s]] * ( (rootcat7 == 2) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 2) )  +
                                d34[[s]] * ( (rootcat7 == 3) * (rootcat2 == 4) + (rootcat7 == 4) * (rootcat2 == 3) )  +
                                d35[[s]] * ( (rootcat7 == 3) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 3) )  +
                                d45[[s]] * ( (rootcat7 == 4) * (rootcat2 == 5) + (rootcat7 == 5) * (rootcat2 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat7 == 1) * (rootcat3 == 2) + (rootcat7 == 2) * (rootcat3 == 1) )    +
                                d13[[s]] * ( (rootcat7 == 1) * (rootcat3 == 3) + (rootcat7 == 3) * (rootcat3 == 1) )  +
                                d14[[s]] * ( (rootcat7 == 1) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 1) )  +
                                d15[[s]] * ( (rootcat7 == 1) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 1) )  +
                                d23[[s]] * ( (rootcat7 == 2) * (rootcat3 == 3) + (rootcat7 == 3) * (rootcat3 == 2) )  +
                                d24[[s]] * ( (rootcat7 == 2) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 2) )  +
                                d25[[s]] * ( (rootcat7 == 2) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 2) )  +
                                d34[[s]] * ( (rootcat7 == 3) * (rootcat3 == 4) + (rootcat7 == 4) * (rootcat3 == 3) )  +
                                d35[[s]] * ( (rootcat7 == 3) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 3) )  +
                                d45[[s]] * ( (rootcat7 == 4) * (rootcat3 == 5) + (rootcat7 == 5) * (rootcat3 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat7 == 1) * (rootcat4 == 2) + (rootcat7 == 2) * (rootcat4 == 1) )    +
                                d13[[s]] * ( (rootcat7 == 1) * (rootcat4 == 3) + (rootcat7 == 3) * (rootcat4 == 1) )  +
                                d14[[s]] * ( (rootcat7 == 1) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 1) )  +
                                d15[[s]] * ( (rootcat7 == 1) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 1) )  +
                                d23[[s]] * ( (rootcat7 == 2) * (rootcat4 == 3) + (rootcat7 == 3) * (rootcat4 == 2) )  +
                                d24[[s]] * ( (rootcat7 == 2) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 2) )  +
                                d25[[s]] * ( (rootcat7 == 2) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 2) )  +
                                d34[[s]] * ( (rootcat7 == 3) * (rootcat4 == 4) + (rootcat7 == 4) * (rootcat4 == 3) )  +
                                d35[[s]] * ( (rootcat7 == 3) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 3) )  +
                                d45[[s]] * ( (rootcat7 == 4) * (rootcat4 == 5) + (rootcat7 == 5) * (rootcat4 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat7 == 1) * (rootcat5 == 2) + (rootcat7 == 2) * (rootcat5 == 1) )    +
                                d13[[s]] * ( (rootcat7 == 1) * (rootcat5 == 3) + (rootcat7 == 3) * (rootcat5 == 1) )  +
                                d14[[s]] * ( (rootcat7 == 1) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 1) )  +
                                d15[[s]] * ( (rootcat7 == 1) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 1) )  +
                                d23[[s]] * ( (rootcat7 == 2) * (rootcat5 == 3) + (rootcat7 == 3) * (rootcat5 == 2) )  +
                                d24[[s]] * ( (rootcat7 == 2) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 2) )  +
                                d25[[s]] * ( (rootcat7 == 2) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 2) )  +
                                d34[[s]] * ( (rootcat7 == 3) * (rootcat5 == 4) + (rootcat7 == 4) * (rootcat5 == 3) )  +
                                d35[[s]] * ( (rootcat7 == 3) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 3) )  +
                                d45[[s]] * ( (rootcat7 == 4) * (rootcat5 == 5) + (rootcat7 == 5) * (rootcat5 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat7 == 1) * (rootcat6 == 2) + (rootcat7 == 2) * (rootcat6 == 1) )    +
                                d13[[s]] * ( (rootcat7 == 1) * (rootcat6 == 3) + (rootcat7 == 3) * (rootcat6 == 1) )  +
                                d14[[s]] * ( (rootcat7 == 1) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 1) )  +
                                d15[[s]] * ( (rootcat7 == 1) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 1) )  +
                                d23[[s]] * ( (rootcat7 == 2) * (rootcat6 == 3) + (rootcat7 == 3) * (rootcat6 == 2) )  +
                                d24[[s]] * ( (rootcat7 == 2) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 2) )  +
                                d25[[s]] * ( (rootcat7 == 2) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 2) )  +
                                d34[[s]] * ( (rootcat7 == 3) * (rootcat6 == 4) + (rootcat7 == 4) * (rootcat6 == 3) )  +
                                d35[[s]] * ( (rootcat7 == 3) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 3) )  +
                                d45[[s]] * ( (rootcat7 == 4) * (rootcat6 == 5) + (rootcat7 == 5) * (rootcat6 == 4) )  ,
                              
                              0,0 ) , 
                  
                  list(       d12[[s]] * ( (rootcat8 == 1) * (rootcat1 == 2) + (rootcat8 == 2) * (rootcat1 == 1) )    +
                                d13[[s]] * ( (rootcat8 == 1) * (rootcat1 == 3) + (rootcat8 == 3) * (rootcat1 == 1) )  +
                                d14[[s]] * ( (rootcat8 == 1) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 1) )  +
                                d15[[s]] * ( (rootcat8 == 1) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 1) )  +
                                d23[[s]] * ( (rootcat8 == 2) * (rootcat1 == 3) + (rootcat8 == 3) * (rootcat1 == 2) )  +
                                d24[[s]] * ( (rootcat8 == 2) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 2) )  +
                                d25[[s]] * ( (rootcat8 == 2) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 2) )  +
                                d34[[s]] * ( (rootcat8 == 3) * (rootcat1 == 4) + (rootcat8 == 4) * (rootcat1 == 3) )  +
                                d35[[s]] * ( (rootcat8 == 3) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 3) )  +
                                d45[[s]] * ( (rootcat8 == 4) * (rootcat1 == 5) + (rootcat8 == 5) * (rootcat1 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat8 == 1) * (rootcat2 == 2) + (rootcat8 == 2) * (rootcat2 == 1) )    +
                                d13[[s]] * ( (rootcat8 == 1) * (rootcat2 == 3) + (rootcat8 == 3) * (rootcat2 == 1) )  +
                                d14[[s]] * ( (rootcat8 == 1) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 1) )  +
                                d15[[s]] * ( (rootcat8 == 1) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 1) )  +
                                d23[[s]] * ( (rootcat8 == 2) * (rootcat2 == 3) + (rootcat8 == 3) * (rootcat2 == 2) )  +
                                d24[[s]] * ( (rootcat8 == 2) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 2) )  +
                                d25[[s]] * ( (rootcat8 == 2) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 2) )  +
                                d34[[s]] * ( (rootcat8 == 3) * (rootcat2 == 4) + (rootcat8 == 4) * (rootcat2 == 3) )  +
                                d35[[s]] * ( (rootcat8 == 3) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 3) )  +
                                d45[[s]] * ( (rootcat8 == 4) * (rootcat2 == 5) + (rootcat8 == 5) * (rootcat2 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat8 == 1) * (rootcat3 == 2) + (rootcat8 == 2) * (rootcat3 == 1) )    +
                                d13[[s]] * ( (rootcat8 == 1) * (rootcat3 == 3) + (rootcat8 == 3) * (rootcat3 == 1) )  +
                                d14[[s]] * ( (rootcat8 == 1) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 1) )  +
                                d15[[s]] * ( (rootcat8 == 1) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 1) )  +
                                d23[[s]] * ( (rootcat8 == 2) * (rootcat3 == 3) + (rootcat8 == 3) * (rootcat3 == 2) )  +
                                d24[[s]] * ( (rootcat8 == 2) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 2) )  +
                                d25[[s]] * ( (rootcat8 == 2) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 2) )  +
                                d34[[s]] * ( (rootcat8 == 3) * (rootcat3 == 4) + (rootcat8 == 4) * (rootcat3 == 3) )  +
                                d35[[s]] * ( (rootcat8 == 3) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 3) )  +
                                d45[[s]] * ( (rootcat8 == 4) * (rootcat3 == 5) + (rootcat8 == 5) * (rootcat3 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat8 == 1) * (rootcat4 == 2) + (rootcat8 == 2) * (rootcat4 == 1) )    +
                                d13[[s]] * ( (rootcat8 == 1) * (rootcat4 == 3) + (rootcat8 == 3) * (rootcat4 == 1) )  +
                                d14[[s]] * ( (rootcat8 == 1) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 1) )  +
                                d15[[s]] * ( (rootcat8 == 1) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 1) )  +
                                d23[[s]] * ( (rootcat8 == 2) * (rootcat4 == 3) + (rootcat8 == 3) * (rootcat4 == 2) )  +
                                d24[[s]] * ( (rootcat8 == 2) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 2) )  +
                                d25[[s]] * ( (rootcat8 == 2) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 2) )  +
                                d34[[s]] * ( (rootcat8 == 3) * (rootcat4 == 4) + (rootcat8 == 4) * (rootcat4 == 3) )  +
                                d35[[s]] * ( (rootcat8 == 3) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 3) )  +
                                d45[[s]] * ( (rootcat8 == 4) * (rootcat4 == 5) + (rootcat8 == 5) * (rootcat4 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat8 == 1) * (rootcat5 == 2) + (rootcat8 == 2) * (rootcat5 == 1) )    +
                                d13[[s]] * ( (rootcat8 == 1) * (rootcat5 == 3) + (rootcat8 == 3) * (rootcat5 == 1) )  +
                                d14[[s]] * ( (rootcat8 == 1) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 1) )  +
                                d15[[s]] * ( (rootcat8 == 1) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 1) )  +
                                d23[[s]] * ( (rootcat8 == 2) * (rootcat5 == 3) + (rootcat8 == 3) * (rootcat5 == 2) )  +
                                d24[[s]] * ( (rootcat8 == 2) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 2) )  +
                                d25[[s]] * ( (rootcat8 == 2) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 2) )  +
                                d34[[s]] * ( (rootcat8 == 3) * (rootcat5 == 4) + (rootcat8 == 4) * (rootcat5 == 3) )  +
                                d35[[s]] * ( (rootcat8 == 3) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 3) )  +
                                d45[[s]] * ( (rootcat8 == 4) * (rootcat5 == 5) + (rootcat8 == 5) * (rootcat5 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat8 == 1) * (rootcat6 == 2) + (rootcat8 == 2) * (rootcat6 == 1) )    +
                                d13[[s]] * ( (rootcat8 == 1) * (rootcat6 == 3) + (rootcat8 == 3) * (rootcat6 == 1) )  +
                                d14[[s]] * ( (rootcat8 == 1) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 1) )  +
                                d15[[s]] * ( (rootcat8 == 1) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 1) )  +
                                d23[[s]] * ( (rootcat8 == 2) * (rootcat6 == 3) + (rootcat8 == 3) * (rootcat6 == 2) )  +
                                d24[[s]] * ( (rootcat8 == 2) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 2) )  +
                                d25[[s]] * ( (rootcat8 == 2) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 2) )  +
                                d34[[s]] * ( (rootcat8 == 3) * (rootcat6 == 4) + (rootcat8 == 4) * (rootcat6 == 3) )  +
                                d35[[s]] * ( (rootcat8 == 3) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 3) )  +
                                d45[[s]] * ( (rootcat8 == 4) * (rootcat6 == 5) + (rootcat8 == 5) * (rootcat6 == 4) )  ,
                              
                              d12[[s]] * ( (rootcat8 == 1) * (rootcat7 == 2) + (rootcat8 == 2) * (rootcat7 == 1) )    +
                                d13[[s]] * ( (rootcat8 == 1) * (rootcat7 == 3) + (rootcat8 == 3) * (rootcat7 == 1) )  +
                                d14[[s]] * ( (rootcat8 == 1) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 1) )  +
                                d15[[s]] * ( (rootcat8 == 1) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 1) )  +
                                d23[[s]] * ( (rootcat8 == 2) * (rootcat7 == 3) + (rootcat8 == 3) * (rootcat7 == 2) )  +
                                d24[[s]] * ( (rootcat8 == 2) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 2) )  +
                                d25[[s]] * ( (rootcat8 == 2) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 2) )  +
                                d34[[s]] * ( (rootcat8 == 3) * (rootcat7 == 4) + (rootcat8 == 4) * (rootcat7 == 3) )  +
                                d35[[s]] * ( (rootcat8 == 3) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 3) )  +
                                d45[[s]] * ( (rootcat8 == 4) * (rootcat7 == 5) + (rootcat8 == 5) * (rootcat7 == 4) )  ,
                              0 ) 
    )
    
    
    
    emdc_settings$utilityOutside = 0
    emdc_settings$utilities = V
    emdc_settings$gamma = gamma
    emdc_settings$delta = delta
    emdc_settings$componentName = paste0("Class_",s)
    
    
    # 
    ### Compute within-class choice probabilities 
    P[[paste0("Class_",s)]] = apollo_emdc(emdc_settings, functionality)
    
    ### Take product across observation for same individual
    P[[paste0("Class_",s)]] = apollo_panelProd(P[[paste0("Class_",s)]], apollo_inputs ,functionality)
    
    s=s+1}
  
  # ### Compute latent class model probabilities
  lc_settings   = list(inClassProb = P, classProb=pi_values)
  
  
  P[["model"]] = apollo_lc(lc_settings, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(-sum(log(P)))
}
# ################################################################# #
#### MODEL ESTIMATION & OUTPUT                                   ####
# ################################################################# #



t1<-Sys.time()
t1

test = optim(apollo_beta, apollo_probabilities,  hessian = T  )


t2<-Sys.time()

timeprint = t2-t1
timeprint


conj1choice = test
stderrconj1choice=sqrt(abs(diag(solve(conj1choice$hessian))))   

writeout=cbind(conj1choice$par,stderrconj1choice)
writeout
