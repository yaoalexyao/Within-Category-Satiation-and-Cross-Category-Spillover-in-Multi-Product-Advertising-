# Within-Category-Satiation-and-Cross-Category-Spillover-in-Multi-Product-Advertising-
Code sharing for paper Within-Category Satiation and Cross-Category Spillover in Multi-Product Advertising 

Demand code has two parts:

1. Demand model 1 eMDC.R to estimate an extended MDC (eMDC) model
2. Demand model 2 multichoice to estimate a multichoice model

Supply code has three parts:

1. Supply 1 interest score and ad rank Table 7.R to obtain results in Table 7
2. Supply 2 vpc_real.R to infer advertiser's (daily) value-per-click based on the FOC of their profit-maximization problems.
3. Supply 3 VPC model Table 8.R to regress inferred vpc on advertiser characteristics

Counterfactual code:

For each counterfactual that contains multipl scenarios, we provide code for one scenario (using 1 day of consumer history), and the code is easy to be generalized (e.g., use 3,5,7, etc days of consumer history). 

1. Counterfactual 1_1 bideqm_1day.R to predict equilibrium bidding based on the inferred value-per-click based on 1 day of consumer history in the matching process.
2. Counterfactual 1_2 aucsim_1day.R to predict equilibrium outcomes based on based on 1 day of consumer history. 
3. Counterfactual 2 CTR and bid in allocating.R to incorporate both bids and CTR in the allocation process
4. Counterfactual 3 aucsim_minbid_1decile. R to change the reservation price (minimum bid) to be 1st deciles of the weighted bids in the ranking process. 

Due to an NDA signed with the company, we cannot share the original dataset. Following JM's alternative disclosure guidelines, we provided detailed summary statistics of important variables from Tables 1 to 5 and Figure 3 in the paper. 
