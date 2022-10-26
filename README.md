# RCode_transportability
The folders and the R code help to reproducing Tables 2-4 and Figures 1-2  in the article "Transporting randomized trial results to estimate counterfactual survival functions in target populations" by Zhiqiang Cao, Youngjoo Cho and Fan Li (under review) 

For questions or comments about the code, please contact Zhiqiang Cao <zcaoae@connect.ust.hk> or Fan Li at <fan.f.li@yale.edu>. 
You will need to change the directory to use the example code in script example_analysis.R This folder includes the following functions:

 I. gen_true_value_for_delta.R is to determine ture value of delta(t) for three specified survival times in simulations 

2. sim_code_taste_finla.R are to compute DR1 and DR2 estimators as well as influence functions. Typically, survres_z is to compute DR1, survres_z2 is to compute DR2, survres_z2app is DR2 used in real data analysis 

3. simulation_covariate_dependent_censor.R and simulation_random_censor.R are to compare weighted KM, two IPWs and two DRs for transporting under weak and strong sampling cases
 
4. generate_figure1.R is to produce a figure similar to Figure 1 of manuscript  5. real_data_analysis.R is for real data analysis
