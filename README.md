
# GIN_validation 
Graphical independence network (GIN) analysis documentation (Pimenoff V. and Cleries R, 2021).

A GIN model was fitted to an observed data (OD) (N=126) to generate a synthetic target population (TP) of 10,000 individuals (hpvsimTP.txt), and resampling from this TP we created the simulated large database. The TP included 17 HPV types, the variables. For validation purposes, each HPV type of the TP with those obtained from generating a synthetic dataset obtained by applying SynSamGIN to a 500 subset of data sampled from the TP in the three sample size scenarios (SCs). 

That is, for SC1 (N=50) we derived an aggregated SynD of size N=50x500=25,000, for SC2 (N=100) the SynD was N=100*500=50,000, and for SC3 (N=250) it was N=250*500=125,000. For each SynD we calculated the expected frequencies of each HPV type and derived a standardized measure for comparability reasons: the expected number of individuals out of 10,000, that is, multiplying the frequencies by 100,000.

The performance of the GIN was also compared with seven commonly used oversampling algorithms: SMOTE (Synthetic Minority Over-Sampling Technique, SM) 14 and its variant Safe-Level-SMOTE (SLS) 15, Borderline-SMOTE variants 1 (BDLS1) and 2 (BDLS2), MWMOTE (Majority Weighted Minority Over-Sampling Technique, MW) 16, ADASYN (Adaptive Synthetic Sampling Approach for Imbalanced Learning) 17 and Random Oversampling (ROS).

hpvsimTP.txt -> simulated synthetic target population (TP) of 10,000 individuals.

Sim 500 Samples.RData -> Files generated from the GIN analysis comparison with seven commonly used oversampling algorithms listed.

Simulation Study Paper_v1.R -> R code for GIN analysis comparison with seven commonly used oversampling algorithms listed.
