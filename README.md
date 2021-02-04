# A Bonus-Malus Framework for Cyber Risk Insurance and Optimal Cybersecurity Provisioning

+ By Qikun Xiang, Ariel Neufeld, Gareth W. Peters, Ido Nevat, and Anwitaman Datta

# Description of files

+ func/main/          contains the core functions  
    - optimal_provision.m:             function 
    - provision_sim.m:                 function 
    
+ func/bonusmalus/    contains auxiliary functions for the Bonus-Malus system
    - bonusmalus_forget.m:             X
    - bonusmalus_inverse.m:            X
    - bonusmalus_transition.m:         X

+ func/g_and_h/       contains functions to compute quantities related to the truncated g-and-h distribution
    - g_and_h_cdf.m:                   X
    - g_and_h_inverse.m:               X
    - trunc_g_and_h_cdf.m:             X
    - trunc_g_and_h_discretize.m:      X
    - trunc_g_and_h_invcdf.m:          X
    - trunc_g_and_h_rand.m:            X
    - trunc_g_and_h_uppexp.m:          X
    
+ func/approx/        contains functions to compute the FFT approximation of the LDA model
    - compound_fft_approx.m:           X
    - compound_nbin_defense_rand.m:    X
    - compound_nbin_rand.m:            X
    - compound_specexp_approx.m:       X
    - compound_specprob_approx.m:      X

+ exp/                contains the scripts to run the experiment (see later)

# Instruction to run the experiments

## Configurations

+ All folders and subfolders must be added to the search path. 


## Experiment

+ Run exp/exp_w_BM.m to generate the result file for the case with Bonus-Malus.
+ Run exp/exp_wo_BM.m to generate the result file for the case without Bonus-Malus.
+ Run exp/exp_w_BM_plot.mlx to plot the figures for the case with Bonus-Malus.
+ Run exp/exp_wo_BM_plot.mlx to plot the figures for the case without Bonus-Malus.
