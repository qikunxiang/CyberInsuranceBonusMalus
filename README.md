# A Bonus-Malus Framework for Cyber Risk Insurance and Optimal Cybersecurity Provisioning

+ By Qikun Xiang, Ariel Neufeld, Gareth W. Peters, Ido Nevat, and Anwitaman Datta
+ Article link (arXiV): https://arxiv.org/abs/2102.05568

# Description of files

+ func/main/          contains the core functions  
    - optimal_provision.m:             dynamic programming-based algorithm for optimal cybersecurity provisioning with cyber risk insurance contract as well as self-mitigation measures 
    - provision_sim.m:                 Monte Carlo simulation algorithm of the cybersecurity provisioning process when the policy is given 
    
+ func/bonusmalus/    contains auxiliary functions for the Bonus-Malus system
    - bonusmalus_forget.m:             execute the rules for updating the Bonus-Malus level when the insurance contract is inactive (corresponding to the BM_0 function in the paper)
    - bonusmalus_inverse.m:            computes the pre-image of the Bonus-Malus transition rule, i.e. the set of claim amounts that would result in a given transition
    - bonusmalus_transition.m:         execute the rules for updating the Bonus-Malus level when the insurance contract is active (corresponding to the BM function in the paper)

+ func/g_and_h/       contains functions to compute quantities related to the truncated g-and-h distribution
    - g_and_h_cdf.m:                   computes the cumulative distribution function of the g-and-h distribution
    - g_and_h_inverse.m:               computes the inverse g-and-h transform using bisection and Newton's method
    - trunc_g_and_h_cdf.m:             computes the cumulative distribution function of the truncated g-and-h distribution
    - trunc_g_and_h_discretize.m:      discretizes the truncated g-and-h distribution into a discrete distribution with equally-spaced atoms
    - trunc_g_and_h_invcdf.m:          computes the inverse of the cumulative distribution function of the truncated g-and-h distribution
    - trunc_g_and_h_rand.m:            randomly generates samples from the truncated g-and-h distribution
    - trunc_g_and_h_uppexp.m:          computes the upper expectation of the form E(max{X-threshold, 0})
    
+ func/approx/        contains functions to compute the FFT approximation of the LDA model
    - compound_fft_approx.m:           approximates the probability mass function of a discrete compound distribution using the Fast Fourier Transform algorithm
    - compound_nbin_defense_rand.m:    randomly generates samples from a compound distribution with negative binomial frequence and given severity distribution, when a self-mitigation measure is present
    - compound_nbin_rand.m:            randomly generates samples from a compound distribution with negative binomial frequence and given severity distribution
    - compound_specexp_approx.m:       computes a specific form of expectation from the FFT-approximated discrete compound distribution
    - compound_specprob_approx.m:      computes a specific form of probability from the FFT-approximated discrete compound distribution

+ exp/                contains the scripts to run the experiment (see later)

+ tight_subplot/:     used for creating figures with narrow margins

# Instruction to run the experiments

## Configurations

+ All folders and subfolders must be added to the search path. 


## Experiment

+ Run exp/exp_w_BM.m to generate the result file for the case with Bonus-Malus.
+ Run exp/exp_wo_BM.m to generate the result file for the case without Bonus-Malus.
+ Run exp/exp_w_BM_plot.mlx to plot the figures for the case with Bonus-Malus.
+ Run exp/exp_wo_BM_plot.mlx to plot the figures for the case without Bonus-Malus.
