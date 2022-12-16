# A Bonus-Malus Framework for Cyber Risk Insurance and Optimal Cybersecurity Provisioning

+ By Qikun Xiang, Ariel Neufeld, Gareth W. Peters, Ido Nevat, and Anwitaman Datta
+ Article link (arXiV): https://arxiv.org/abs/2102.05568

# Description of files

+ func/main/          contains the core functions  
    - optimal\_provision.m:             dynamic programming-based algorithm for optimal cybersecurity provisioning with cyber risk insurance contract as well as self-mitigation measures 
    - provision\_sim.m:                 Monte Carlo simulation algorithm of the cybersecurity provisioning process when the policy is given 
    
+ func/bonusmalus/    contains auxiliary functions for the Bonus-Malus system
    - bonusmalus\_forget.m:             execute the rules for updating the Bonus-Malus level when the insurance contract is inactive (corresponding to the BM\_0 function in the paper)
    - bonusmalus\_inverse.m:            computes the pre-image of the Bonus-Malus transition rule, i.e. the set of claim amounts that would result in a given transition
    - bonusmalus\_transition.m:         execute the rules for updating the Bonus-Malus level when the insurance contract is active (corresponding to the BM function in the paper)

+ func/g\_and\_h/       contains functions to compute quantities related to the truncated g-and-h distribution
    - g\_and\_h\_cdf.m:                   computes the cumulative distribution function of the g-and-h distribution
    - g\_and\_h\_inverse.m:               computes the inverse g-and-h transform using bisection and Newton's method
    - trunc\_g\_and\_h\_cdf.m:             computes the cumulative distribution function of the truncated g-and-h distribution
    - trunc\_g\_and\_h\_discretize.m:      discretizes the truncated g-and-h distribution into a discrete distribution with equally-spaced atoms
    - trunc\_g\_and\_h\_invcdf.m:          computes the inverse of the cumulative distribution function of the truncated g-and-h distribution
    - trunc\_g\_and\_h\_rand.m:            randomly generates samples from the truncated g-and-h distribution
    - trunc\_g\_and\_h\_uppexp.m:          computes the upper expectation of the form E[max{X-threshold, 0}] where X has the truncated-g-and-h distribution
    - trunc\_g\_and\_h\_var.m:			   computes the variance of the truncated g-and-h distribution

+ func/lognorm/       contains functions to compute quantities related to the log-normal distribution
    - lognorm\_discretize.m:      discretizes the log-normal distribution into a discrete distribution with equally-spaced atoms
    - lognorm\_invcdf.m:          computes the inverse of the cumulative distribution function of the log-normal distribution
    - lognorm\_mm\_trunc\_g\_and\_h:		   computes the parameters of the log-normal distribution such that the first two moments match with a given truncated g-and-h distribution
    - lognorm\_partialexp.m: 		 computes the expectation of the form E[(aX+b)I_{k1<=X<=k2}] where X has the log-normal distribution
    - lognorm\_rand.m:            randomly generates samples from the log-normal distribution
    - lognorm\_uppexp.m:          computes the upper expectation of the form E[max{X-threshold, 0}] where X has the log-normal distribution
    
+ func/approx/        contains functions to compute the FFT approximation of the LDA model
    - compound\_fft\_approx.m:           approximates the probability mass function of a discrete compound distribution using the Fast Fourier Transform algorithm
    - compound\_nbin\_mitigation\_rand.m:    randomly generates samples from a compound distribution with negative binomial frequence and given severity distribution, when a self-mitigation measure is present
    - compound\_nbin\_rand.m:            randomly generates samples from a compound distribution with negative binomial frequence and given severity distribution
    - compound\_specexp\_approx.m:       computes a specific form of expectation from the FFT-approximated discrete compound distribution
    - compound\_specprob\_approx.m:      computes a specific form of probability from the FFT-approximated discrete compound distribution

+ utils/tight\_subplot: 		used for creating figures with narrow margins
+ utils/parfor\_progress:		used for showing a progress bar

+ exp/                contains the scripts to run the experiment (see below)

+ truncated\_g\_and\_h.pdf: 	contains additional notes about the truncated g-and-h distribution
+ online\_appendix.pdf: 		contains the experimental results under the modified settings and discussions

# Instruction to run the experiments

## Configurations

+ All folders and subfolders must be added to the search path. 


## Experiment in the paper

+ Run exp/exp\_original/exp\_w\_BM.m to generate the result file for the case with Bonus-Malus.
+ Run exp/exp\_original/exp\_wo\_BM.m to generate the result file for the case without Bonus-Malus.
+ Run exp/exp\_original/exp\_w\_BM\_plot.mlx to plot the figures for the case with Bonus-Malus.
+ Run exp/exp\_original/exp\_wo\_BM\_plot.mlx to plot the figures for the case without Bonus-Malus.

## Experiment under the first modified setting (changing the h parameter in the truncated g-and-h distribution)

### Case where h = 0.10

+ Run exp/exp\_diff\_h/exp10\_w\_BM.m to generate the result file for the case with Bonus-Malus.
+ Run exp/exp\_diff\_h/exp10\_wo\_BM.m to generate the result file for the case without Bonus-Malus.
+ Run exp/exp\_diff\_h/exp10\_w\_BM\_plot.mlx to plot the figures for the case with Bonus-Malus.
+ Run exp/exp\_diff\_h/exp10\_wo\_BM\_plot.mlx to plot the figures for the case without Bonus-Malus.

### Case where h = 0.20

+ Run exp/exp\_diff\_h/exp20\_w\_BM.m to generate the result file for the case with Bonus-Malus.
+ Run exp/exp\_diff\_h/exp20\_wo\_BM.m to generate the result file for the case without Bonus-Malus.
+ Run exp/exp\_diff\_h/exp20\_w\_BM\_plot.mlx to plot the figures for the case with Bonus-Malus.
+ Run exp/exp\_diff\_h/exp20\_wo\_BM\_plot.mlx to plot the figures for the case without Bonus-Malus.

### Case where h = 0.25

+ Run exp/exp\_diff\_h/exp25\_w\_BM.m to generate the result file for the case with Bonus-Malus.
+ Run exp/exp\_diff\_h/exp25\_wo\_BM.m to generate the result file for the case without Bonus-Malus.
+ Run exp/exp\_diff\_h/exp25\_w\_BM\_plot.mlx to plot the figures for the case with Bonus-Malus.
+ Run exp/exp\_diff\_h/exp25\_wo\_BM\_plot.mlx to plot the figures for the case without Bonus-Malus.

## Experiment under the second modified setting (changing the loss severity distribution to log-normal)

+ Run exp/exp\_lognorm/exp\_ln\_w\_BM.m to generate the result file for the case with Bonus-Malus.
+ Run exp/exp\_lognorm/exp\_ln\_wo\_BM.m to generate the result file for the case without Bonus-Malus.
+ Run exp/exp\_lognorm/exp\_ln\_w\_BM\_plot.mlx to plot the figures for the case with Bonus-Malus.
+ Run exp/exp\_lognorm/exp\_ln\_wo\_BM\_plot.mlx to plot the figures for the case without Bonus-Malus.
