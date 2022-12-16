function [samp, prevented] = compound_nbin_mitigation_rand(n, ...
    rate_mean, rate_var, sev_gen, mitigation)
%COMPOUND_POISSON_MITIGATION_RAND Randomly generate samples from the
%compound negative binomial distribution with a severity distribution with
%mitigation
% Inputs: 
%       n: the number of samples
%       rate_mean: the mean of the negative binomial distribution 
%       rate_var: the variance of the negative binomial distribution
%       sev_gen: a function that randomly generates samples from the
%       	severity distribution, takes a single input which is the number
%           of samples, returns the untruncated loss
%       mitigation: the loss reduction effect of mitigation
% Outputs: 
%       samp: samples in a vector
%       prevented: prevented losses in a vector

samp = zeros(n, 1);
prevented = zeros(n, 1);

if abs(rate_mean - rate_var) < eps
    N = poissrnd(rate_mean, n, 1);
else
    nbin_p = rate_mean / rate_var;
    nbin_r = rate_mean * nbin_p / (1 - nbin_p);
    N = nbinrnd(nbin_r, nbin_p, n, 1);
end

S_raw = sev_gen(sum(N));
S = max(S_raw - mitigation, 0);
S_p = min(S_raw, mitigation);
list_pos = N > 0;

if sum(list_pos) > 0
    samp_pos = accumarray(repelem((1:sum(list_pos))', N(list_pos), 1), S);
    samp(list_pos) = samp_pos;
    pr_pos = accumarray(repelem((1:sum(list_pos))', N(list_pos), 1), S_p);
    prevented(list_pos) = pr_pos;
end

end

