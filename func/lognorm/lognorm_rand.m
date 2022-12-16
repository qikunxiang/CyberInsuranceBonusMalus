function samp = lognorm_rand(n, mu, sig2)
%LOGNORM_RAND Randomly generate samples from the log-normal
%distribution
% Inputs:
%       n: number of samples
%       mu: the mu parameter
%       sig2: the sigma^2 parameter
% Outputs:
%       samp: samples in a vector

samp = exp(rnorm(n, 1) * sqrt(sig2) + mu);

end

