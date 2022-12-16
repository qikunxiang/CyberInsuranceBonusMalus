function v = lognorm_uppexp(thres, mu, sig2)
%LOGNORM_UPPEXP Compute the upper expectation of the form 
%E[max{X-thres,0}] where X has a log-normal distribution
% Inputs: 
%       thres: thresholds (can be a vector)
%       mu: the mu parameter
%       sig2: the sigma^2 parameter
% Outputs: 
%       v: computed upper expected values

v = lognorm_partialexp(mu, sig2, thres, inf, 1, -thres);

end

