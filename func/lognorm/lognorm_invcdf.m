function x = lognorm_invcdf(p, mu, sig2)
%LOGNORM_INVCDF The inverse cumulative distribution function of the 
%log-normal distribution
% Inputs:
%       p: cumulative probability
%       mu: the mu parameter
%       sig2: the sigma^2 parameter
% Outputs:
%       x: corresponding variable value

x = exp(norminv(p) * sqrt(sig2) + mu);

end

