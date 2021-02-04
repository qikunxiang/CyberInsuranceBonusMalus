function [p, iter] = trunc_g_and_h_cdf(x, g, h, loc, sca, tol)
%TRUNC_G_AND_H_CDF The cumulative distribution function of the truncated
%g-and-h distribution (truncated below 0)
% Inputs:
%       x: input to the cdf
%       g: the g parameter
%       h: the h parameter
%       loc: the location parameter, default is 0
%       sca: the scale parameter, default is 1
%       tol: the numerical tolerance, default is 1e-8
% Outputs:
%       p: cumulative probability
%       iter: iteration number from the inversion

if ~exist('loc', 'var') || isempty(loc)
    loc = 0;
end

if ~exist('sca', 'var') || isempty(sca)
    sca = 1;
end

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-8;
end

p = zeros(length(x), 1);
list_pos = x > 0;

[z, iter] = g_and_h_inverse(([0; x(list_pos)] - loc) / sca, g, h, tol);
p_all = normcdf(z);
p(list_pos) = (p_all(2:end) - p_all(1)) / (1 - p_all(1));

end

