function [mu, sig2] = lognorm_mm_trunc_g_and_h(g, h, loc, sca, tol)
%LOGNORM_MM_TRUNC_G_AND_H Fitting a log-normal distribution to a truncated
%g-and-h distribution via moment matching
% Inputs:
%       g: the g parameter
%       h: the h parameter
%       loc: the location parameter, default is 0
%       sca: the scale parameter, default is 1
%       tol: the numerical tolerance, default is 1e-8

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-8;
end

ghmean = trunc_g_and_h_uppexp(0, g, h, loc, sca, tol);
ghvar = trunc_g_and_h_var(g, h, loc, sca, tol);

logm1 = log(ghmean);
logm2 = log(ghvar + ghmean ^ 2);

mu = 2 * logm1 - 0.5 * logm2;
sig2 = logm2 - 2 * logm1;

end

