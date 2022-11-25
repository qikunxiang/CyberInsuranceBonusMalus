function [v, iter] = trunc_g_and_h_var(g, h, loc, sca, tol)
%TRUNC_G_AND_H_UPPEXP Compute the variance of the truncated g-and-h
%distribution 
% Inputs: 
%       g: the g parameter (g > 0)
%       h: the h parameter
%       loc: the location parameter, default is 0
%       sca: the scale parameter, default is 1
%       tol: numerical tolerance, default is 1e-8
% Outputs: 
%       v: computed variance
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

[F_0, iter1] = g_and_h_cdf(0, g, h, loc, sca, tol);
Fneg_0 = 1 - F_0;
[Tinv, iter2] = g_and_h_inverse(-loc / sca, g, h, tol);
sq_h = sqrt(1 - h);
sq_2h = sqrt(1 - 2 * h);

v1 = sca ^ 2 / (Fneg_0 * g ^ 2 * sq_2h) * (exp(2 * g ^ 2 / (1 - 2 * h)) ...
    * normcdf((2 * g / (1 - 2 * h) - Tinv) * sq_2h) ...
    - 2 * exp(g ^ 2 / (2 * (1 - 2 * h))) ...
    * normcdf((g / (1 - 2 * h) - Tinv) * sq_2h) + normcdf(-Tinv * sq_2h));

v2 = sca ^ 2 / (Fneg_0 ^ 2 * g ^ 2 * (1 - h)) ...
    * (exp(g ^ 2 / (2 * (1 - h))) ...
    * normcdf((g / (1 - h) - Tinv) * sq_h) - normcdf(-Tinv * sq_h)) ^ 2;

v = v1 - v2;
iter = iter1 + iter2;

end

