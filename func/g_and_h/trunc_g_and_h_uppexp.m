function [v, iter] = trunc_g_and_h_uppexp(thres, g, h, loc, sca, tol)
%TRUNC_G_AND_H_UPPEXP Compute the upper expectation of the form 
%E[max{X-thres,0}] where X has a truncated g-and-h distribution
% Inputs: 
%       thres: thresholds (can be a vector)
%       g: the g parameter (g > 0)
%       h: the h parameter
%       loc: the location parameter, default is 0
%       sca: the scale parameter, default is 1
%       tol: numerical tolerance, default is 1e-8
% Outputs: 
%       v: computed upper expected values
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

[F_all, iter1] = g_and_h_cdf([0; thres], g, h, loc, sca, tol);
Fneg_0 = 1 - F_all(1);
Fneg_thres = 1 - F_all(2:end);
[Tinv, iter2] = g_and_h_inverse((thres - loc) / sca, g, h, tol);

sq_h = sqrt(1 - h);
inner = exp(0.5 * g ^ 2 / (1 - h)) * normcdf((g / (1 - h) - Tinv) ...
    * sq_h) - normcdf(-Tinv * sq_h);
v = sca / (Fneg_0 * g * sq_h) * inner ...
    + (loc - thres) .* Fneg_thres / Fneg_0;
iter = iter1 + iter2;

end

