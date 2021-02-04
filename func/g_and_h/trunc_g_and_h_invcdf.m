function [x, iter] = trunc_g_and_h_invcdf(p, g, h, loc, sca, tol)
%TRUNC_G_AND_H_INVCDF The inverse cumulative distribution function of the 
%truncated g-and-h distribution (truncated below 0)
% Inputs:
%       p: cumulative probability
%       g: the g parameter
%       h: the h parameter
%       loc: the location parameter, default is 0
%       sca: the scale parameter, default is 1
%       tol: the numerical tolerance, default is 1e-8
% Outputs:
%       x: corresponding variable value
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

if g ~= 0
    t_func = @(v)((exp(g * v) - 1) / g .* exp(0.5 * h * v .^ 2));
else
    t_func = @(v)(v .* exp(0.5 * h * v .^ 2));
end

x = zeros(length(p), 1);
x(p == 1) = inf;
list_pos = p > 0 & p < 1;

[z0, iter] = g_and_h_inverse(-loc / sca, g, h, tol);
p0 = normcdf(z0);

p_adj = p(list_pos) * (1 - p0) + p0;
x(list_pos) = t_func(norminv(p_adj)) * sca + loc;

end

