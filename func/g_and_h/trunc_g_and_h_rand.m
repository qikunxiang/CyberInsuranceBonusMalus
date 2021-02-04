function [samp, iter] = trunc_g_and_h_rand(n, g, h, loc, sca, tol)
%TRUNC_G_AND_H_RAND Randomly generate samples from the truncated g-and-h
%distribution
% Inputs:
%       n: number of samples
%       g: the g parameter
%       h: the h parameter
%       loc: the location parameter, default is 0
%       sca: the scale parameter, default is 1
%       tol: the numerical tolerance, default is 1e-8
% Outputs:
%       samp: samples in a vector
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

[samp, iter] = trunc_g_and_h_invcdf(rand(n, 1), g, h, loc, sca, tol);

end

