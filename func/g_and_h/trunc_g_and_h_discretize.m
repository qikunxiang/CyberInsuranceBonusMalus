function [atoms, P, gap, iter] = trunc_g_and_h_discretize(range, num, ...
    g, h, loc, sca, trunc, tol)
%TRUNC_G_AND_H_DISCRETIZE Discretize the truncated g-and-h distribution
%into a distribution with equally spaced atoms
% Inputs: 
%       range: lower and upper ends of the discretization range
%       num: the number of atoms
%       g: the g parameter
%       h: the h parameter
%       loc: the location parameter, default is 0
%       sca: the scale parameter, default is 1
%       trunc: a vector of truncation points, default is [0]
%       tol: numerical tolerance, default is 1e-8
% Outputs: 
%       atoms: the location of the atoms
%       P: probability mass for each atom (matrix)
%       gap: the gap between neighboring atoms
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

if ~exist('trunc', 'var') || isempty(trunc)
    trunc = 0;
end

gap = (range(2) - range(1)) / (num - 1);

atoms = range(1) + (0:num - 1)' * gap;
edges = [range(1) - gap / 2; atoms + gap / 2];

P = zeros(num, length(trunc));
iter = 0;
for i = 1:length(trunc)
    F = zeros(length(trunc), 1);
    list_pos = edges >= 0;
    [F(list_pos), iter1] = trunc_g_and_h_cdf( ...
        edges(list_pos) + trunc(i), g, h, loc, sca, tol);
    P(:, i) = F(2:end) - F(1:end - 1);
    iter = iter + iter1;
end

end

