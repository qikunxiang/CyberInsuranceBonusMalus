function [atoms, P, gap] = lognorm_discretize(range, num, mu, sig2, trunc)
%LOGNORM_DISCRETIZE Discretize the log-normal distribution
%into a distribution with equally spaced atoms
% Inputs: 
%       range: lower and upper ends of the discretization range
%       num: the number of atoms
%       mu: the mu parameter
%       sig2: the sigma^2 parameter
%       trunc: a vector of truncation points, default is [0]
% Outputs: 
%       atoms: the location of the atoms
%       P: probability mass for each atom (matrix)
%       gap: the gap between neighboring atoms

if ~exist('trunc', 'var') || isempty(trunc)
    trunc = 0;
end

gap = (range(2) - range(1)) / (num - 1);

atoms = range(1) + (0:num - 1)' * gap;
edges = [range(1) - gap / 2; atoms + gap / 2];

P = zeros(num, length(trunc));
for i = 1:length(trunc)
    F = zeros(length(trunc), 1);
    list_pos = edges >= 0;
    F(list_pos) = normcdf((log(edges(list_pos) + trunc(i)) ...
        - mu) ./ sqrt(sig2));
    P(:, i) = F(2:end) - F(1:end - 1);
end

end

