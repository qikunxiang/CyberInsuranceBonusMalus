function v = compound_specexp_approx(atoms, p_cp, ...
    intercept_list, range_list, trunc, cap)
%COMPOUND_SPECEXP_APPROX Approximate a specitic type of expectation of the
%compound distribution, with the form E[max{X-intercept,0}I[l_min,l_max]]
% Inputs: 
%       atoms: the list of atoms in the approximated compound distribution
%       p_cp: the approximated pmf of the compound distribution
%       intercept_list: the intercept value in a vector
%       range_list: the range of integration in a matrix with two columns
%       trunc: lower truncation on X, default is 0
%       cap: upper limit of X, default is inf
% Outputs:
%       v: the computed expected values in a vector

if ~exist('trunc', 'var') || isempty(trunc)
    trunc = 0;
end

if ~exist('cap', 'var') || isempty(cap)
    cap = inf;
end

num = length(intercept_list);
assert(num == size(range_list, 1), 'lengths of input mismatch');

atoms = min(max(atoms - trunc, 0), cap);

v = zeros(num, 1);

for i = 1:num
    i_list = atoms >= max(intercept_list(i), range_list(i, 1)) ...
        & atoms <= max(intercept_list(i), range_list(i, 2));
    v(i) = sum((atoms - intercept_list(i)) .* p_cp .* i_list);
end

end

