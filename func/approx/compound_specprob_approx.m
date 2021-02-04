function p = compound_specprob_approx(atoms, p_cp, range_list, trunc, cap)
%COMPOUND_SPECPROB_APPROX Approximate a specitic type of probability of the
%compound distribution, with the form P[max{X-intercept,0}I[l_min,l_max]]
% Inputs: 
%       atoms: the list of atoms in the approximated compound distribution
%       p_cp: the approximated pmf of the compound distribution
%       range_list: the range of integration in a matrix with two columns
%       trunc: lower truncation on X, default is 0
%       cap: upper limit of X, default is inf
% Outputs:
%       v: the computed probabilities in a vector

if ~exist('trunc', 'var') || isempty(trunc)
    trunc = 0;
end

if ~exist('cap', 'var') || isempty(cap)
    cap = inf;
end

atoms = min(max(atoms - trunc, 0), cap);

num = size(range_list, 1);
p = zeros(num, 1);

for i = 1:num
    i_list = atoms >= range_list(i, 1) & atoms <= range_list(i, 2);
    p(i) = sum(p_cp .* i_list);
end

end

