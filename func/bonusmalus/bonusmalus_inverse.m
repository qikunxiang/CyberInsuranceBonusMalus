function range = bonusmalus_inverse(BM, cur_level, new_levels)
%BONUSMALUS_INVERSE The inverse of the Bonus-Malus transition rule, that
%is, the claim that would result in the transition cur_level -> new_level
% Inputs: 
%       BM: the Bonus-Malus rule table
%           BM(i, j) corresponds to the maximum claim for transition i -> j
%           BM(i, j - 1) corresponds to the minimum claim for such
%           transition (in the case that j = 1, the minimum is 0)
%       cur_level: the current Bonus-Malus level
%       new_levels: the updated Bonus-Malus levels (vector)
% Outputs: 
%       range: two column matrix representing the range of claim

range = [zeros(length(new_levels), 1), BM(cur_level, new_levels)'];
range(new_levels > 1, 1) = max(BM(cur_level, ...
    new_levels(new_levels > 1) - 1)', 0) * (1 + eps) + eps;

% mark all empty intervals as infinite
range(range(:, 2) < range(:, 1), :) = inf;

end

