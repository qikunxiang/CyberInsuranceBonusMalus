function new_level = bonusmalus_transition(rule, cur_level, claim)
%BONUSMALUS_TRANSITION Bonus-Malus transition rule
% Inputs: 
%       rule: the Bonus-Malus rule table
%           rule(i, j) corresponds to the maximum claim for transition 
%           i -> j
%           rule(i, j - 1) corresponds to the minimum claim for such
%           transition (in the case that j = 1, the minimum is 0)
%       cur_level: the current Bonus-Malus level
%       claim: the amount of claim in the current time period
% Outputs: 
%       new_level: the updated Bonus-Malus level

new_level = sum(claim > rule(cur_level, :)) + 1;

end

