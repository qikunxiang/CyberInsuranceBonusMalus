function [new_level, new_wait] = bonusmalus_forget(rule, cur_level, wait)
%BONUSMALUS_FORGET Bonus-Malus forgetting rule
% Inputs: 
%       rule: the Bonus-Malus forgetting rule table
%           rule(i, 1) corresponds to the number of time periods to wait
%           before the transition triggers
%           rule(i, 2) corresponds to the Bonus-Malus level transitioned
%           into after waiting for the specified time periods
%       cur_level: the current Bonus-Malus level
%       wait: the number of years waited
% Outputs: 
%       new_level: the updated Bonus-Malus level
%       new_wait: the updated wait time

if wait == -1
    % if the contract was never active
    new_level = cur_level;
    new_wait = -1;
elseif rule(cur_level, 1) <= wait
    % if a transition rule is triggered, transition to a new level
    new_level = rule(cur_level, 2);
    new_wait = 1;
else
    % if a transition rule is not triggered, increment the counter
    new_level = cur_level;
    new_wait = wait + 1;
end

end

