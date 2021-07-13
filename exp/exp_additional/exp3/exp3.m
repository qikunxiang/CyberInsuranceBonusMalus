rng(1000);

% This script is used to quickly identify an 'optimal' base premium given a
% fixed cyber insurance contract and cyber loss model by using a bisection
% heuristic.
% The procedure starts with the interval [0, 10] and uses a bisection
% strategy to find 'the highest base premium such that on average a
% rational insured has a non-zero probability to adopt the cyber insurance
% contract'

% Bonus-Malus has 6 levels.

BM = struct;

cap_base = 1000;
deductible_base = 2;

params = struct;
params.horizon = 20;
params.discount = 0.95;
params.scale = 1;

% using Bonus-Malus with 6 levels 

BM.init = 4;
BM.cap = ones(6, 1) * cap_base;
BM.rule = [0, 0, 0, 0, 100, BM.cap(1); ...
    0, 0, 0, 0, 100, BM.cap(2); ...
    -1, 0, 0, 0, 100, BM.cap(3); ...
    -1, -1, 0, 0, 100, BM.cap(4); ...
    -1, -1, -1, 0, 100, BM.cap(5); ...
    -1, -1, -1, -1, 0, BM.cap(6)];
BM.inactive_rule = [1, 2; 1, 3; 1, 4; 0, 4; 1, 4; 1, 5];
BM.deductible = repmat([1; 1; 1; 1; 10; 10] * deductible_base, 1, ...
    params.horizon);
BM.deductible(1:4, end) = 5;
BM.penalty_out = linspace(0, 5, params.horizon) + inf;
BM.penalty_in = [zeros(1, params.horizon - 1), inf];
BM.penalty_rejoin = 3;
premium_levels = [0.6; 0.73; 0.87; 1.0; 2; 3];

% not using Bonus-Malus

% BM.init = 1;
% BM.cap = ones(1, 1) * cap_base;
% BM.rule = BM.cap(1);
% BM.inactive_rule = [0, 1];
% BM.deductible = repmat(ones(1, 1) * deductible_base, 1, ...
%     params.horizon);
% BM.penalty_out = linspace(0, 5, params.horizon) + 3;
% BM.penalty_in = [zeros(1, params.horizon - 5), linspace(0, 3, 5)];
% BM.penalty_rejoin = 3;
% premium_levels = 1;

LDA = struct;
LDA.freq_mean = 0.8;
LDA.freq_var = 0.8;
LDA.g = 1.8;
LDA.h = 0.15;
LDA.loc = 0;
LDA.sca = 1;
LDA.granularity = 16;
LDA.approx_max = 10 * max(BM.cap);
LDA.tilt = 20 / (2 ^ LDA.granularity);
sev_mean = trunc_g_and_h_uppexp(0, LDA.g, LDA.h, LDA.loc, LDA.sca);

% risk profiles 1 to 8
risk_profile = 8;
switch risk_profile
    case 1
        % freq mean +, freq var -, sev mean -, sev tail -
        LDA.freq_mean = 0.8;
        LDA.freq_var = 0.8;
        LDA.h = 0.15;
        LDA.sca = 1;
    case 2
        % freq mean +, freq var +, sev mean -, sev tail -
        LDA.freq_mean = 0.8;
        LDA.freq_var = 4;
        LDA.h = 0.15;
        LDA.sca = 1;
    case 3
        % freq mean +, freq var -, sev mean -, sev tail +
        LDA.freq_mean = 0.8;
        LDA.freq_var = 0.8;
        LDA.h = 0.3;
        LDA.sca = 0.5808558239;
    case 4
        % freq mean +, freq var +, sev mean -, sev tail +
        LDA.freq_mean = 0.8;
        LDA.freq_var = 4;
        LDA.h = 0.3;
        LDA.sca = 0.5808558239;
    case 5
        % freq mean -, freq var -, sev mean +, sev tail -
        LDA.freq_mean = 0.08;
        LDA.freq_var = 0.08;
        LDA.h = 0.15;
        LDA.sca = 10;
    case 6
        % freq mean -, freq var +, sev mean +, sev tail -
        LDA.freq_mean = 0.08;
        LDA.freq_var = 0.4;
        LDA.h = 0.15;
        LDA.sca = 10;
    case 7
        % freq mean -, freq var -, sev mean +, sev tail +
        LDA.freq_mean = 0.08;
        LDA.freq_var = 0.08;
        LDA.h = 0.3;
        LDA.sca = 0.5808558239 * 10;
    case 8
        % freq mean -, freq var +, sev mean +, sev tail +
        LDA.freq_mean = 0.08;
        LDA.freq_var = 0.4;
        LDA.h = 0.3;
        LDA.sca = 0.5808558239 * 10;
end

Miti = struct;
Miti.effect = trunc_g_and_h_invcdf([0.0; 0.7], ...
    LDA.g, LDA.h, LDA.loc, LDA.sca) * params.scale;
Miti.cost = [0; 0.5];

saved_compound = [];


premium_range = [0, 10];
current_range = premium_range;
tol = 1e-4;

while current_range(2) - current_range(1) > tol
    premium_base = 0.5 * sum(current_range);
    BM.premium = premium_levels * premium_base;
    
    
    if ~isempty(saved_compound)
        params.saved_compound = saved_compound;
    end
    
    [cost_c, policy_c, init_c, output_c] = optimal_provision(BM, Miti, ...
        LDA, params);
    
    if output_c.mean.premium == 0
        current_range(2) = premium_base;
    else
        current_range(1) = premium_base;
        cost = cost_c;
        policy = policy_c;
        init = init_c;
        output = output_c;
    end
    
end

profit_gap = output.mean.premium + output.mean.penalty ...
    - output.mean.claim;
miti_adoption = mean(output.mitistats(:, 2));

display(premium_base)
display(profit_gap)
display(miti_adoption)