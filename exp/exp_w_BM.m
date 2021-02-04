rng(1000);

% Bonus-Malus has 4 levels.

BM = struct;

cap_base = 1000;
deductible_base = 0.5;

premium_range = [0, 7];

premium_list = (premium_range(1):0.005:premium_range(2))';
test_levels = length(premium_list);

params = struct;
params.horizon = 20;
params.discount = 0.95;
params.scale = 1;

BM.init = 3;
BM.cap = ones(4, 1) * cap_base;
BM.rule = [0, 0, 0, BM.cap(1); ...
    0, 0, 0, BM.cap(2); ...
    -1, 0, 0, BM.cap(3); ...
    -1, -1, 0, BM.cap(4)];
BM.inactive_rule = [1, 2; 1, 3; 0, 3; 1, 3];
BM.deductible = repmat(ones(4, 1) * deductible_base, 1, ...
    params.horizon);
BM.deductible(:, end) = BM.deductible(:, end) * 10;
BM.penalty_out = linspace(0, 5, params.horizon) + 3;
BM.penalty_in = [zeros(1, params.horizon - 5), linspace(0, 3, 5)];
BM.penalty_rejoin = 3;
premium_levels = [0.6; 0.8; 1.0; 1.5];

LDA = struct;
LDA.freq_mean = 0.8;
LDA.freq_var = 0.8;
LDA.g = 1.8;
LDA.h = 0.15;
LDA.loc = 0;
LDA.sca = 1;
LDA.granularity = 20;
LDA.approx_max = 10 * max(BM.cap);
LDA.tilt = 20 / (2 ^ LDA.granularity);
sev_mean = trunc_g_and_h_uppexp(0, LDA.g, LDA.h, LDA.loc, LDA.sca);

Def = struct;
Def.effect = trunc_g_and_h_invcdf([0.0; 0.7], ...
    LDA.g, LDA.h, LDA.loc, LDA.sca) * params.scale;
Def.cost = [0; 0.5];

cost_list = zeros(test_levels, 1);
detail_list = zeros(test_levels, 6);
BM_retention_list = zeros(test_levels, length(BM.cap) + 1);
Def_retention_list = zeros(test_levels, length(Def.effect));

sim_num = 1e5;
cost_samples = zeros(test_levels, sim_num);

saved_compound = [];

for test_id = 1:test_levels
    premium_base = premium_list(test_id);
    BM.premium = premium_levels * premium_base;
    
    
    if ~isempty(saved_compound)
        params.saved_compound = saved_compound;
    end
    
    [cost, policy, init, output] = optimal_provision(BM, Def, ...
        LDA, params);
    
    occupancy = output.occupancy;
    defstats = output.defstats;
    claimstats = output.claimstats;
    
    cost_list(test_id) = cost(init.BM, init.Ins);
    detail_list(test_id, 1) = output.mean.defense;
    detail_list(test_id, 2) = output.mean.premium;
    detail_list(test_id, 3) = output.mean.penalty;
    detail_list(test_id, 4) = output.mean.loss;
    detail_list(test_id, 5) = output.mean.claim;
    detail_list(test_id, 6) = output.mean.prevented;
    
    BM_retention_list(test_id, :) = sum(output.BMstats, 1);
    Def_retention_list(test_id, :) = sum(defstats, 1);
    
    % Uncomment the part below to perform Monte-Carlo simulation
    
    % em_costs = provision_sim(BM, Def, LDA, params, ...
    %     init, policy, sim_num, false, true);
    % cost_samples(test_id, :) = em_costs;
    
    if isempty(saved_compound)
        saved_compound = output.saved_compound;
    end
    
end

save('exp/data/exp_w_BM.mat', 'premium_list', 'detail_list', ...
    'BM_retention_list', 'Def_retention_list', 'cost_list', ...
    'BM', 'Def', 'LDA', 'params', 'sev_mean', 'cost_samples');