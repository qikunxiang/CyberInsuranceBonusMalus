rng(1000);

% No Bonus-Malus.

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

BM.init = 1;
BM.cap = ones(1, 1) * cap_base;
BM.rule = BM.cap(1);
BM.inactive_rule = [0, 1];
BM.deductible = repmat(ones(1, 1) * deductible_base, 1, ...
    params.horizon);
BM.deductible(:, end) = BM.deductible(:, end) * 10;
BM.penalty_out = linspace(0, 5, params.horizon) + 3;
BM.penalty_in = [zeros(1, params.horizon - 5), linspace(0, 3, 5)];
BM.penalty_rejoin = 3;
premium_levels = 1;

LDA = struct;
LDA.freq_mean = 0.8;
LDA.freq_var = 0.8;
LDA.g = 1.8;
LDA.h = 0.20;
LDA.loc = 0;
LDA.sca = 1;
LDA.granularity = 20;
LDA.approx_max = 10 * max(BM.cap);
LDA.tilt = 20 / (2 ^ LDA.granularity);
sev_mean = trunc_g_and_h_uppexp(0, LDA.g, LDA.h, LDA.loc, LDA.sca);

Miti = struct;
Miti.effect = trunc_g_and_h_invcdf([0.0; 0.7], ...
    LDA.g, LDA.h, LDA.loc, LDA.sca) * params.scale;
Miti.cost = [0; 0.5];

cost_list = zeros(test_levels, 1);
detail_list = zeros(test_levels, 6);
BM_retention_list = zeros(test_levels, length(BM.cap) + 1);
Miti_retention_list = zeros(test_levels, length(Miti.effect));

sim_num = 1e5;
cost_samples = zeros(test_levels, sim_num);

saved_compound = [];

parfor_progress(test_levels);
for test_id = 1:test_levels
    premium_base = premium_list(test_id);
    BM.premium = premium_levels * premium_base;
    
    
    if ~isempty(saved_compound)
        params.saved_compound = saved_compound;
    end
    
    [cost, policy, init, output] = optimal_provision(BM, Miti, ...
        LDA, params);
    
    occupancy = output.occupancy;
    mitistats = output.mitistats;
    claimstats = output.claimstats;
    
    cost_list(test_id) = cost(init.BM, init.Ins);
    detail_list(test_id, 1) = output.mean.mitigation;
    detail_list(test_id, 2) = output.mean.premium;
    detail_list(test_id, 3) = output.mean.penalty;
    detail_list(test_id, 4) = output.mean.loss;
    detail_list(test_id, 5) = output.mean.claim;
    detail_list(test_id, 6) = output.mean.prevented;
    
    BM_retention_list(test_id, :) = sum(output.BMstats, 1);
    Miti_retention_list(test_id, :) = sum(mitistats, 1);
    
    % Uncomment the part below to perform Monte-Carlo simulation
    
    % em_costs = provision_sim(BM, Miti, LDA, params, ...
    %     init, policy, sim_num, false, true);
    % cost_samples(test_id, :) = em_costs;
    
    if isempty(saved_compound)
        saved_compound = output.saved_compound;
    end
    
    parfor_progress;
end
parfor_progress(0);

save('exp/exp_diff_h/exp20_wo_BM.mat', 'premium_list', 'detail_list', ...
    'BM_retention_list', 'Miti_retention_list', 'cost_list', ...
    'BM', 'Miti', 'LDA', 'params', 'sev_mean', 'cost_samples');