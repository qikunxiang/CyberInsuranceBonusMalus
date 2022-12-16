function [cost, policy, init, output] = optimal_provision(BM, Miti, ...
    LDA, params)
%OPTIMAL_PROVISION Optimal cyber security provisioning with dynamic
%programming
% Inputs: 
%       BM: cyber insurance with Bonus-Malus (struct)
%           BM.rule: matrix representing transition between BM levels
%           BM.inactive_rule: transition between BM levels when inactive
%           BM.premium: premium for each BM level
%           BM.cap: cap on the amount of claim for each BM level and each
%               time step
%           BM.deductible: deductible for each BM level and each time step
%           BM.penalty_out: penalty for withdrawing the contract
%           BM.penalty_in: penalty for signing the contract late
%           BM.penalty_rejoin: penalty for rejoining the contract
%       Miti: cyber mitigation (struct)
%           Miti.cost: cost of each mitigation level
%           Miti.effect: effect of each mitigation level
%       LDA: LDA model which is a compound Poisson with g-and-h severity
%           distribution
%           LDA.freq_mean: mean of the frequency distribution
%           LDA.freq_var: variance of the frequency distribution
%           LDA.g: the g parameter of the g-and-h distribution
%           LDA.h: the h parameter of the g-and-h distribution
%           LDA.loc: the location parameter of the g-and-h distribution
%           LDA.sca: the scale parameter of the g-and-h distribution
%           LDA.tol: the numerical tolerance for computation with g-and-h
%           	distribution, default is 1e-8
%           LDA.lognorm_sev: boolean indicating whether the severity
%               distribution is changed to log-normal (default is false)
%           LDA.mu: the mu parameter of the log-normal distribution
%           LDA.sig2: the sigma^2 parameter of the log-normal distribution
%           LDA.granularity: specifies the number of atoms (2^granularity)
%           	in the approximation of the compound distribution
%           LDA.approx_max: specifies the range beyond BM.cap to
%               perform the discretization and FFT, default is max(BM.cap)
%           LDA.tilt: the tilting parameter used in the FFT approximation,
%               default is (20 / (2 ^ LDA.granularity))
%       params: other parameters
%           params.horizon: the total number of time periods
%           params.discount: the discount factor, default is 1
%           params.scale: the scale of the company, default is 1
%           params.saved_compound: saved discrete approximation of the
%               compound loss distribution that is reused
% Outputs: 
%       cost: the optimal expected cost for each initial state
%       policy: the optimal policy (cell), each cell contains a struct
%       init: preferred initial state
%       output: struct containing additional outputs:
%           output.costs: the cost matrix at each time step
%           output.occupancy: the occupancy probability of states, which is
%               a cell array of T+1 tables
%           output.transition: the transition probability under the optimal
%               control, which is a cell array containing cell arrarys that
%               contain struct
%           output.BMstats: the table of statistics about Bonus-Malus
%               levels
%           output.mitistats: the table of statistics about mitigation
%           output.claimstats: the table of statistics about claim
%           output.mean: struct containing the expected value of some
%               quantities
%           output.saved_compound: struct containing the discrete
%               approximation of the compound loss distribution with
%               various mitigation, can be reused
%           output.risk_premium: the risk premium ignoring Bonus-Malus for
%               each mitigation level

if size(BM.deductible, 2) == 1
    BM.deductible = repmat(BM.deductible, 1, params.horizon);
end

if size(BM.cap, 2) == 1
    BM.cap = repmat(BM.cap, 1, params.horizon);
end

if ~isfield(LDA, 'tol')
    LDA.tol = 1e-8;
end

if ~isfield(LDA, 'lognorm_sev')
    LDA.lognorm_sev = false;
end

if ~isfield(LDA, 'approx_max')
    LDA.approx_max = max(max(BM.cap));
end

if ~isfield(LDA, 'tilt')
    LDA.tilt = 20 / (2 ^ LDA.granularity);
end

if ~isfield(params, 'discount')
    params.discount = 1;
end

if ~isfield(params, 'scale')
    params.scale = 1;
end

T = params.horizon;
discount = params.discount;

BM_num = length(BM.premium);
assert(BM_num == size(BM.rule, 1) && BM_num == size(BM.rule, 2), ...
    'mis-specified Bonus-Malus');
assert(BM_num == size(BM.inactive_rule, 1) ...
    && 2 == size(BM.inactive_rule, 2), ...
    'mis-specified Bonue-Malus inactive rules');
BM_wait_max = max(max(BM.inactive_rule(:, 1)), 1);
Miti_num = length(Miti.cost);
assert(Miti_num == length(Miti.effect), 'mis-specified mitigation');

if ~isfield(params, 'saved_compound')
    % compute the upper expectations
    if ~LDA.lognorm_sev
        upper_exp = trunc_g_and_h_uppexp(Miti.effect, LDA.g, LDA.h, ...
            LDA.loc, LDA.sca, LDA.tol) * LDA.freq_mean;
    else
        upper_exp = lognorm_uppexp(Miti.effect, LDA.mu, LDA.sig2) ...
            * LDA.freq_mean;
    end
    
    % compute the approximate compound distribution
    atom_num = 2 ^ LDA.granularity;
    if ~LDA.lognorm_sev
        [atoms, P_trunc] = trunc_g_and_h_discretize( ...
            [0, LDA.approx_max], atom_num, LDA.g, LDA.h, ...
            LDA.loc, LDA.sca, Miti.effect, LDA.tol);
    else
        [atoms, P_trunc] = lognorm_discretize([0, LDA.approx_max], ...
            atom_num, LDA.mu, LDA.sig2, Miti.effect);
    end
    
    if abs(LDA.freq_mean - LDA.freq_var) < eps
        freq_pgf = @(s)(exp(LDA.freq_mean * (s - 1)));
    else
        nbin_p = 1 - LDA.freq_mean / LDA.freq_var;
        nbin_r = LDA.freq_mean * (1 - nbin_p) / nbin_p;
        freq_pgf = @(s)(((1 - nbin_p) ./ (1 - nbin_p * s)) .^ nbin_r);
    end
    
    P_cp = compound_fft_approx(P_trunc, freq_pgf, LDA.tilt);
    P_cp(end, :) = 1 - sum(P_cp(1:end - 1, :), 1);
    
    saved_compound = struct;
    saved_compound.atoms = atoms;
    saved_compound.P = P_cp;
    saved_compound.upper_exp = upper_exp;
else
    saved_compound = params.saved_compound;
    atoms = saved_compound.atoms;
    P_cp = saved_compound.P;
    upper_exp = saved_compound.upper_exp;
end

upper_exp = upper_exp * params.scale;
atoms = atoms * params.scale;

% compute the risk premium
risk_premium = zeros(Miti_num, 1);
for Miti_id = 1:Miti_num
    risk_premium(Miti_id) = compound_specexp_approx(atoms, ...
        P_cp(:, Miti_id), 0, [0, inf], BM.deductible(BM.init, 1), ...
        BM.cap(BM.init, 1));
end

% start dynamic programming
policy = cell(T, 2);
costs = cell(T + 1, 1);
transition = cell(T, 1);

% initialize the objective function
costs{1} = zeros(BM_num, BM_wait_max + 2);

for t = 1:T
    costs{t + 1} = zeros(BM_num, BM_wait_max + 2);
    policy{T - t + 1, 2} = cell(BM_num, 1);
    policy{T - t + 1, 1} = cell(BM_num, BM_wait_max + 2);
    transition{T - t + 1} = cell(BM_num, BM_wait_max + 2);
    
    % for each Bonus-Malus level
    for BM_id = 1:BM_num
        % when insurance is active (otherwise there is no control)
        BM_rule_t = BM.rule;
        BM_rule_t = min(BM_rule_t, BM.cap(:, t));
        BM_lowest = bonusmalus_transition(BM_rule_t, BM_id, 0);
        BM_highest = bonusmalus_transition(BM_rule_t, BM_id, ...
            BM.cap(BM_id, T - t + 1));
        BM_wait = max(BM.inactive_rule(BM_id, 1), 1);
        BM_inverse_ranges = bonusmalus_inverse(BM_rule_t, BM_id, ...
            BM_lowest:BM_highest);
        claim_intercept = discount * (costs{t}(BM_lowest:BM_highest, 2) ...
            - costs{t}(BM_lowest, 2));
        
        % list of possible Bonus-Malus transitions
        BM_possible_trans = (BM_lowest:BM_highest)';
        
        % list of claiming ranges (some may be empty)
        claim_ranges = [max(BM_inverse_ranges(:, 1), claim_intercept ...
            * (1 + eps) + eps), BM_inverse_ranges(:, 2)];
        
        % remove the empty intervals
        claim_empty_list = claim_ranges(:, 2) < claim_ranges(:, 1) ...
            | all(isinf(claim_ranges), 2);
        claim_ranges = claim_ranges(~claim_empty_list, :);
        % BM level transitioning into after making the claim
        BM_claim_trans = BM_possible_trans(~claim_empty_list);
        if claim_ranges(end, 2) == BM.cap(BM_id, T - t + 1)
            % if the optimal control is to claim when the loss reaches the
            % cap, then the control should be to claim with even larger
            % losses
            claim_ranges(end, 2) = inf;
        end
        
        % shift the intervals due to the deductible
        claim_ranges = claim_ranges + BM.deductible(BM_id, T - t + 1);
        
        policy{T - t + 1, 2}{BM_id} = struct('range', ...
            claim_ranges, 'BM', BM_claim_trans);
        
        % compute expected cost for solving the provision stage
        % the control to activate insurance
        expected_cost1 = zeros(Miti_num, 1);
        
        for Miti_id = 1:Miti_num
            expected_cost1(Miti_id) = discount * costs{t}(BM_lowest, 2) ...
                - sum(compound_specexp_approx(atoms, ...
                P_cp(:, Miti_id), claim_intercept, BM_inverse_ranges, ...
                BM.deductible(BM_id, T - t + 1), ...
                BM.cap(BM_id, T - t + 1)));
        end
        
        % compute the combined cost for each control
        % iterate through all possible insurance states (1 denotes that
        % insurance was never active, 2 denotes that insurance is active)
        for Ins_id = 1:BM_wait + 2
            wait_id = Ins_id - 2;
            expected_future_cost = zeros(Miti_num, 2);
            expected_future_cost(:, 1) = expected_cost1;
            % the control to withdraw insurance
            [BM_up, wait_up] = bonusmalus_forget(BM.inactive_rule, ...
                BM_id, wait_id);
            expected_future_cost(:, 2) = discount ...
                * costs{t}(BM_up, wait_up + 2);
            
            combined_cost = expected_future_cost + Miti.cost + upper_exp;
            combined_cost(:, 1) = combined_cost(:, 1) + BM.premium(BM_id);
            if Ins_id == 2
                % if insurance is active
                combined_cost(:, 2) = combined_cost(:, 2) ...
                    + BM.penalty_out(T - t + 1);
            else
                % if insurance is inactive
                combined_cost(:, 1) = combined_cost(:, 1) ...
                    + BM.penalty_in(T - t + 1);
                
                if Ins_id > 2
                    % if insurance has been withdrawn before provision 
                    % stage
                    combined_cost(:, 1) = combined_cost(:, 1) ...
                        + BM.penalty_rejoin;
                end
            end
            
            % choose the optimal control
            [min_cost, best_Miti] = min(combined_cost(:));
            if best_Miti > Miti_num
                best_Ins = 2;
                best_Miti = best_Miti - Miti_num;
            else
                best_Ins = 1;
            end
            
            costs{t + 1}(BM_id, Ins_id) = min_cost;
            policy{T - t + 1, 1}{BM_id, Ins_id} = ...
                struct('Miti', best_Miti, 'Ins', best_Ins);
            
            if best_Ins == 1
                % if the optimal control is to activate insurance, compute
                % the transition probabilities under the optimal control
                BM_trans_prob = compound_specprob_approx(atoms, ...
                    P_cp(:, best_Miti), claim_ranges);
                no_loss_prob = P_cp(1, best_Miti);
                no_claim_prob = 1 - sum(BM_trans_prob);
                
                BM_lowest_index = find(BM_claim_trans == BM_lowest, 1);
                if ~isempty(BM_lowest_index)
                    % if the lowest possible Bonus-Malus transition is
                    % already considered, add the remaining probability to
                    % this case
                    BM_all_trans = BM_claim_trans;
                    BM_trans_prob(BM_lowest_index) = ...
                        BM_trans_prob(BM_lowest_index) ...
                        + 1 - sum(BM_trans_prob);
                else
                    % otherwise, create a new case corresponding to the
                    % lowest possible Bonus-Malus transition
                    BM_all_trans = [BM_lowest; BM_claim_trans];
                    BM_trans_prob = [1 - sum(BM_trans_prob); ...
                        BM_trans_prob]; %#ok<AGROW>
                end
                
                claim_amount_ranges = claim_ranges ...
                    - BM.deductible(BM_id, T - t + 1);
                mean_claim_amount = sum(compound_specexp_approx(atoms, ...
                    P_cp(:, best_Miti), zeros(size(claim_amount_ranges, ...
                    1), 1), claim_amount_ranges, ...
                    BM.deductible(BM_id, T - t + 1), ...
                    BM.cap(BM_id, T - t + 1)));
                
                transition{T - t + 1}{BM_id, Ins_id} = ...
                    struct('BM', BM_all_trans, ...
                    'Ins', 2 * ones(length(BM_all_trans), 1), ...
                    'prob', BM_trans_prob, ...
                    'no_loss_prob', no_loss_prob, ...
                    'no_claim_prob', no_claim_prob, ...
                    'mean_claim', mean_claim_amount);
            else
                % if the optimal control is to not activate insurance, the
                % transition is deterministic (according to the rule for
                % inactive contract)
                transition{T - t + 1}{BM_id, Ins_id} = ...
                    struct('BM', BM_up, 'Ins', wait_up + 2, 'prob', 1);
            end
        end
    end
end

cost = costs{end};

init = struct('BM', BM.init, 'Ins', 1);

output.transition = transition;

occupancy_before = cell(T + 1, 1);
BMstats = zeros(T, BM_num + 1);
mitistats = zeros(T, Miti_num);
claimstats = zeros(T, 4);
mean_mitigation = 0;
mean_premium = 0;
mean_penalty = 0;
mean_loss = 0;
mean_claim = 0;
mean_prevented = 0;

if ~LDA.lognorm_sev
    unprevented_exp = trunc_g_and_h_uppexp(0, LDA.g, LDA.h, ...
        LDA.loc, LDA.sca, LDA.tol) * params.scale * LDA.freq_mean;
else
    unprevented_exp = lognorm_uppexp(0, LDA.mu, LDA.sig2) ...
        * params.scale * LDA.freq_mean;
end

occupancy_before{1} = zeros(BM_num, BM_wait_max + 2);
occupancy_before{1}(init.BM, init.Ins) = 1;

for t = 1:T
    occupancy_before{t + 1} = zeros(BM_num, BM_wait_max + 2);
    
    for BM_id = 1:BM_num
        for Ins_id = 1:BM_wait_max + 2
            occupancy_prob_t = occupancy_before{t}(BM_id, Ins_id);
            if occupancy_prob_t == 0
                continue;
            end
            
            trans_t = transition{t}{BM_id, Ins_id};
            BM_list = trans_t.BM;
            Ins_list = trans_t.Ins;
            prob_list = trans_t.prob;
            
            assert(length(BM_list) == length(Ins_list) ...
                && length(BM_list) == length(prob_list), ...
                'the transition cell arrary is invalid');
            
            % compute occupancy probability
            for j = 1:length(BM_list)
                occupancy_before{t + 1}(BM_list(j), Ins_list(j)) = ...
                    occupancy_before{t + 1}(BM_list(j), Ins_list(j)) ...
                    + occupancy_prob_t * prob_list(j);
            end
            
            % compute mitigation statistics
            policy_p = policy{t, 1}{BM_id, Ins_id};
            if ~isempty(policy_p)
                mitistats(t, policy_p.Miti) = mitistats(t, ...
                    policy_p.Miti) + occupancy_prob_t;
            end
            
            % compute claim statistics
            if isfield(trans_t, 'no_claim_prob')
                claimstats(t, 3) = claimstats(t, 3) ...
                    + occupancy_prob_t ...
                    * (trans_t.no_claim_prob - trans_t.no_loss_prob);
                claimstats(t, 1) = claimstats(t, 1) ...
                    + occupancy_prob_t * trans_t.no_loss_prob;
                claimstats(t, 2) = claimstats(t, 2) ...
                    + occupancy_prob_t * (1 - trans_t.no_claim_prob);
            end
            
            % compute expected values
            mean_mitigation = mean_mitigation + occupancy_prob_t ...
                * discount ^ (t - 1) * Miti.cost(policy_p.Miti);
            
            if policy_p.Ins == 1
                mean_premium = mean_premium + occupancy_prob_t ...
                    * discount ^ (t - 1) * BM.premium(BM_id);
                
                if Ins_id == 1
                    % joining for the first time
                    mean_penalty = mean_penalty + occupancy_prob_t ...
                        * discount ^ (t - 1) * BM.penalty_in(t);
                end
                
                if Ins_id > 2
                    % rejoining
                    mean_penalty = mean_penalty + occupancy_prob_t ...
                        * discount ^ (t - 1) * BM.penalty_rejoin;
                end
                
                mean_claim = mean_claim + occupancy_prob_t ...
                    * discount ^ (t - 1) ...
                    * transition{t}{BM_id, Ins_id}.mean_claim;
                BMstats(t, BM_id + 1) = ...
                    BMstats(t, BM_id + 1) + occupancy_prob_t;
            else
                BMstats(t, 1) = BMstats(t, 1) ...
                    + occupancy_prob_t;
                
                if Ins_id == 2
                    % withdrawing
                    mean_penalty = mean_penalty + occupancy_prob_t ...
                        * discount ^ (t - 1) * BM.penalty_out(t);
                end
            end
            
            mean_loss = mean_loss + occupancy_prob_t ...
                * discount ^ (t - 1) * upper_exp(policy_p.Miti);
            mean_prevented = mean_prevented + occupancy_prob_t ...
                * discount ^ (t - 1) * (unprevented_exp ...
                - upper_exp(policy_p.Miti));
        end
    end
end

claimstats(:, 4) = max(1 - sum(claimstats(:, 1:3), 2), 0);

output.costs = costs;
output.occupancy = occupancy_before;
output.BMstats = BMstats;
output.mitistats = mitistats;
output.claimstats = claimstats;
output.mean = struct;
output.mean.mitigation = mean_mitigation;
output.mean.premium = mean_premium;
output.mean.penalty = mean_penalty;
output.mean.loss = mean_loss;
output.mean.claim = mean_claim;
output.mean.prevented = mean_prevented;

output.risk_premium = risk_premium;
output.saved_compound = saved_compound;

end

