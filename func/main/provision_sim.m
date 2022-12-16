function [costs, output] = provision_sim(BM, Miti, LDA, params, init, ...
    policy, sim_num, additional_stats, parallel)
%PROVISION_SIM Simulation of the cyber security provisioning process
% Inputs: 
%       BM: same as in optimal_provision
%       Miti: same as in optimal_provision
%       LDA: same as in optimal_provision
%       params: same as in optimal_provision
%       init: initial state (struct)
%       policy: same as in the output of optimal_provision
%       sim_num: number of simulations
%       additional_stats: boolean indicating whether additional statistics
%           are needed (default is true)
%       parallel: boolean indicating whether simulations should be run in
%           parallel (default is false). If so then additional_stats =
%           false
% Outputs: 
%       costs: the actual cost of each simulation
%       output: struct containing additional outputs
%           output.occupancy: the empirical occupancy frequencies
%           output.BMstats: the table of statistics about Bonus-Malus
%               levels
%           output.mitistats: the table of statistics about mitigation
%           output.claimstats: the table of statistics about claim
%           output.details: struct containing other detailed distributions

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

if ~isfield(params, 'discount')
    params.discount = 1;
end

if ~isfield(params, 'scale')
    params.scale = 1;
end

if ~exist('additional_stats', 'var') || isempty(additional_stats)
    additional_stats = true;
end

if ~exist('parallel', 'var') || isempty(parallel)
    parallel = false;
end

if parallel
    additional_stats = false;
end

T = params.horizon;
discount = params.discount;
comp_scale = params.scale;

BM_num = length(BM.premium);
assert(BM_num == size(BM.rule, 1) && BM_num == size(BM.rule, 2), ...
    'mis-specified Bonus-Malus');
assert(BM_num == size(BM.inactive_rule, 1) ...
    && 2 == size(BM.inactive_rule, 2), ...
    'mis-specified Bonue-Malus inactive rules');
BM_wait_max = max(max(BM.inactive_rule(:, 1)), 1);
Miti_num = length(Miti.cost);
assert(Miti_num == length(Miti.effect), 'mis-specified mitigation');

assert(T == size(policy, 1) && 2 == size(policy, 2), ...
    'mis-specified policy');

if ~LDA.lognorm_sev
    sev_gen = @(n)(trunc_g_and_h_rand(n, LDA.g, LDA.h, ...
        LDA.loc, LDA.sca, LDA.tol));
else
    sev_gen = @(n)(lognorm_rand(n, LDA.mu, LDA.sig2));
end

costs = zeros(sim_num, 1);

if additional_stats
    occupancy = cell(T + 1, 1);
    BMstats = zeros(T, BM_num + 1);
    mitistats = zeros(T, Miti_num);
    claimstats = zeros(T, 4);
    details = struct;
    details.mitigation = zeros(sim_num, 1);
    details.premium = zeros(sim_num, 1);
    details.penalty = zeros(sim_num, 1);
    details.loss = zeros(sim_num, 1);
    details.claim = zeros(sim_num, 1);
    details.prevented = zeros(sim_num, 1);
    
    for t = 1:T + 1
        occupancy{t} = zeros(BM_num, BM_wait_max + 2);
    end
end

% pre-generate losses
LDA_samples = cell(Miti_num, 1);
LDA_prevented = cell(Miti_num, 1);
for Miti_id = 1:Miti_num
    [comp_samp, prev_samp] = compound_nbin_mitigation_rand(sim_num * T, ...
        LDA.freq_mean, LDA.freq_var, sev_gen, Miti.effect(Miti_id));
    LDA_samples{Miti_id} = reshape(comp_samp * comp_scale, sim_num, T);
    LDA_prevented{Miti_id} = reshape(prev_samp * comp_scale, sim_num, T);
end

if ~parallel
    for sim_id = 1:sim_num
        % initialize the state
        state = init;
        cost = 0;
        cur_discount = 1;
        
        if additional_stats
            % increment the occupancy counter
            occupancy{1}(state.BM, state.Ins) = occupancy{1}(state.BM, ...
                state.Ins) + 1;
        end
        
        % iterate for each time period
        for t = 1:T
            % Provision Stage
            
            % select control for the provision stage
            cont_p = policy{t, 1}{state.BM, state.Ins};
            
            % randomly generate loss
            loss_samp = LDA_samples{cont_p.Miti}(sim_id, t);
            prev_samp = LDA_prevented{cont_p.Miti}(sim_id, t);
            
            % compute the cost
            cost = cost + cur_discount * (Miti.cost(cont_p.Miti) ...
                + (cont_p.Ins == 1) * BM.premium(state.BM) + loss_samp);
            
            penalty = 0;
            
            if cont_p.Ins == 1 && state.Ins ~= 2
                penalty = penalty + BM.penalty_in(t);
            end
            
            if cont_p.Ins == 2 && state.Ins == 2
                penalty = penalty + BM.penalty_out(t);
            end
            
            if cont_p.Ins == 1 && state.Ins > 2
                penalty = penalty + BM.penalty_rejoin;
            end
            
            cost = cost + cur_discount * penalty;
            
            if additional_stats
                % update mitigation statistics
                mitistats(t, cont_p.Miti) = mitistats(t, cont_p.Miti) + 1;
                
                % cost break-down
                details.mitigation(sim_id) = details.mitigation(sim_id) ...
                    + cur_discount * Miti.cost(cont_p.Miti);
                details.premium(sim_id) = details.premium(sim_id) ...
                    + cur_discount * (cont_p.Ins == 1) ...
                    * BM.premium(state.BM);
                details.penalty(sim_id) = details.penalty(sim_id) ...
                    + cur_discount * penalty;
                details.loss(sim_id) = details.loss(sim_id) ...
                    + cur_discount * loss_samp;
                details.prevented(sim_id) = details.prevented(sim_id) ...
                    + cur_discount * prev_samp;
            end
            
            
            % update the insurance status for the claim stage
            if cont_p.Ins == 1
                % if the optimal control is to activate insurance
                state.Ins = 2;
            elseif state.Ins == 2
                % if the optimal control is to deactivate insurance and the
                % insurance was previously active
                state.Ins = -1;
            end
            state.Loss = loss_samp;
            
            if additional_stats
                if state.Ins == 2
                    BMstats(t, state.BM + 1) = BMstats(t, state.BM + 1) ...
                        + 1;
                else
                    BMstats(t, 1) = BMstats(t, 1) + 1;
                end
            end
            
            % Claim Stage
            if state.Ins == 2
                % if the insurance is active (otherwise there is no 
                % control)
                % select control for the claim stage
                cont_c = policy{t, 2}{state.BM};
                if ~isempty(cont_c.range) ...
                        && any(state.Loss >= cont_c.range(:, 1) ...
                        & state.Loss <= cont_c.range(:, 2))
                    % choose to make a claim
                    claim = min(max(loss_samp ...
                        - BM.deductible(state.BM, t), 0), ...
                        BM.cap(state.BM, t));
                    
                    if additional_stats
                        claimstats(t, 2) = claimstats(t, 2) + 1;
                    end
                else
                    % choose not to make a claim
                    claim = 0;
                    if additional_stats
                        if loss_samp == 0
                            claimstats(t, 1) = claimstats(t, 1) + 1;
                        else
                            claimstats(t, 3) = claimstats(t, 3) + 1;
                        end
                    end
                end
                
                if additional_stats
                    details.claim(sim_id) = details.claim(sim_id) ...
                        + cur_discount * claim;
                end
                
                % compute the cost
                cost = cost - cur_discount * claim;
                
                % update the state
                state.BM = bonusmalus_transition(BM.rule, state.BM, claim);
                state.Loss = 0;
            else
                if state.Ins == -1
                    % if the insurance is freshly deactivated
                    wait_cur = 0;
                else
                    wait_cur = state.Ins - 2;
                end
                % if the insurance has been deactivated before
                [BM_up, wait_up] = bonusmalus_forget(BM.inactive_rule, ...
                    state.BM, wait_cur);
                state.BM = BM_up;
                state.Ins = wait_up + 2;
                state.Loss = 0;
                
                if additional_stats
                    claimstats(t, 4) = claimstats(t, 4) + 1;
                end
            end
            
            if additional_stats
                % increment the occupancy counter
                occupancy{t + 1}(state.BM, state.Ins) = ...
                    occupancy{t + 1}(state.BM, state.Ins) + 1;
            end
            
            % update the discount factor
            cur_discount = cur_discount * discount;
        end
        
        costs(sim_id) = cost;
    end
else
    parfor sim_id = 1:sim_num
        % initialize the state
        state = init;
        cost = 0;
        cur_discount = 1;
        
        % iterate for each time period
        for t = 1:T
            % Provision Stage
            
            % select control for the provision stage
            cont_p = policy{t, 1}{state.BM, state.Ins}; %#ok<PFBNS> 
            
            % randomly generate loss
            loss_samp = LDA_samples{cont_p.Miti}(sim_id, t); %#ok<PFBNS> 
            
            % compute the cost
            cost = cost + cur_discount * (Miti.cost(cont_p.Miti) ...
                + (cont_p.Ins == 1) * BM.premium(state.BM) ...
                + loss_samp); %#ok<PFBNS> 
            
            penalty = 0;
            
            if cont_p.Ins == 1 && state.Ins ~= 2
                penalty = penalty + BM.penalty_in(t);
            end
            
            if cont_p.Ins == 2 && state.Ins == 2
                penalty = penalty + BM.penalty_out(t);
            end
            
            if cont_p.Ins == 1 && state.Ins > 2
                penalty = penalty + BM.penalty_rejoin;
            end
            
            cost = cost + cur_discount * penalty;
            
            % update the insurance status for the claim stage
            if cont_p.Ins == 1
                % if the optimal control is to activate insurance
                state.Ins = 2;
            elseif state.Ins == 2
                % if the optimal control is to deactivate insurance and the
                % insurance was previously active
                state.Ins = -1;
            end
            state.Loss = loss_samp;
            
            % Claim Stage
            if state.Ins == 2
                % if the insurance is active (otherwise there is no control)
                % select control for the claim stage
                cont_c = policy{t, 2}{state.BM};
                if ~isempty(cont_c.range) ...
                        && any(state.Loss >= cont_c.range(:, 1) ...
                        & state.Loss <= cont_c.range(:, 2))
                    % choose to make a claim
                    claim = min(max(loss_samp - BM.deductible(state.BM, ...
                        t), 0), BM.cap(state.BM, t));
                else
                    % choose not to make a claim
                    claim = 0;
                end
                
                % compute the cost
                cost = cost - cur_discount * claim;
                
                % update the state
                state.BM = bonusmalus_transition(BM.rule, state.BM, claim);
                state.Loss = 0;
            else
                if state.Ins == -1
                    % if the insurance is freshly deactivated
                    wait_cur = 0;
                else
                    wait_cur = state.Ins - 2;
                end
                % if the insurance has been deactivated before
                [BM_up, wait_up] = bonusmalus_forget(BM.inactive_rule, ...
                    state.BM, wait_cur);
                state.BM = BM_up;
                state.Ins = wait_up + 2;
                state.Loss = 0;
            end
            
            % update the discount factor
            cur_discount = cur_discount * discount;
        end
        
        costs(sim_id) = cost;
    end
end

if additional_stats
    % normalize the occupancy
    for t = 1:T + 1
        occupancy{t} = occupancy{t} / sim_num;
    end

    BMstats = BMstats / sim_num;
    mitistats = mitistats / sim_num;
    claimstats = claimstats / sim_num;

    output.occupancy = occupancy;
    output.BMstats = BMstats;
    output.mitistats = mitistats;
    output.claimstats = claimstats;
    output.details = details;
else
    output = [];
end
    
end

