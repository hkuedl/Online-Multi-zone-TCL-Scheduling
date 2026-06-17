clc; close all; clear all;

U_pen = 0.008;

% Four robustness cases for RC-based multi-zone thermal dynamics.
% Strong coupling: smaller inter-zone thermal resistance R_ij.
% Weak coupling: larger inter-zone thermal resistance R_ij.
case_list = struct([]);
case_list(1).name = 'SC_200_strong_coupling_Rij0p8_20zone_10floor'; case_list(1).Rij_scale = 0.8; case_list(1).floor_num = 10; case_list(1).seed = 21;
case_list(2).name = 'SC_100_strong_coupling_Rij0p8_20zone_5floor';  case_list(2).Rij_scale = 0.8; case_list(2).floor_num = 5;  case_list(2).seed = 22;
case_list(3).name = 'WC_200_weak_coupling_Rij1p2_20zone_10floor';   case_list(3).Rij_scale = 1.2; case_list(3).floor_num = 10; case_list(3).seed = 23;
case_list(4).name = 'WC_100_weak_coupling_Rij1p2_20zone_5floor';    case_list(4).Rij_scale = 1.2; case_list(4).floor_num = 5;  case_list(4).seed = 24;
for kk = 1:numel(case_list), case_list(kk).noise_level = 0.05; end

N_days = 10;
T = 288;
zones_per_floor = 20;
num_clusters = 8;              % fixed cluster number; no cluster-number sensitivity analysis
% eval_days = 1;               % uncomment this line for a quick smoke test
eval_days = 1:N_days;          % default: compare GT and proposed method over all sampled days

output_root = 'Scalability';
if ~exist(output_root, 'dir'), mkdir(output_root); end
all_case_results = struct([]);
summary_rows = {};

for case_idx = 1:numel(case_list)
    case_cfg = case_list(case_idx);
    rng(case_cfg.seed, 'twister');
    fprintf('\n================ Running case %d/%d: %s ================\n', case_idx, numel(case_list), case_cfg.name);
    floor_num = case_cfg.floor_num;
    Z_num = zones_per_floor * floor_num;

    data_in = readmatrix('90zone_15min.csv');
    T_o0 = data_in(2:end, 2);
    T_occ0 = data_in(2:end, 3:92)/1000;

    T_o1 = zeros(T, 28);
    T_occ1 = zeros(T, 28, Z_num);
    original_indices = 0:95;
    new_indices = linspace(0, 95, T);
    for day_i = 0:28-1
        for z_i = 1:Z_num
            T_occ1(:,day_i+1,z_i) = interp1(original_indices, T_occ0(day_i*96+1:(day_i+1)*96, mod(z_i-1,90)+1), new_indices, 'linear')';
        end
        T_o1(:,day_i+1) = interp1(original_indices, T_o0(day_i*96+1:(day_i+1)*96, 1), new_indices, 'linear')';
    end

    C0 = readmatrix('C.csv');
    R0 = readmatrix('R.csv');
    R0(isnan(R0)) = 10000;
    for ii = 2:20
        for jj = 1:ii-1
            R0(ii,jj) = R0(jj,ii);
        end
    end

    R0_non_diag_mask = ~eye(size(R0, 1));
    R0_mask_less = (R0 <= 2.3) & R0_non_diag_mask;
    R0_mask_greater = (R0 > 2.3) & R0_non_diag_mask;
    R0(R0_mask_less) = R0(R0_mask_less)*0.8;
    R0(R0_mask_greater) = R0(R0_mask_greater)*2.5;
    for rr = 1:20
        R0(rr,rr) = R0(rr,rr)*1.8;
    end
    C_lambda = 1.0;
    C0_mean = mean(C0);
    C0 = C_lambda * C0 + (1 - C_lambda) * C0_mean;

    C = repmat(C0, floor_num, 1);
    R = 10000*ones(Z_num,Z_num);
    R_block_case = R0;
    Rij_mask = (~eye(zones_per_floor)) & (R_block_case < 1000);
    rng(case_cfg.seed, 'twister');
    Rij_noise = 1 + case_cfg.noise_level * randn(zones_per_floor, zones_per_floor);
    Rij_noise = (Rij_noise + Rij_noise') / 2;
    Rij_noise = max(Rij_noise, 0.10);
    R_block_case(Rij_mask) = R_block_case(Rij_mask) .* case_cfg.Rij_scale .* Rij_noise(Rij_mask);
    for ii = 1:floor_num
        idx = (ii-1)*zones_per_floor+1:ii*zones_per_floor;
        R(idx,idx) = R_block_case;
    end
    C = 3*C;
    R = 2*R;

    COP = 5;
    beta_ = COP./C;
    eta_ = zeros(Z_num,1);
    alpha_ = zeros(Z_num,Z_num);
    for z_i1 = 1:Z_num
        eta_(z_i1,1) = 1;
        for z_i2 = 1:Z_num
            eta_(z_i1,1) = eta_(z_i1,1) - 1/(R(z_i1,z_i2)*C(z_i1,1));
            if z_i2 ~= z_i1
                alpha_(z_i1,z_i2) = 1/(R(z_i1,z_i2)*C(z_i1,1));
            end
        end
    end

    delta_ = zeros(T, 28, Z_num);
    delta_new = zeros(T, 28, Z_num);
    Sref = [23, repmat(23, 1, 12 * 8), repmat(23, 1, 12 * 4), repmat(23, 1, 12 * 4), repmat(23, 1, 12 * 8)];
    for day_i = 0:28-1
        for z_i = 1:Z_num
            delta_(:,day_i+1,z_i) = T_o1(:,day_i+1)/(R(z_i,z_i)*C(z_i,1)) + T_occ1(:,day_i+1,z_i)/C(z_i,1);
            for t = 1:T
                delta_new(t,day_i+1,z_i) = delta_(t,day_i+1,z_i) - Sref(t+1) + eta_(z_i,1)*Sref(t) + sum(alpha_(z_i,:)*Sref(t));
            end
        end
    end
    delta_new = delta_new(:,[11, 12, 14, 15, 16, 18, 19, 20, 21, 22],:);

    P_min = zeros(Z_num,T);
    P_max = 4*ones(Z_num,T);
    T_min = 20*ones(Z_num,T);
    T_max = 26*ones(Z_num,T);
    zone_group1 = [1,13,14,6,7,8,15,16,17];
    zone_group2 = [4,5,3,9,10,18,19];
    zone_group3 = [2,11,12,20];
    for ii = 1:floor_num
        T_min(20*(ii-1)+zone_group1,:) = 21;   T_max(20*(ii-1)+zone_group1,:) = 25;
        T_min(20*(ii-1)+zone_group2,:) = 20.5; T_max(20*(ii-1)+zone_group2,:) = 25.5;
        T_min(20*(ii-1)+zone_group3,:) = 20;   T_max(20*(ii-1)+zone_group3,:) = 26;
    end
    T_ref_min = ones(Z_num,T);
    T_ref_max = ones(Z_num,T);
    for tt = 1:T
        T_ref_min(:,tt) = T_min(:,tt) - Sref(tt+1);
        T_ref_max(:,tt) = T_max(:,tt) - Sref(tt+1);
    end

    R_1 = R(1:20,1:20);
    alpha_1 = alpha_(1:20,1:20);
    mask_upper_triangle = logical(triu(ones(size(R_1)), 1));
    mask_value_lt_100 = (R_1 < 1000);
    final_mask = mask_upper_triangle & mask_value_lt_100;
    [rows, cols] = find(final_mask);
    E = [rows, cols];
    alpha_ij = ones(size(rows,1),1);
    for ii = 1:size(rows,1)
        alpha_ij(ii,1) = 0.5*(alpha_1(rows(ii,1),cols(ii,1)) + alpha_1(cols(ii,1),rows(ii,1)));
    end
    [b_result,e_result] = Zone_cluster(E,num_clusters,alpha_ij);
    [Eta_,Beta_,Delta_new_,P_c_min,P_c_max,T_c_ref_min,T_c_ref_max] = Zone_aggre(num_clusters,T,N_days,C,R_1,b_result,eta_,beta_,delta_new,P_max,T_ref_min,T_ref_max,COP);

    %% Ground-truth vs proposed comparison; no DP and no cluster-number sensitivity
    price_data = load('Data.mat');
    Pri_rt_Feb = price_data.Pri_rt_Feb/1000;
    ops = sdpsettings('solver','gurobi','verbose',0);
    num_equiv_zones = floor_num * num_clusters;

    cost.GT.ele = zeros(numel(eval_days),1);
    cost.GT.comfort = zeros(numel(eval_days),1);
    cost.GT.total = zeros(numel(eval_days),1);
    cost.Proposed.model_ele = zeros(numel(eval_days),1);
    cost.Proposed.model_comfort = zeros(numel(eval_days),1);
    cost.Proposed.model_total = zeros(numel(eval_days),1);
    cost.Proposed.real_ele = zeros(numel(eval_days),1);
    cost.Proposed.real_comfort = zeros(numel(eval_days),1);
    cost.Proposed.real_total = zeros(numel(eval_days),1);
    solution_time.GT = zeros(numel(eval_days),1);
    solution_time.Proposed = zeros(numel(eval_days),1);
    diagnostics_GT = cell(numel(eval_days),1);
    diagnostics_Proposed = cell(numel(eval_days),1);
    proposed_solution = struct([]);

    for dd = 1:numel(eval_days)
        day_i = eval_days(dd);
        fprintf('\n--- Case %s | day %d ---\n', case_cfg.name, day_i);
        price_vec = Pri_rt_Feb((day_i-1)*T+1:day_i*T,1);

        % Ground-truth optimization on the original coupled model
        v_p = sdpvar(Z_num,T,'full');
        v_tem = sdpvar(Z_num,T+1,'full');
        cost_ele = (1/12)*sum(price_vec.*sum(v_p,1)');
        cost_comfort = U_pen * sum(v_tem(:,2:end).^2, 'all');
        Objective = cost_ele + cost_comfort;
        Constraints = [v_tem(:,1) == 0, P_min <= v_p, P_max >= v_p, T_ref_min <= v_tem(:,2:end), T_ref_max >= v_tem(:,2:end)];
        Constraints_model = [];
        for t = 1:T
            Constraints_model = [Constraints_model, v_tem(:,t+1) == eta_.*v_tem(:,t) + beta_.*v_p(:,t) + reshape(delta_new(t,day_i,:), [], 1) + alpha_(:,:)*v_tem(:,t)];
        end
        tic;
        diagnostics = optimize([Constraints,Constraints_model], Objective, ops);
        solution_time.GT(dd) = toc;
        diagnostics_GT{dd} = diagnostics;
        GT_v_p_opt = value(v_p);
        GT_v_tem_opt = value(v_tem);
        cost.GT.ele(dd) = value(cost_ele);
        cost.GT.comfort(dd) = value(cost_comfort);
        cost.GT.total(dd) = value(Objective);

        % Proposed method: fixed clustering + equivalent-zone optimization + zone-level recovery
        v_c_p = sdpvar(num_equiv_zones,T,'full');
        v_c_tem = sdpvar(num_equiv_zones,T+1,'full');
        cost_c_ele = (1/12)*sum(price_vec.*sum(v_c_p,1)');
        cost_c_comfort = U_pen * sum(v_c_tem(:,2:end).^2, 'all');
        Objective_c = cost_c_ele + cost_c_comfort;
        Constraints_c = [v_c_tem(:,1) == 0, P_c_min <= v_c_p, P_c_max >= v_c_p, T_c_ref_min <= v_c_tem(:,2:end), T_c_ref_max >= v_c_tem(:,2:end)];
        Constraints_c_model = [];
        for t = 1:T
            Constraints_c_model = [Constraints_c_model, v_c_tem(:,t+1) == Eta_.*v_c_tem(:,t) + Beta_.*v_c_p(:,t) + reshape(Delta_new_(t,day_i,:), [], 1)];
        end
        tic;
        diagnostics_c = optimize([Constraints_c,Constraints_c_model], Objective_c, ops);
        solution_time.Proposed(dd) = toc;
        diagnostics_Proposed{dd} = diagnostics_c;
        PR_v_c_p_opt = value(v_c_p);
        PR_v_c_tem_opt = value(v_c_tem);
        cost.Proposed.model_ele(dd) = value(cost_c_ele);
        cost.Proposed.model_comfort(dd) = value(cost_c_comfort);
        cost.Proposed.model_total(dd) = value(Objective_c);

        [PR_v_p_final_proposed,PR_v_tem_final_proposed,~,~] = Zone_recovery(Z_num,num_clusters,b_result,PR_v_c_p_opt,T_ref_max,T_ref_min,P_max,day_i,eta_,beta_,delta_new,alpha_,PR_v_c_p_opt(:,1),1);
        cost.Proposed.real_ele(dd) = (1/12) * sum(price_vec .* sum(PR_v_p_final_proposed, 1)');
        cost.Proposed.real_comfort(dd) = U_pen * sum(PR_v_tem_final_proposed(:, 2:end).^2, 'all');
        cost.Proposed.real_total(dd) = cost.Proposed.real_ele(dd) + cost.Proposed.real_comfort(dd);

        proposed_solution(dd).day = day_i;
        proposed_solution(dd).cluster_power = PR_v_c_p_opt;
        proposed_solution(dd).cluster_temperature = PR_v_c_tem_opt;
        proposed_solution(dd).recovered_power = PR_v_p_final_proposed;
        proposed_solution(dd).recovered_temperature = PR_v_tem_final_proposed;

        fprintf('GT:       ele = %.6f, comfort = %.6f, total = %.6f, time = %.3fs\n', cost.GT.ele(dd), cost.GT.comfort(dd), cost.GT.total(dd), solution_time.GT(dd));
        fprintf('Proposed: ele = %.6f, comfort = %.6f, total = %.6f, time = %.3fs\n', cost.Proposed.real_ele(dd), cost.Proposed.real_comfort(dd), cost.Proposed.real_total(dd), solution_time.Proposed(dd));
    end

    metrics.cost_gap_percent = 100*(cost.Proposed.real_total - cost.GT.total)./cost.GT.total;
    metrics.ele_gap_percent = 100*(cost.Proposed.real_ele - cost.GT.ele)./cost.GT.ele;
    metrics.comfort_gap_percent = 100*(cost.Proposed.real_comfort - cost.GT.comfort)./cost.GT.comfort;
    metrics.speedup_GT_over_Proposed = solution_time.GT ./ solution_time.Proposed;
    metrics.mean_GT_total = mean(cost.GT.total);
    metrics.mean_Proposed_total = mean(cost.Proposed.real_total);
    metrics.mean_cost_gap_percent = mean(metrics.cost_gap_percent);
    metrics.mean_GT_ele = mean(cost.GT.ele);
    metrics.mean_Proposed_ele = mean(cost.Proposed.real_ele);
    metrics.mean_GT_comfort = mean(cost.GT.comfort);
    metrics.mean_Proposed_comfort = mean(cost.Proposed.real_comfort);
    metrics.mean_GT_time = mean(solution_time.GT);
    metrics.mean_Proposed_time = mean(solution_time.Proposed);
    metrics.mean_speedup_GT_over_Proposed = mean(metrics.speedup_GT_over_Proposed);

    all_case_results(case_idx).name = case_cfg.name;
    all_case_results(case_idx).Rij_scale = case_cfg.Rij_scale;
    all_case_results(case_idx).floor_num = floor_num;
    all_case_results(case_idx).Z_num = Z_num;
    all_case_results(case_idx).num_clusters = num_clusters;
    all_case_results(case_idx).eval_days = eval_days;
    all_case_results(case_idx).cost = cost;
    all_case_results(case_idx).metrics = metrics;
    all_case_results(case_idx).solution_time = solution_time;

    summary_rows(end+1,:) = {case_cfg.name, floor_num, Z_num, case_cfg.Rij_scale, metrics.mean_GT_ele, metrics.mean_GT_comfort, metrics.mean_GT_total, metrics.mean_Proposed_ele, metrics.mean_Proposed_comfort, metrics.mean_Proposed_total, metrics.mean_cost_gap_percent, metrics.mean_GT_time, metrics.mean_Proposed_time, metrics.mean_speedup_GT_over_Proposed};
end

summary_table = cell2table(summary_rows, 'VariableNames', {'case_name','floor_num','Z_num','Rij_scale','GT_ele','GT_comfort','GT_total','Proposed_ele','Proposed_comfort','Proposed_total','cost_gap_percent','GT_time','Proposed_time','speedup_GT_over_Proposed'});
save(fullfile(output_root, 'all_case_results.mat'), 'all_case_results', 'case_list', 'summary_table');
disp(summary_table);

%%
load(fullfile('Scalability','all_case_results.mat'));

case_titles = {'SC-200','SC-100','WC-200','WC-100'};

figure('Position',[100,100,1000,650]);

for i = 1:numel(all_case_results)
    subplot(2,2,i);
    
    days = all_case_results(i).eval_days;
    GT_total = all_case_results(i).cost.GT.total;
    PR_total = all_case_results(i).cost.Proposed.real_total;
    gap_percent = 100 * (PR_total - GT_total) ./ GT_total;
    mean_gap = mean(gap_percent);
    
    plot(days, GT_total, '-o', 'LineWidth', 1.8, 'MarkerSize', 5); hold on;
    plot(days, PR_total, '--s', 'LineWidth', 1.8, 'MarkerSize', 5);
    
    grid on; box on;
    xlabel('Day');
    ylabel('Total cost ($)');
    title(case_titles{i}, 'FontWeight', 'bold');
    
    legend({'GT','Proposed'}, 'Location', 'best');
    
    %text(0.05, 0.92, sprintf('Mean gap = %.2f%%', mean_gap), ...
    %    'Units', 'normalized', 'FontSize', 10);
    
    xlim([min(days), max(days)]);
end

%sgtitle('Daily Total Cost Comparison between GT and Proposed Method');

set(findall(gcf,'-property','FontName'),'FontName','Arial');
set(findall(gcf,'-property','FontSize'),'FontSize',11);

exportgraphics(gcf, fullfile('Scalability','New_fig_scalability.pdf'), 'ContentType', 'vector');

%%

load(fullfile('Scalability','all_case_results.mat'));

case_short_names = {'SC-200','SC-100','WC-200','WC-100'};

GT_mean = zeros(1,4);
GT_std = zeros(1,4);
PR_mean = zeros(1,4);
PR_std = zeros(1,4);

for i = 1:4
    GT_total = all_case_results(i).cost.GT.total;
    PR_total = all_case_results(i).cost.Proposed.real_total;
    
    GT_mean(i) = mean(GT_total);
    GT_std(i) = std(GT_total);
    PR_mean(i) = mean(PR_total);
    PR_std(i) = std(PR_total);
end

fprintf('\\begin{table}[t]\n');
fprintf('\\centering\n');
fprintf('\\caption{Cost comparison between GT and Proposed methods under four building cases. The cost values are reported as the 10-day mean $\\pm$ standard deviation.}\n');
fprintf('\\label{tab:cost_comparison_robustness}\n');
fprintf('\\begin{tabular}{lcccc}\n');
fprintf('\\hline\n');
fprintf('Method & SC-200 & SC-100 & WC-200 & WC-100 \\\\\n');
fprintf('\\hline\n');

fprintf('GT ');
for i = 1:4
    fprintf('& %.4f $\\pm$ %.4f ', GT_mean(i), GT_std(i));
end
fprintf('\\\\\n');

fprintf('Proposed ');
for i = 1:4
    fprintf('& %.4f $\\pm$ %.4f ', PR_mean(i), PR_std(i));
end
fprintf('\\\\\n');

fprintf('\\hline\n');
fprintf('\\end{tabular}\n');
fprintf('\\end{table}\n');

%% Functions
function [PR_v_p_final_proposed,PR_v_tem_final_proposed,toc1,toc2] = Zone_recovery(Z_num,num_clusters,b_result,PR_v_c_p_opt,T_ref_max,T_ref_min,P_max,day_i,eta_,beta_,delta_new_in,alpha_,tem_ini,t_now)
    rng(42, 'twister');
    zones_per_floor = 20;
    floor_num = Z_num / zones_per_floor;
    T = size(PR_v_c_p_opt,2);
    PR_v_p_final_proposed = zeros(Z_num, T);
    PR_v_tem_final_proposed = zeros(Z_num, T+1);
    for t = 1:T
        tic1 = tic;
        for cc = 1:num_clusters
            for cci = 1:floor_num
                cluster_idx_in_40 = (cci-1)*num_clusters + cc;
                if num_clusters == 1
                    zone_list = (cci-1)*20 + (1:20);
                else
                    zone_list = (cci-1)*20 + find(b_result(:, cc) >= 0.99);
                end
                
                n_sub = length(zone_list);
                p_cluster_target = PR_v_c_p_opt(cluster_idx_in_40, t);
                
                if n_sub == 1
                    PR_v_p_final_proposed(zone_list, t) = p_cluster_target;
                else
                    % 2. QP: min sum( (s_{i,t+1})^2 )
                    % s_{i,t+1} = beta_i * P_i + K_i
                    % K_i = eta_i*s_{i,t} + delta_i + sum(alpha_ij * s_{j,t})
                    if T >= 2
                        current_temps = PR_v_tem_final_proposed(zone_list, t);
                        current_temps_ = PR_v_tem_final_proposed(:, t);
                        K = eta_(zone_list) .* current_temps + ...
                            reshape(delta_new_in(t, day_i, zone_list), [], 1) + ...
                            alpha_(zone_list, :) * current_temps_  + 0.1*randn(n_sub,1);
                    else
                        current_temps = tem_ini(zone_list, t);
                        current_temps_ = tem_ini(:, t);
                    
                        K = eta_(zone_list) .* current_temps + ...
                            reshape(delta_new_in(t_now, day_i, zone_list), [], 1) + ...
                            alpha_(zone_list, :) * current_temps_  + 0.1*randn(n_sub,1);
                    end
                    B = beta_(zone_list);

                    % sum( (B_i*P_i + K_i)^2 ) = P'*(diag(B^2))*P + (2*B.*K)'*P + const
                    
                    model.Q = sparse(diag(B.^2));     
                    model.obj = 2 * (B .* K);         
                    
                    % sum(P_i) = P_cluster
                    model.A = sparse(ones(1, n_sub));
                    model.rhs = p_cluster_target;
                    model.sense = '=';

                    model.lb = zeros(n_sub, 1);
                    model.ub = P_max(zone_list, t);

                    params.OutputFlag = 0;
                    results = gurobi(model, params);
                    
                    PR_v_p_final_proposed(zone_list, t) = results.x;
                end
            end
        end
        toc1 = toc(tic1);
        tic2 = tic;
     
        if T >= 2
            PR_v_tem_final_proposed(:, t+1) = eta_ .* PR_v_tem_final_proposed(:, t) + ...
                                       beta_ .* PR_v_p_final_proposed(:, t) + ...
                                       reshape(delta_new_in(t, day_i, :), [], 1) + ...
                                       alpha_(:, :) * PR_v_tem_final_proposed(:, t);
        else
            a1 = eta_ .* tem_ini;
            a2 = beta_ .* PR_v_p_final_proposed(:, t);
            a3 = reshape(delta_new_in(t_now, day_i, :), [], 1);
            a4 = alpha_(:, :) * tem_ini;
            PR_v_tem_final_proposed(:, t+1) = a1 + ...
                                       a2 + ...
                                       a3 + ...
                                       a4;
            PR_v_tem_final_proposed(:, t+1) = max(PR_v_tem_final_proposed(:, t+1), T_ref_min(:, t));
            PR_v_tem_final_proposed(:, t+1) = min(PR_v_tem_final_proposed(:, t+1), T_ref_max(:, t));
        end
        toc2 = toc(tic2);
    end
end


function [b_result,e_result] = Zone_cluster(E,num_clusters,alpha_ij)
    num_nodes = 20;
    B = 1:num_nodes;
    
    num_edges = size(E, 1);
    
    Clu = 1:num_clusters;
    
    fprintf("Node number: %d\n", num_nodes);
    fprintf("Edge number: %d\n", num_edges);
    fprintf("Cluster number: %d\n", num_clusters);
    
    b = binvar(num_nodes, num_clusters, 'full');   % b_ic
    r = binvar(num_edges, num_clusters, 'full');   % r_ijc
    e = binvar(num_edges, 1);                      % e_ij
    
    Constraints = [];
    
    % a
    for i = B
        Constraints = [Constraints, sum(b(i, :)) == 1];
    end
    
    % b
    for c = Clu
        Constraints = [Constraints, sum(b(:, c)) >= 1];
    end
    
    for k = 1:num_edges
        i = E(k, 1);
        j = E(k, 2);
        %f
        Constraints = [Constraints, e(k) == sum(r(k, :))];
        
        for c = Clu
            %c
            Constraints = [Constraints, r(k, c) <= b(i, c)];
            %d
            Constraints = [Constraints, r(k, c) <= b(j, c)];
            %e
            Constraints = [Constraints, b(i, c) + b(j, c) <= r(k, c) + 1];
        end
    end
    
    min_cluster_size = 1;
    if num_clusters == 1
        max_cluster_size = 20;
    elseif num_clusters == 2
        max_cluster_size = 10;
    else
        max_cluster_size = 7;
    end
    for c = Clu
        Constraints = [Constraints, sum(b(:,c)) >= min_cluster_size];
        Constraints = [Constraints, sum(b(:,c)) <= max_cluster_size];
    end
    
    omega = 0.5;
    Objective = sum((1-e).*alpha_ij) - omega*sum(sum(r.*repmat(alpha_ij,1,num_clusters)));
    ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
    sol = optimize(Constraints, Objective, ops);
    fprintf('Obj: %.2f\n\n', value(Objective));
    b_result = value(b);
    e_result = value(e);
    for c = 1:num_clusters
        if num_clusters == 1
            cluster_assignment = 1:20;
        else
            cluster_assignment = find(b_result(:, c) >= 0.99);
        end
        fprintf('cluster %d: zones %s \n', c, num2str(cluster_assignment(:)'));
    end
end

function [Eta_,Beta_,Delta_new_,P_c_min,P_c_max,T_c_ref_min,T_c_ref_max] = Zone_aggre(num_clusters,T,N_days,C,R_1,b_result,eta_,beta_,delta_new,P_max,T_ref_min,T_ref_max,COP)
    zones_per_floor = 20;
    floor_num = size(C,1) / zones_per_floor;
    Eta_ = ones(floor_num*num_clusters,1);
    Beta_ = ones(floor_num*num_clusters,1);
    Delta_new_ = zeros(T, N_days, floor_num*num_clusters);
    P_c_min = zeros(floor_num*num_clusters, T);
    P_c_max = zeros(floor_num*num_clusters, T);
    T_c_ref_min = zeros(floor_num*num_clusters, T);
    T_c_ref_max = zeros(floor_num*num_clusters, T);
    
    C_1 = C(1:20,1);
    for cc = 1:num_clusters
        if num_clusters == 1
            zone_list = 1:20;
        else
            zone_list = find(b_result(:, cc) >= 0.99);
        end
        if isscalar(zone_list)
            % case1：only one zone in a cluster
            for cci = 1:floor_num
                Eta_((cci-1)*num_clusters+cc,1) = eta_((cci-1)*20+zone_list,1);
                Beta_((cci-1)*num_clusters+cc,1) = beta_((cci-1)*20+zone_list,1);
                Delta_new_(:,:,(cci-1)*num_clusters+cc) = delta_new(:,:,(cci-1)*20+zone_list);
                P_c_max((cci-1)*num_clusters+cc,:) = P_max((cci-1)*20+zone_list,:);
                T_c_ref_min((cci-1)*num_clusters+cc,:) = T_ref_min((cci-1)*20+zone_list,:);
                T_c_ref_max((cci-1)*num_clusters+cc,:) = T_ref_max((cci-1)*20+zone_list,:);
            end
        else
            % case2：multiple zones in a cluster
            for cci = 1:floor_num
                zone_C = C_1(zone_list);
                C_sum = sum(zone_C);
                Beta_((cci-1)*num_clusters+cc,1) = COP/C_sum;
                
                G_outdoor_sum = 0;
                for i = 1:length(zone_list)
                    z_idx = zone_list(i);
                    G_outdoor_sum = G_outdoor_sum + 1 / R_1(z_idx, z_idx);
                end
                Eta_((cci-1)*num_clusters+cc,1) = 1 - G_outdoor_sum/C_sum;
                
                weighted_delta = bsxfun(@times, delta_new(:,:,(cci-1)*20+zone_list), reshape(zone_C, 1, 1, []));
                Delta_new_(:,:,(cci-1)*num_clusters+cc) = sum(weighted_delta, 3) / C_sum;
        
                P_c_max((cci-1)*num_clusters+cc,:) = sum(P_max((cci-1)*20+zone_list,:), 1);
                T_c_ref_min((cci-1)*num_clusters+cc,:) = sum(C((cci-1)*20+zone_list,1).*T_ref_min((cci-1)*20+zone_list,:))/sum(C((cci-1)*20+zone_list,1));
                T_c_ref_max((cci-1)*num_clusters+cc,:) = sum(C((cci-1)*20+zone_list,1).*T_ref_max((cci-1)*20+zone_list,:))/sum(C((cci-1)*20+zone_list,1));
            end
        end
    end
end