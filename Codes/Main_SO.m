clc; close all; clear all;

U_pen = 0.008;
% Data_quan: day-ahead price + distritbution of (real-time minor day-ahead) in Jan.
% Data_quan_20: real-time plus normal noise with sigma = 20

% RC of multi-zone thermal dynamics
% Based on which building topology
% C*ds/dt = (s^out-s)/R_i + sum(s-s)/R_ij + pi + s_in

rng(42, 'twister');
N_days = 10;
Z_num = 100;
T = 288;

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

C0 = readmatrix('C.csv'); % C: 2-5 J/K
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
    if R0(rr,rr) > 2.0
        R0(rr,rr) = R0(rr,rr)*1.8;
    else
        R0(rr,rr) = R0(rr,rr)*1.8;
    end
end
C_lambda = 1.0;
C0_mean = mean(C0);
C0 = C_lambda * C0 + (1 - C_lambda) * C0_mean;

C = repmat(C0, 5, 1);
R = 10000*ones(Z_num,Z_num);
for ii = 1:5
    R((ii-1)*20+1:ii*20,(ii-1)*20+1:ii*20) = R0;
end
C = 3*C;
R = 2*R;

COP = 5;

beta_ = COP./C; %2.2;
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
S0 = 23*ones(Z_num,1);
%Sref = [23, repmat(23, 1, 12 * 8), repmat(23.5, 1, 12 * 4), repmat(24, 1, 12 * 4), repmat(23, 1, 12 * 8)];
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

P_min = zeros(Z_num,T);  %KW
P_max = 4*ones(Z_num,T); %KW

T_min = 20*ones(Z_num,T);
T_max = 26*ones(Z_num,T);
zone_group1 = [1,13,14,6,7,8,15,16,17]; %meeting room and office room
zone_group2 = [4,5,3,9,10,18,19];%print room and lounge
zone_group3 = [2,11,12,20];%corridor and finance room

for ii = 1:5
    T_min(20*(ii-1)+zone_group1,:) = 21;
    T_max(20*(ii-1)+zone_group1,:) = 25;
    T_min(20*(ii-1)+zone_group2,:) = 20.5;
    T_max(20*(ii-1)+zone_group2,:) = 25.5;
    T_min(20*(ii-1)+zone_group3,:) = 20;
    T_max(20*(ii-1)+zone_group3,:) = 26;
end

T_ref_min = ones(Z_num,T);
T_ref_max = ones(Z_num,T);
for tt = 1:T
    T_ref_min(:,tt) = T_min(:,tt) - Sref(tt+1);
    T_ref_max(:,tt) = T_max(:,tt) - Sref(tt+1);
end


% Graph and partitioning
% G = (B, E)
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
rng(42, 'twister');
alpha_ij_sigma = sqrt(0.2 * alpha_ij);
alpha_ij_sampling = alpha_ij + alpha_ij_sigma.*randn(size(alpha_ij,1), 100);

num_clusters = 8;
[b_result,e_result] = Zone_cluster(E,num_clusters,alpha_ij);

% Build equivalent single-zone model
[Eta_,Beta_,Delta_new_,P_c_min,P_c_max,T_c_ref_min,T_c_ref_max] = Zone_aggre(num_clusters,T,N_days,C,R_1,b_result,eta_,beta_,delta_new,P_max,T_ref_min,T_ref_max,COP);


%% Online SO - MPC
rng(42, 'twister');
N_scen = 5;
H = 24;

mu_0 = 0;
sigma_0 = 0.01;

SO_Cost = zeros(N_days,3);
SO_Time = zeros(N_days,1);

for day_i = 7:10
    quantiles_all = load('Data_quan_Final_10.mat');
    Price_samss = [20, 50, 100];
    Price_sam_i = 1;
    Price_sam = Price_samss(Price_sam_i);
    quantile_y = quantiles_all.(['y_' num2str(day_i-1) '_' num2str(Price_sam)])/1000; 
    
    scen_price = zeros(T, N_scen);
    scen_error = zeros(5*num_clusters, T, N_scen);
    
    for s = 1:N_scen
        idx = randi(size(quantile_y, 2));
        scen_price(:, s) = quantile_y(:, idx);
        scen_error(:, :, s) = mu_0 + sigma_0 * randn(5*num_clusters, T);
    end
    
    v_real_100z = zeros(Z_num, 1);
    all_p_real_100z = zeros(Z_num, T);
    all_tem_real_100z = zeros(Z_num, T+1);
    
    data_in = load('Data_10days.mat');
    Pri_rt_Feb = data_in.Pri_rt_Feb/1000;
    Pri_rt_online = Pri_rt_Feb((day_i-1)*T+1 : day_i*T, 1);
    real_total_cost_ele = 0;
    tic;
    for t_step = 1:T
        if mod(t_step, 24) == 0
            fprintf('Solving Step %d/%d...\n', t_step, T);
        end
        % (Mapping 100 Zones -> Clusters) ---
        v_agg_feedback = zeros(5*num_clusters, 1);
        for cc = 1:num_clusters
            for cci = 1:5
                cluster_idx = (cci-1)*num_clusters + cc;
                zone_list = (cci-1)*20 + find(b_result(:, cc) >= 0.99);
                v_agg_feedback(cluster_idx) = sum(C(zone_list).*v_real_100z(zone_list)) / sum(C(zone_list));
            end
        end
    
        % --- 2. (MPC) ---
        current_H = min(H, T - t_step + 1);
        v_c_p = sdpvar(5*num_clusters, current_H, 'full');
        v_c_tem = cell(N_scen, 1);
        for s = 1:N_scen
            v_c_tem{s} = sdpvar(5*num_clusters, current_H + 1, 'full');
        end
    
        Constraints = [];
        Objective = 0;
        for s = 1:N_scen
            Constraints = [Constraints, v_c_tem{s}(:, 1) == v_agg_feedback];
            for k = 1:current_H
                Constraints = [Constraints, v_c_tem{s}(:, k+1) == ...
                    Eta_ .* v_c_tem{s}(:, k) + Beta_ .* v_c_p(:, k) + ...
                    reshape(Delta_new_(t_step + k - 1, day_i, :), [], 1) + scen_error(:, t_step + k - 1, s)];
            end
            Constraints = [Constraints, P_c_min(:, t_step:t_step+current_H-1) <= v_c_p <= P_c_max(:, t_step:t_step+current_H-1)];
            Constraints = [Constraints, T_c_ref_min(:, t_step:t_step+current_H-1) <= v_c_tem{s}(:, 2:end) <= T_c_ref_max(:, t_step:t_step+current_H-1)];
            
            cost_ele_s = (1/12) * scen_price(t_step:t_step+current_H-1, s)' * sum(v_c_p, 1)';
            cost_comfort_s = U_pen * sum(v_c_tem{s}(:, 2:end).^2, 'all');
            Objective = Objective + (1/N_scen) * (cost_ele_s + cost_comfort_s);
        end
        optimize(Constraints, Objective, sdpsettings('solver','gurobi','verbose',0));
        
        opt_p_agg_t = value(v_c_p(:, 1));
        
        [p_real_t, v_next_100z] = Online_Step_Evolution(Z_num, num_clusters, b_result, opt_p_agg_t, ...
                                    v_real_100z, T_ref_max(:, t_step), T_ref_min(:, t_step), P_max(:, t_step), ...
                                    day_i, t_step, eta_, beta_, delta_new, alpha_);
        
        all_p_real_100z(:, t_step) = p_real_t;
        v_real_100z = v_next_100z; 
        all_tem_real_100z(:, t_step+1) = v_real_100z;
        
        real_total_cost_ele = real_total_cost_ele + (1/12) * Pri_rt_online(t_step) * sum(p_real_t);
    end
    time_elapsed_online = toc;

    real_total_cost_comfort = U_pen * sum(all_tem_real_100z(:, 2:end).^2, 'all');
    fprintf('Online SO Real Ele Cost: %.4f\n', real_total_cost_ele);
    fprintf('Online SO Real Comfort Cost: %.4f\n', real_total_cost_comfort);
    fprintf('Sum: %.4f\n', real_total_cost_ele + real_total_cost_comfort);
    SO_Cost(day_i,1) = real_total_cost_ele;
    SO_Cost(day_i,2) = real_total_cost_comfort;
    SO_Cost(day_i,3) = real_total_cost_ele + real_total_cost_comfort;
    SO_Time(day_i,1) = time_elapsed_online; %/(5*num_clusters); can not divide by zone number due to its central optimization
    if day_i == 4
        aa = all_tem_real_100z(:,2:end)';
        save('Tem_SO_new.mat','aa');
    end
end

%%
function [p_real_t, v_next_100z] = Online_Step_Evolution(Z_num, num_clusters, b_result, p_agg_t, v_curr_100z, T_ref_max_t, T_ref_min_t, P_max_t, day_i, t, eta_, beta_, delta_new, alpha_)
    p_real_t = zeros(Z_num, 1);
    for cc = 1:num_clusters
        for cci = 1:5
            cluster_idx_in_40 = (cci-1)*num_clusters + cc;
            if num_clusters == 1
                zone_list = (cci-1)*20 + (1:20);
            else
                zone_list = (cci-1)*20 + find(b_result(:, cc) >= 0.99);
            end
            
            n_sub = length(zone_list);
            p_cluster_target = p_agg_t(cluster_idx_in_40);
            
            if n_sub == 1
                p_real_t(zone_list) = p_cluster_target;
            else
                % 2. : min sum( (s_{i,t+1})^2 )
                % s_{i,t+1} = beta_i * P_i + K_i
                % K_i = eta_i*s_{i,t} + delta_i + sum(alpha_ij * s_{j,t})
                current_temps = v_curr_100z(zone_list, :);
                current_temps_ = v_curr_100z;
                K = eta_(zone_list) .* current_temps + ...
                    reshape(delta_new(t, day_i, zone_list), [], 1) + ...
                    alpha_(zone_list, :) * current_temps_  + 0.1*randn(n_sub,1);
                
                B = beta_(zone_list);
    
                % sum( (B_i*P_i + K_i)^2 ) = P'*(diag(B^2))*P + (2*B.*K)'*P + const
                % 
                model.Q = sparse(diag(B.^2));      
                model.obj = 2 * (B .* K);          
                
                % sum(P_i) = P_cluster
                model.A = sparse(ones(1, n_sub));
                model.rhs = p_cluster_target;
                model.sense = '=';
                
                % 
                model.lb = zeros(n_sub,1);
                model.ub = P_max_t(zone_list,1);
                
                % 
                params.OutputFlag = 0;
                results = gurobi(model, params);
                
                p_real_t(zone_list) = results.x;
            end
        end
    end
    % 
    v_next_100z = eta_ .* v_curr_100z + beta_ .* p_real_t + ...
                  reshape(delta_new(t, day_i, :), [], 1) + alpha_ * v_curr_100z;
    v_next_100z = max(v_next_100z, T_ref_min_t);
    v_next_100z = min(v_next_100z, T_ref_max_t);
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
    Eta_ = ones(5*num_clusters,1);
    Beta_ = ones(5*num_clusters,1);
    Delta_new_ = zeros(T, N_days, 5*num_clusters);
    P_c_min = zeros(5*num_clusters, T);
    P_c_max = zeros(5*num_clusters, T);
    T_c_ref_min = zeros(5*num_clusters, T);
    T_c_ref_max = zeros(5*num_clusters, T);
    
    C_1 = C(1:20,1);
    for cc = 1:num_clusters
        if num_clusters == 1
            zone_list = 1:20;
        else
            zone_list = find(b_result(:, cc) >= 0.99);
        end
        if isscalar(zone_list)
            % case1：only one zone in a cluster
            for cci = 1:5
                Eta_((cci-1)*num_clusters+cc,1) = eta_((cci-1)*20+zone_list,1);
                Beta_((cci-1)*num_clusters+cc,1) = beta_((cci-1)*20+zone_list,1);
                Delta_new_(:,:,(cci-1)*num_clusters+cc) = delta_new(:,:,(cci-1)*20+zone_list);
                P_c_max((cci-1)*num_clusters+cc,:) = P_max((cci-1)*20+zone_list,:);
                T_c_ref_min((cci-1)*num_clusters+cc,:) = T_ref_min((cci-1)*20+zone_list,:);
                T_c_ref_max((cci-1)*num_clusters+cc,:) = T_ref_max((cci-1)*20+zone_list,:);
            end
        else
            % case2：multiple zones in a cluster
            for cci = 1:5
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