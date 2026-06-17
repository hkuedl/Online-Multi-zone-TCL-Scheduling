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
%%
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

%% Sensitivity analysis of cluster numbers
% Ground Truth (GT) and Decoupled (DE) results
day_i = 1;

data_in = load('Data.mat');
Pri_rt_Feb = data_in.Pri_rt_Feb/1000; %%Unit: $/MWh

v_p = sdpvar(Z_num,T,'full');
v_tem = sdpvar(Z_num,T+1,'full');

cost_ele = (1/12)*sum(Pri_rt_Feb((day_i-1)*T+1:day_i*T,1).*sum(v_p,1)');
cost_comfort = U_pen * sum(v_tem(:,2:end).^2, 'all');

Objective = cost_ele + cost_comfort;

Constraints = [];
Constraints = [Constraints, v_tem(:,1) == 0];

Constraints_model = [];
for t = 1:T
    Constraints_model = [Constraints_model, v_tem(:,t+1) == eta_.*v_tem(:,t) + beta_.*v_p(:,t) + reshape(delta_new(t,day_i,:), [], 1) + alpha_(:,:)*v_tem(:,t)];
end

Constraints = [Constraints, P_min <= v_p];
Constraints = [Constraints, P_max >= v_p];
Constraints = [Constraints, T_ref_min <= v_tem(:,2:end)];
Constraints = [Constraints, T_ref_max >= v_tem(:,2:end)];

ops = sdpsettings('solver','gurobi','verbose',0);
diagnostics = optimize([Constraints,Constraints_model], Objective, ops);

GT_v_p_opt = value(v_p);
GT_v_tem_opt = value(v_tem);
GT_v_cost_total = value(Objective);
GT_v_cost_ele = value(cost_ele);
GT_v_cost_comfort = value(cost_comfort);
% disp('------------');
% disp(['ele cost_Ground truth: ', num2str(GT_v_cost_ele)]);
% disp(['ele comfort_Ground truth: ', num2str(GT_v_cost_comfort)]);
% disp(['total cost_Ground truth: ', num2str(GT_v_cost_total)]);

Constraints_model_decouple = [];
for t = 1:T
    Constraints_model_decouple = [Constraints_model_decouple, v_tem(:,t+1) == eta_.*v_tem(:,t) + beta_.*v_p(:,t) + reshape(delta_new(t,day_i,:), [], 1)];
end

diagnostics_decouple = optimize([Constraints,Constraints_model_decouple], Objective, ops);

DE_v_p_opt = value(v_p);
DE_v_tem_opt = value(v_tem);
DE_v_tem_opt_real = zeros(Z_num,T+1);
for t = 1:T
    DE_v_tem_opt_real(:, t+1) = eta_.*DE_v_tem_opt_real(:, t) + beta_.*DE_v_p_opt(:, t) + reshape(delta_new(t,day_i,:), [], 1) + alpha_(:,:)*DE_v_tem_opt_real(:, t);
end

DE_v_cost_total = value(Objective);
DE_v_cost_ele = value(cost_ele);
DE_v_cost_comfort = value(cost_comfort);
DE_v_cost_comfort_real = U_pen * sum(DE_v_tem_opt_real(:,2:end).^2, 'all');
DE_v_cost_total_real = DE_v_cost_ele + DE_v_cost_comfort_real;
% disp(['ele cost (Decouple): ', num2str(DE_v_cost_ele)]);
% disp(['ele comfort (Decouple): ', num2str(DE_v_cost_comfort)]);
% disp(['total cost (Decouple): ', num2str(DE_v_cost_total)]);
% disp(['ele comfort_real (Decouple): ', num2str(DE_v_cost_comfort_real)]);
% disp(['total cost_real (Decouple): ', num2str(DE_v_cost_total_real)]);

PR_v_cost_total_real_sum = ones(20,1);
for num_clusters_i = 1:20
    [b_result,e_result] = Zone_cluster(E,num_clusters_i,alpha_ij);
    [Eta_,Beta_,Delta_new_,P_c_min,P_c_max,T_c_ref_min_i,T_c_ref_max_i] = Zone_aggre(num_clusters_i,T,N_days,C,R_1,b_result,eta_,beta_,delta_new,P_max,T_ref_min,T_ref_max,COP);
    
    % Proposed method (PR)
    v_c_p = sdpvar(5*num_clusters_i,T,'full');
    v_c_tem = sdpvar(5*num_clusters_i,T+1,'full');
    
    cost_c_ele = (1/12)*sum(Pri_rt_Feb((day_i-1)*T+1:day_i*T,1).*sum(v_c_p,1)');
    cost_c_comfort = U_pen * sum(v_c_tem(:,2:end).^2, 'all');
    
    Objective_c = cost_c_ele + cost_c_comfort;
    
    Constraints_c = [];
    Constraints_c = [Constraints_c, v_c_tem(:,1) == 0];
    
    Constraints_c_model = [];
    for t = 1:T
        Constraints_c_model = [Constraints_c_model, v_c_tem(:,t+1) == Eta_.*v_c_tem(:,t) + Beta_.*v_c_p(:,t) + reshape(Delta_new_(t,day_i,:), [], 1)];
    end
    
    Constraints_c = [Constraints_c, P_c_min <= v_c_p];
    Constraints_c = [Constraints_c, P_c_max >= v_c_p];
    Constraints_c = [Constraints_c, T_c_ref_min_i <= v_c_tem(:,2:end)];
    Constraints_c = [Constraints_c, T_c_ref_max_i >= v_c_tem(:,2:end)];
    
    ops = sdpsettings('solver','gurobi','verbose',0);
    diagnostics = optimize([Constraints_c,Constraints_c_model], Objective_c, ops);
    
    PR_v_c_p_opt = value(v_c_p);
    PR_v_c_tem_opt = value(v_c_tem);
    PR_v_c_cost_total = value(Objective_c);
    PR_v_c_cost_ele = value(cost_c_ele);
    PR_v_c_cost_comfort = value(cost_c_comfort);
    fprintf('\n--- Final Comparison (Real Performance) ---\n'); 
    disp(['ele cost: (Proposed)', num2str(PR_v_c_cost_ele)]);
    disp(['ele comfort: (Proposed)', num2str(PR_v_c_cost_comfort)]);
    disp(['total cost: (Proposed)', num2str(PR_v_c_cost_total)]);
    
    % Power Recovery and Real Performance Calculation (Proposed)
    [PR_v_p_final_proposed,PR_v_tem_final_proposed,~,~] = Zone_recovery(Z_num,num_clusters_i,b_result,PR_v_c_p_opt,T_ref_max,T_ref_min,P_max,day_i,eta_,beta_,delta_new,alpha_,PR_v_c_p_opt(:,1),1);
    
    PR_v_cost_ele_real = (1/12) * sum(Pri_rt_Feb((day_i-1)*T+1:day_i*T, 1) .* sum(PR_v_p_final_proposed, 1)');
    
    PR_v_cost_comfort_real = U_pen * sum(PR_v_tem_final_proposed(:, 2:end).^2, 'all');
    
    PR_v_cost_total_real = PR_v_cost_ele_real + PR_v_cost_comfort_real;

    PR_v_cost_total_real_sum(num_clusters_i) = PR_v_cost_total_real;
    
    disp(['ele cost_real (Partitioning): ', num2str(PR_v_cost_ele_real)]);
    disp(['comfort cost_real (Partitioning): ', num2str(PR_v_cost_comfort_real)]);
    disp(['Total cost_real (Partitioning): ', num2str(PR_v_cost_total_real)]);
    disp(['Total cost_real (Decouple): ', num2str(DE_v_cost_total_real)]);
    disp(['Total cost (Ground-truth): ', num2str(GT_v_cost_total)]);
end

disp('To sum up');
disp(['Total cost_real (Decouple): ', num2str(DE_v_cost_total_real)]);
disp(['Total cost (Ground-truth): ', num2str(GT_v_cost_total)]);
disp(PR_v_cost_total_real_sum);

%% Two-stage deterministic optimization (First stage)

gamma = 0.05;  %1-gamma = 95%
mu_0 = 0;
sigma_0 = 0.01; 

num_equiv_zones = 5 * num_clusters;
T_c_ref_min_CC = zeros(num_equiv_zones, T); 
T_c_ref_max_CC = zeros(num_equiv_zones, T); 

z_score = norminv(1 - gamma, mu_0, sigma_0);

for t = 1:T
    xi_over_t = mu_0 + t * sigma_0 * z_score;  %  \overline{\xi}_{i,t}
    xi_under_t = mu_0 - t * sigma_0 * z_score; %  \underline{\xi}_{i,t}
    
    T_c_ref_min_CC(:, t) = T_c_ref_min(:, t) - xi_under_t; 
    T_c_ref_max_CC(:, t) = T_c_ref_max(:, t) - xi_over_t;
end

%% Two-stage deterministic optimization (Second stage)
rng(42, 'twister');

Price_samss = [20,50,100];
Price_sam_i = 1;
Price_sam = Price_samss(Price_sam_i);
data_in = load('Data_10days.mat');
Pri_rt_Feb = data_in.Pri_rt_Feb/1000;

% Proposed method: using percentiles to represent the price uncertainty
quantiles_all = load('Data_quan_Final_10.mat');
s_sam_number = 801;
quantiles_x = quantiles_all.(['x_' num2str(Price_sam)])(1,2:end-1);
s_sam = -4 + 0.01*(0:s_sam_number-1)';
u = zeros(s_sam_number, 2);
u(:,1) = s_sam;
u(:,2) = U_pen*2*u(:,1);

% offline: calculate value function
v = zeros(N_days, T, num_equiv_zones, s_sam_number);
v_3 = zeros(N_days, T, num_equiv_zones, s_sam_number);

% online scheduling
DP_Cost = zeros(N_days,3);
DP_Time = zeros(N_days,1);

DP_v = zeros(N_days,T,num_equiv_zones);
DP_tem = zeros(N_days,T,num_equiv_zones);
DP_tem_recover = zeros(N_days,T,Z_num);
DP_v_recover = zeros(N_days,T,Z_num);
DP_tem_aggre = zeros(N_days,T,num_equiv_zones);

for opt_day = 1:10
    disp(['days: ',num2str(opt_day)]);
    for zone = 1:num_equiv_zones
        for t = T:-1:1
            if t == T
                v(opt_day,t,zone,:) = zeros(1, 1, s_sam_number);
            else
                quantile_y = quantiles_all.(['y_' num2str(opt_day-1) '_' num2str(Price_sam)])/1000;
                quantiles_rt = quantile_y(t+1,:);

                for ii = 1:s_sam_number
                    [v(opt_day,t,zone,ii), w_x1, w_x2, v_3(opt_day,t,zone,ii)] = ...
                        F_value(s_sam, quantiles_x, quantiles_rt, s_sam(ii), ...
                               Eta_(zone,1), Beta_(zone,1), u, squeeze(v(opt_day,t+1,zone,:)), P_c_max(zone,t), ...
                               Delta_new_(t,opt_day,zone));
                end
            end
        end
        if mod(zone , 10) == 0
            disp(['Finish ', num2str(zone)]);
        end
    end

    time_elapsed_online = zeros(T,1);
    for t = 1:288
        tic1 = tic;
        for zone = 1:num_equiv_zones
            if t == 1
                tem_ini = 0;
            else
                tem_ini = DP_tem_aggre(opt_day,t-1,zone);
                %tem_ini = DP_tem(opt_day,t-1,zone);
                %tem_ini = 0.3*DP_tem_aggre(opt_day,t-1,zone) + 0.7*DP_tem(opt_day,t-1,zone);
            end

            v_dt = squeeze(v(opt_day,t,zone,:));
            v_dt_0 = diff(v_dt);
            v_dt1 = [v_dt(1); v_dt_0];

            price_rt = Pri_rt_Feb((opt_day-1)*288+t,1);
            %s_max = Eta_(zone,1)*tem_ini + Beta_(zone,1)*P_c_max(zone,t) + Delta_new_(t,opt_day,zone);
            %s_min = Eta_(zone,1)*tem_ini + Beta_(zone,1)*P_c_min(zone,t) + Delta_new_(t,opt_day,zone);
            %[~, min_ind] = min(abs(s_sam - s_min));
            %[~, max_ind] = min(abs(s_sam - s_max));
            
            s_phys_max = Eta_(zone,1)*tem_ini + Beta_(zone,1)*P_c_max(zone,t) + Delta_new_(t,opt_day,zone);
            s_phys_min = Eta_(zone,1)*tem_ini + Beta_(zone,1)*P_c_min(zone,t) + Delta_new_(t,opt_day,zone);
            s_constr_min = T_c_ref_min_CC(zone,t);
            s_constr_max = T_c_ref_max_CC(zone,t);
            
            s_valid_min = max(s_phys_min, s_constr_min);
            s_valid_max = min(s_phys_max, s_constr_max);
            
            if s_valid_min > s_valid_max
                disp('Wrong !')
            end
            [~, min_ind] = min(abs(s_sam - s_valid_min));
            [~, max_ind] = min(abs(s_sam - s_valid_max));

            s = s_sam(min_ind:max_ind);
            p = (s - Eta_(zone,1)*tem_ini - Delta_new_(t,opt_day,zone)) / Beta_(zone,1);

            if p(1) < P_c_min(zone,t)
                s = s(2:end);
                p = p(2:end);
            end
            if p(end) > P_c_max(zone,t)
                s = s(1:end-1);
                p = p(1:end-1);
            end

            c = zeros(length(s), 4);
            for i = 1:length(s)
                c(i,1) = (1/12)*price_rt*p(i);
                c(i,2) = U_pen*(s(i)^2);
                c(i,3) = V_cal(v_dt1, s_sam, s(i));
                c(i,4) = sum(c(i,1:3));
            end
            
            [~, idx] = min(c(:,4));
            DP_v(opt_day,t,zone) = p(idx);
            DP_tem(opt_day,t,zone) = s(idx);
        end
        toc1 = toc(tic1);
        
        if t == 1
            DP_tem_ini = zeros(Z_num,1);
        else
            DP_tem_ini = squeeze(DP_tem_recover(opt_day,t-1,:));
        end

        [DP_v_recover_t,DP_tem_recover_t,tt1,tt2] = Zone_recovery(Z_num,num_clusters,b_result,reshape(DP_v(opt_day,t,:), [], 1),T_ref_max,T_ref_min,P_max,opt_day,eta_,beta_,delta_new,alpha_,DP_tem_ini,t);

        tic2 = tic;
        DP_v_recover(opt_day,t,:) = DP_v_recover_t(:,1);
        DP_tem_recover(opt_day,t,:) = DP_tem_recover_t(:,2:end);
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
                    DP_tem_aggre(opt_day,t,(cci-1)*num_clusters+cc) = DP_tem_recover(opt_day,t,(cci-1)*20+zone_list);
                    %DP_tem_aggre(opt_day,t,(cci-1)*num_clusters+cc) = min(DP_tem_aggre(opt_day,t,(cci-1)*num_clusters+cc), T_c_ref_max_CC((cci-1)*num_clusters+cc,t));
                    %DP_tem_aggre(opt_day,t,(cci-1)*num_clusters+cc) = max(DP_tem_aggre(opt_day,t,(cci-1)*num_clusters+cc), T_c_ref_min_CC((cci-1)*num_clusters+cc,t));
                end
            else
                for cci = 1:5
                    zone_C = C_1(zone_list);
                    C_sum = sum(zone_C);
                    DP_tem_aggre(opt_day,t,(cci-1)*num_clusters+cc) = sum(C_1(zone_list,1).*reshape(DP_tem_recover(opt_day,t,(cci-1)*20+zone_list), [], 1))/C_sum;
                    %DP_tem_aggre(opt_day,t,(cci-1)*num_clusters+cc) = min(DP_tem_aggre(opt_day,t,(cci-1)*num_clusters+cc), T_c_ref_max_CC((cci-1)*num_clusters+cc,t));
                    %DP_tem_aggre(opt_day,t,(cci-1)*num_clusters+cc) = max(DP_tem_aggre(opt_day,t,(cci-1)*num_clusters+cc), T_c_ref_min_CC((cci-1)*num_clusters+cc,t));
                end
            end
        end
        toc2 = toc(tic2);
        time_elapsed_online(t,1) = (toc1+tt1)/num_equiv_zones + tt2 + toc2;
    end
    
    DP_Cost(opt_day,1) = (1/12) * sum(Pri_rt_Feb((opt_day-1)*T+1:opt_day*T, 1) .* reshape(sum(DP_v_recover(opt_day,:,:), 3),[],1));
    DP_Cost(opt_day,2) = U_pen * sum(DP_tem_recover(opt_day, :, :).^2, 'all');
    DP_Cost(opt_day,3) = DP_Cost(opt_day,1) + DP_Cost(opt_day,2);
    disp(['ele cost_real (Proposed): ', num2str(DP_Cost(opt_day,1))]);
    disp(['comfort cost_real (Proposed): ', num2str(DP_Cost(opt_day,2))]);
    disp(['Total cost_real (Proposed): ', num2str(DP_Cost(opt_day,3))]);
    DP_Time(opt_day,1) = sum(time_elapsed_online);
    %if opt_day == 4
    %    aa = squeeze(DP_tem_recover(opt_day,:,:));
    %    save('Tem_PR.mat','aa');
    %end
end

%% Functions
function [F_val, w_x1, w_x2, num_filtered] = F_value(s_sam, quantiles_x, quantiles_rt, s, eta_, beta_, u, v, Pmax, delta)
    x1 = eta_*s+ beta_*Pmax + delta;
    x2 = eta_*s + delta;

    [w_x1, w_x2] = deal(v_u(s_sam, v, u, x1), v_u(s_sam, v, u, x2));

    P1 = eta_*w_x1*F_price(quantiles_x, quantiles_rt, -12*beta_*w_x1); % modify
    P2 = eta_*w_x2*(1 - F_price(quantiles_x, quantiles_rt, -12*beta_*w_x2)); % modify

    weights0 = diff(quantiles_x);
    weights = [weights0(1), weights0];

    mask = (quantiles_rt >= -12*beta_*w_x1) & (quantiles_rt <= -12*beta_*w_x2); % modify
    filtered_values = quantiles_rt(mask);
    filtered_weights = weights(mask);

    P3 = (eta_/(12*beta_))*sum(filtered_values .* filtered_weights);
    F_val = P1 + P2 - P3; % modify
    num_filtered = length(filtered_values);
end

function p = F_price(quantiles_x, quantiles_rt, x)
    if x < quantiles_rt(1)
        p = 0;
    elseif x > quantiles_rt(end)
        p = 1;
    elseif x == quantiles_rt(end)
        p = quantiles_x(end);
    else
        for i = 1:length(quantiles_x)-1
            if x >= quantiles_rt(i) && x < quantiles_rt(i+1)
                p = quantiles_x(i);
                return;
            end
        end
    end
end

function diff_vu = v_u(s_sam, v, u, x)
    if x < u(1,1)
        index_u = 1;
    elseif x > u(end,1)
        index_u = size(u,1);
    else
        [~, index_u] = min(abs(u(:,1) - x));
    end
    ux = u(index_u,2);

    if x < s_sam(1)
        index_v = 1;
    elseif x > s_sam(end)
        index_v = length(s_sam);
    else
        [~, index_v] = min(abs(s_sam - x));
    end
    vx = v(index_v);

    diff_vu = vx + ux; % modify
end

function Vx = V_cal(v_dt1, s_sam, x)
    Vx = 0;
    for i = 1:length(s_sam)
        Vx = Vx + v_dt1(i)*max(x - s_sam(i), 0);
    end
end

function [PR_v_p_final_proposed,PR_v_tem_final_proposed,toc1,toc2] = Zone_recovery(Z_num,num_clusters,b_result,PR_v_c_p_opt,T_ref_max,T_ref_min,P_max,day_i,eta_,beta_,delta_new_in,alpha_,tem_ini,t_now)
    rng(42, 'twister');
    T = size(PR_v_c_p_opt,2);
    PR_v_p_final_proposed = zeros(Z_num, T);
    PR_v_tem_final_proposed = zeros(Z_num, T+1);
    for t = 1:T
        tic1 = tic;
        for cc = 1:num_clusters
            for cci = 1:5
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