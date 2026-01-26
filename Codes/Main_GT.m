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

%% Ground Truth (GT) and Decoupled (DE) results
GT_cost = zeros(N_days,3);

for day_i = 1:10
    data_in = load('Data_10days.mat');
    Pri_rt_Feb = data_in.Pri_rt_Feb/1000; %%Unit: $/MWh
    Pri_da_Feb = data_in.Pri_da_Feb/1000;
    
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
    disp(['ele cost_Ground truth: ', num2str(GT_v_cost_ele)]);
    disp(['ele comfort_Ground truth: ', num2str(GT_v_cost_comfort)]);
    disp(['total cost_Ground truth: ', num2str(GT_v_cost_total)]);
    GT_cost(day_i,1) = GT_v_cost_ele;
    GT_cost(day_i,2) = GT_v_cost_comfort;
    GT_cost(day_i,3) = GT_v_cost_total;
    if day_i == 4
        aa = GT_v_tem_opt(:,2:end)';
        save('Tem_GT.mat','aa');
    end
end
