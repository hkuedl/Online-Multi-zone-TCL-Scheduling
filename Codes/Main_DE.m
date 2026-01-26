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

%% Two-stage deterministic optimization (First stage)

gamma = 0.05;  %1-gamma = 95%
mu_0 = 0;
sigma_0 = 0.01;

T_ref_min_CC = zeros(Z_num, T);
T_ref_max_CC = zeros(Z_num, T);

z_score = norminv(1 - gamma, mu_0, sigma_0);

for t = 1:T
    xi_over_t = mu_0 + t * sigma_0 * z_score;  % \overline{\xi}_{i,t}
    xi_under_t = mu_0 - t * sigma_0 * z_score; % \underline{\xi}_{i,t}
    
    T_ref_min_CC(:, t) = T_ref_min(:, t) - xi_under_t; 
    T_ref_max_CC(:, t) = T_ref_max(:, t) - xi_over_t;
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
v = zeros(N_days, T, Z_num, s_sam_number);
v_3 = zeros(N_days, T, Z_num, s_sam_number);

% online scheduling
DP_Cost = zeros(N_days,Z_num);
DP_Cost_ele = zeros(N_days,Z_num);
DP_Cost_comfort = zeros(N_days,Z_num);
DP_Cost_sum = zeros(N_days,3);
DP_Time_sum = zeros(N_days,1);

DP_v = zeros(N_days,T,Z_num);
DP_tem = zeros(N_days,T,Z_num);
DP_tem_recover = zeros(N_days,T,Z_num);

for opt_day = 1:10
    disp(['days: ',num2str(opt_day)]);
    for zone = 1:Z_num        
        for t = T:-1:1
            if t == T
                v(opt_day,t,zone,:) = zeros(1, 1, s_sam_number);
            else
                quantile_y = quantiles_all.(['y_' num2str(opt_day-1) '_' num2str(Price_sam)])/1000;
                quantiles_rt = quantile_y(t+1,:);

                for ii = 1:s_sam_number
                    [v(opt_day,t,zone,ii), w_x1, w_x2, v_3(opt_day,t,zone,ii)] = ...
                        F_value(s_sam, quantiles_x, quantiles_rt, s_sam(ii), ...
                               eta_(zone,1), beta_(zone,1), u, squeeze(v(opt_day,t+1,zone,:)), P_max(zone,t), ...
                               delta_new(t,opt_day,zone));
                end
            end
        end
        if mod(zone , 10) == 0
            disp(['Finish ', num2str(zone)]);
        end
    end

    time_elapsed_online = zeros(T,1);
    for t = 1:T
        tic1 = tic;
        for zone = 1:Z_num
            if t == 1
                tem_ini = 0;
            else
                tem_ini = DP_tem_recover(opt_day,t-1,zone);
                %tem_ini = DP_tem(opt_day,t-1,zone);
                %tem_ini = 0.3*DP_tem_recover(opt_day,t-1,zone) + 0.7*DP_tem(opt_day,t-1,zone);
                
            end

            v_dt = squeeze(v(opt_day,t,zone,:));
            v_dt_0 = diff(v_dt);
            v_dt1 = [v_dt(1); v_dt_0];

            price_rt = Pri_rt_Feb((opt_day-1)*288+t,1);
            s_phys_max = eta_(zone,1)*tem_ini + beta_(zone,1)*P_max(zone,t) + delta_new(t,opt_day,zone);
            s_phys_min = eta_(zone,1)*tem_ini + beta_(zone,1)*P_min(zone,t) + delta_new(t,opt_day,zone);
            
            s_constr_min = T_ref_min_CC(zone,t);
            s_constr_max = T_ref_max_CC(zone,t);
            
            s_valid_min = max(s_phys_min, s_constr_min);
            s_valid_max = min(s_phys_max, s_constr_max);

            [~, min_ind] = min(abs(s_sam - s_valid_min));
            [~, max_ind] = min(abs(s_sam - s_valid_max));

            if max_ind-min_ind <= 1 && abs(s_valid_max - s_constr_max) <= 0.01
                min_ind = min_ind-2;
            elseif max_ind-min_ind <= 1 && abs(s_valid_min - s_constr_min) <= 0.01
                max_ind = max_ind+2;
            end
            
            s = s_sam(min_ind:max_ind);
            p = (s - eta_(zone,1)*tem_ini - delta_new(t,opt_day,zone)) / beta_(zone,1);

            if p(1) < P_min(zone,t)
                s = s(2:end);
                p = p(2:end);
            end
            if p(end) > P_max(zone,t)
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
        tic2 = tic;
        if t == 1
            tem_ini_recover = zeros(Z_num,1);
        else
            tem_ini_recover = squeeze(DP_tem_recover(opt_day,t-1,:));
        end
        DP_tem_recover(opt_day,t,:) = eta_ .* tem_ini_recover + ...
                                      beta_ .* squeeze(DP_v(opt_day,t,:)) + ...
                                      reshape(delta_new(t, opt_day, :), [], 1) + ...
                                      alpha_(:, :) * tem_ini_recover;
        DP_tem_recover(opt_day,t,:) = max(squeeze(DP_tem_recover(opt_day,t,:)), T_ref_min(:, t));
        DP_tem_recover(opt_day,t,:) = min(squeeze(DP_tem_recover(opt_day,t,:)), T_ref_max(:, t));
        toc2 = toc(tic2);
        time_elapsed_online(t,1) = toc1/Z_num + toc2;
    end
    
    DP_Cost_ele(opt_day, :) = (1/12) * sum(Pri_rt_Feb((opt_day-1)*T+1:opt_day*T, 1) .* squeeze(DP_v(opt_day,:,:)), 1);
    DP_Cost_comfort(opt_day,:) = U_pen * sum(squeeze(DP_tem_recover(opt_day, :, :)).^2, 1);
    DP_Cost(opt_day,:) = DP_Cost_ele(opt_day,:) + DP_Cost_comfort(opt_day,:);
    disp(['ele cost_real (Proposed): ', num2str(sum(DP_Cost_ele(opt_day,:)))]);
    disp(['comfort cost_real (Proposed): ', num2str(sum(DP_Cost_comfort(opt_day,:)))]);
    disp(['Total cost_real (Proposed): ', num2str(sum(DP_Cost(opt_day,:)))]);
    DP_Cost_sum(opt_day,1) = sum(DP_Cost_ele(opt_day, :));
    DP_Cost_sum(opt_day,2) = sum(DP_Cost_comfort(opt_day, :));
    DP_Cost_sum(opt_day,3) = sum(DP_Cost(opt_day, :));
    DP_Time_sum(opt_day,1) = sum(time_elapsed_online);
    if opt_day == 4
        aa = squeeze(DP_tem_recover(opt_day,:,:));
        save('Tem_DE.mat','aa');
    end
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
