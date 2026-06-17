clc; close all; clear all;

% Robustness analysis under ambient-temperature random fluctuation
% Outputs cases:
% 1) d1 original;
% 2-6) d1 with five zero-mean random ambient-temperature fluctuations;
% 7) d2 original;
% 8) d2 extreme-low-temperature case.
%
% The fluctuation cases are defined as:
%   T_out_case(t) = T_out_base(t) + epsilon(t),
%   epsilon(t) ~ N(0, sigma_T^2).
% The generated epsilon trajectory is demeaned to make its sample mean exactly zero.

rng(42, 'twister');
U_pen = 0.008;
N_days = 10;
Z_num = 100;
T = 288;
num_clusters = 8;
gamma = 0.05;
mu_0 = 0;
sigma_0 = 0.01;
Price_samss = [20, 50, 100];
Price_sam_i = 1;
Price_sam = Price_samss(Price_sam_i);
s_sam_number = 801;
s_sam = -4 + 0.01*(0:s_sam_number-1)';

% ===== User-defined robustness cases =====
d1 = 1;                          % selected day for ambient fluctuation
d2 = 4;                          % selected day for extreme low-temperature test
fluctuation_variances = [0.25, 1, 4, 9, 16];  % five variances of zero-mean ambient-temperature noise, unit: degC^2
extreme_low_shift = -10;          % extreme low-temperature case, unit: degC
noise_seed_base = 20;           % fixed seed for reproducible random fluctuation trajectories

% Load common system data and build base thermal model
[eta_, beta_, alpha_, delta_new_base, C, R, P_min, P_max, T_ref_min, T_ref_max, R_1, E, alpha_ij, COP] = build_base_model(T, N_days, Z_num);

% Same partitioning for all robustness cases
[b_result, ~] = Zone_cluster(E, num_clusters, alpha_ij);

% Load prices and quantile forecasts
price_data = load('Data_10days.mat');
Pri_rt_Feb = price_data.Pri_rt_Feb/1000;
quantiles_all = load('Data_quan_Final_10.mat');
quantiles_x = quantiles_all.(['x_' num2str(Price_sam)])(1,2:end-1);

%% Define cases
case_name = {};
case_day = [];
case_variance = [];
case_extreme_shift = [];
case_noise_seed = [];
case_type = {};

case_name{end+1,1} = sprintf('d%d_original', d1);
case_day(end+1,1) = d1; case_variance(end+1,1) = 0; case_extreme_shift(end+1,1) = 0; case_noise_seed(end+1,1) = NaN; case_type{end+1,1} = 'original';

for i = 1:length(fluctuation_variances)
    var_i = fluctuation_variances(i);
    case_name{end+1,1} = sprintf('d%d_zero_mean_noise_var_%gC2', d1, var_i);
    case_day(end+1,1) = d1; case_variance(end+1,1) = var_i; case_extreme_shift(end+1,1) = 0; case_noise_seed(end+1,1) = noise_seed_base + i; case_type{end+1,1} = 'zero_mean_fluctuation';
end

case_name{end+1,1} = sprintf('d%d_original', d2);
case_day(end+1,1) = d2; case_variance(end+1,1) = 0; case_extreme_shift(end+1,1) = 0; case_noise_seed(end+1,1) = NaN; case_type{end+1,1} = 'original';

case_name{end+1,1} = sprintf('d%d_extreme_low_%+gC', d2, extreme_low_shift);
case_day(end+1,1) = d2; case_variance(end+1,1) = 0; case_extreme_shift(end+1,1) = extreme_low_shift; case_noise_seed(end+1,1) = NaN; case_type{end+1,1} = 'extreme_low';

n_cases = length(case_day);
Results = table('Size',[n_cases,15], ...
    'VariableTypes', {'string','double','string','double','double','double','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Case','Day','Type','TempNoiseVariance_C2','TempNoiseStd_C','ExtremeShift_C','NoiseSeed','GT_Ele','GT_Comfort','GT_Total','PR_Ele','PR_Comfort','PR_Total','AbsError_Total','RelError_Total_pct'});

%% Run robustness cases
for ci = 1:n_cases
    opt_day = case_day(ci);
    noise_variance = case_variance(ci);
    noise_std = sqrt(noise_variance);
    extreme_shift_C = case_extreme_shift(ci);
    noise_seed = case_noise_seed(ci);
    fprintf('\n==============================\n');
    fprintf('Case %d/%d: %s\n', ci, n_cases, case_name{ci});
    fprintf('Day = %d, zero-mean noise variance = %.4f degC^2, std = %.4f degC, extreme shift = %+g degC\n', opt_day, noise_variance, noise_std, extreme_shift_C);

    delta_new_case = apply_temperature_random_fluctuation(delta_new_base, R, C, opt_day, noise_std, noise_seed);
    delta_new_case = apply_temperature_constant_shift(delta_new_case, R, C, opt_day, extreme_shift_C);

    % Ground truth: full coupled deterministic optimization
    GT_cost = solve_GT_one_day(opt_day, T, Z_num, U_pen, Pri_rt_Feb, eta_, beta_, alpha_, delta_new_case, P_min, P_max, T_ref_min, T_ref_max);

    % Proposed method: partitioning + deterministic transformation + ASDP
    PR_cost = solve_PR_one_day(opt_day, T, N_days, Z_num, num_clusters, U_pen, Pri_rt_Feb, quantiles_all, quantiles_x, Price_sam, ...
        s_sam, s_sam_number, gamma, mu_0, sigma_0, C, R_1, b_result, eta_, beta_, alpha_, delta_new_case, P_max, T_ref_min, T_ref_max, COP);

    abs_err = abs(PR_cost(3) - GT_cost(3));
    rel_err = 100 * abs_err / max(abs(GT_cost(3)), 1e-9);

    Results.Case(ci) = string(case_name{ci});
    Results.Day(ci) = opt_day;
    Results.Type(ci) = string(case_type{ci});
    Results.TempNoiseVariance_C2(ci) = noise_variance;
    Results.TempNoiseStd_C(ci) = noise_std;
    Results.ExtremeShift_C(ci) = extreme_shift_C;
    Results.NoiseSeed(ci) = noise_seed;
    Results.GT_Ele(ci) = GT_cost(1);
    Results.GT_Comfort(ci) = GT_cost(2);
    Results.GT_Total(ci) = GT_cost(3);
    Results.PR_Ele(ci) = PR_cost(1);
    Results.PR_Comfort(ci) = PR_cost(2);
    Results.PR_Total(ci) = PR_cost(3);
    Results.AbsError_Total(ci) = abs_err;
    Results.RelError_Total_pct(ci) = rel_err;

    fprintf('GT total = %.4f, PR total = %.4f, abs error = %.4f, rel error = %.2f%%\n', GT_cost(3), PR_cost(3), abs_err, rel_err);
end

disp(Results);
%writetable(Results, 'Robustness_Temperature_Results.csv');
%save('Robustness_Temperature_Results.mat', 'Results');
%fprintf('\nSaved results to Robustness_Temperature_Results.csv and Robustness_Temperature_Results.mat\n');

%% Figures
figure('Color','w','Position',[100 100 1050 420]);

Y = [Results.GT_Total, Results.PR_Total];
b = bar(Y, 'grouped', 'LineWidth', 0.8);
b(1).FaceColor = [0.72 0.72 0.72];
b(2).FaceColor = [0.25 0.45 0.75];

n_cases = height(Results);
x = 1:n_cases;

case_labels = {'Day1 original', '\sigma^2=0.25', '\sigma^2=1', '\sigma^2=4', ...
               '\sigma^2=9', '\sigma^2=16', 'Day4 original', 'Extreme low'};

set(gca, 'XTick', x, 'XTickLabel', case_labels, ...
    'FontName','Arial', 'FontSize',11, 'LineWidth',1, 'Box','off');
ylabel('Total cost ($)', 'FontName','Arial', 'FontSize',12);

legend(b, {'GT','PR-ASDP'}, 'Location','northwest', 'Orientation','horizontal', 'Box','off');

hold on;

yl = ylim;
xline(6.5, '-', 'Color','black', 'LineWidth', 1.2, 'HandleVisibility','off');
ylim(yl);

% Annotate absolute difference for each
% group
abs_diff = abs(Results.PR_Total - Results.GT_Total);
y_max = max(Y, [], 2);
yrange = max(Y(:)) - min(Y(:));
if yrange == 0
    yrange = max(Y(:));
end

for i = 1:n_cases
    text(i, y_max(i) + 0.035*yrange, sprintf('|\\Delta|=%.3f', abs_diff(i)), ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontName','Arial', 'FontSize',11);
end

grid on;
set(gca, 'GridAlpha',0.15, 'TickDir','out');
ylim([0, max(y_max + 0.12*yrange)]);

% Optional export
exportgraphics(gcf, 'New_fig_temperature_bar.pdf', 'ContentType', 'vector');

%% =========================== Functions ===========================
function delta_new_case = apply_temperature_random_fluctuation(delta_new_base, R, C, day_i, noise_std, noise_seed)
    delta_new_case = delta_new_base;
    if noise_std <= 0
        return;
    end
    if ~isnan(noise_seed)
        rng(noise_seed, 'twister');
    end
    T = size(delta_new_base, 1);
    Z_num = size(delta_new_base, 3);

    % Common ambient-temperature fluctuation shared by all zones.
    % Demeaning makes the sample mean exactly zero for the selected day.
    eps_T = noise_std * randn(T, 1);
    eps_T = eps_T - mean(eps_T);

    add_term = zeros(T,1,Z_num);
    for z = 1:Z_num
        add_term(:,1,z) = eps_T/(R(z,z)*C(z));
    end
    delta_new_case(:, day_i, :) = delta_new_case(:, day_i, :) + add_term;
end

function delta_new_case = apply_temperature_constant_shift(delta_new_base, R, C, day_i, shift_C)
    delta_new_case = delta_new_base;
    if abs(shift_C) < 1e-12
        return;
    end
    T = size(delta_new_base, 1);
    Z_num = size(delta_new_base, 3);
    add_term = zeros(1,1,Z_num);
    for z = 1:Z_num
        add_term(1,1,z) = shift_C/(R(z,z)*C(z));
    end
    delta_new_case(:, day_i, :) = delta_new_case(:, day_i, :) + repmat(add_term, T, 1, 1);
end

function GT_cost = solve_GT_one_day(day_i, T, Z_num, U_pen, Pri_rt_Feb, eta_, beta_, alpha_, delta_new, P_min, P_max, T_ref_min, T_ref_max)
    v_p = sdpvar(Z_num,T,'full');
    v_tem = sdpvar(Z_num,T+1,'full');
    cost_ele = (1/12)*sum(Pri_rt_Feb((day_i-1)*T+1:day_i*T,1).*sum(v_p,1)');
    cost_comfort = U_pen * sum(v_tem(:,2:end).^2, 'all');
    Objective = cost_ele + cost_comfort;
    Constraints = [v_tem(:,1) == 0, P_min <= v_p, v_p <= P_max, T_ref_min <= v_tem(:,2:end), v_tem(:,2:end) <= T_ref_max];
    Constraints_model = [];
    for t = 1:T
        Constraints_model = [Constraints_model, v_tem(:,t+1) == eta_.*v_tem(:,t) + beta_.*v_p(:,t) + reshape(delta_new(t,day_i,:), [], 1) + alpha_*v_tem(:,t)];
    end
    ops = sdpsettings('solver','gurobi','verbose',0);
    diagnostics = optimize([Constraints, Constraints_model], Objective, ops);
    if diagnostics.problem ~= 0
        warning('GT optimization warning on day %d: %s', day_i, diagnostics.info);
    end
    GT_cost = [value(cost_ele), value(cost_comfort), value(Objective)];
end

function PR_cost = solve_PR_one_day(opt_day, T, N_days, Z_num, num_clusters, U_pen, Pri_rt_Feb, quantiles_all, quantiles_x, Price_sam, ...
        s_sam, s_sam_number, gamma, mu_0, sigma_0, C, R_1, b_result, eta_, beta_, alpha_, delta_new, P_max, T_ref_min, T_ref_max, COP)

    [Eta_, Beta_, Delta_new_, P_c_min, P_c_max, T_c_ref_min, T_c_ref_max] = Zone_aggre(num_clusters,T,N_days,C,R_1,b_result,eta_,beta_,delta_new,P_max,T_ref_min,T_ref_max,COP);
    num_equiv_zones = 5*num_clusters;

    % Chance-constraint deterministic bounds
    T_c_ref_min_CC = zeros(num_equiv_zones, T);
    T_c_ref_max_CC = zeros(num_equiv_zones, T);
    z_score = norminv(1 - gamma, mu_0, sigma_0);
    for t = 1:T
        xi_over_t = mu_0 + t * sigma_0 * z_score;
        xi_under_t = mu_0 - t * sigma_0 * z_score;
        T_c_ref_min_CC(:, t) = T_c_ref_min(:, t) - xi_under_t;
        T_c_ref_max_CC(:, t) = T_c_ref_max(:, t) - xi_over_t;
    end

    % Offline value-function calculation for the selected day only
    u = zeros(s_sam_number, 2);
    u(:,1) = s_sam;
    u(:,2) = U_pen*2*u(:,1);
    v = zeros(T, num_equiv_zones, s_sam_number);
    quantile_y = quantiles_all.(['y_' num2str(opt_day-1) '_' num2str(Price_sam)])/1000;
    for zone = 1:num_equiv_zones
        for t = T:-1:1
            if t == T
                v(t,zone,:) = zeros(1,1,s_sam_number);
            else
                quantiles_rt = quantile_y(t+1,:);
                for ii = 1:s_sam_number
                    v(t,zone,ii) = F_value(s_sam, quantiles_x, quantiles_rt, s_sam(ii), Eta_(zone,1), Beta_(zone,1), u, squeeze(v(t+1,zone,:)), P_c_max(zone,t), Delta_new_(t,opt_day,zone));
                end
            end
        end
        if mod(zone,10) == 0
            fprintf('  Finish value function: zone %d/%d\n', zone, num_equiv_zones);
        end
    end

    % Online scheduling
    DP_v = zeros(T,num_equiv_zones);
    DP_tem = zeros(T,num_equiv_zones);
    DP_tem_recover = zeros(T,Z_num);
    DP_v_recover = zeros(T,Z_num);
    DP_tem_aggre = zeros(T,num_equiv_zones);

    for t = 1:T
        for zone = 1:num_equiv_zones
            if t == 1
                tem_ini = 0;
            else
                tem_ini = DP_tem_aggre(t-1,zone);
            end
            v_dt = squeeze(v(t,zone,:));
            v_dt_0 = diff(v_dt);
            v_dt1 = [v_dt(1); v_dt_0];
            price_rt = Pri_rt_Feb((opt_day-1)*T+t,1);
            s_phys_max = Eta_(zone,1)*tem_ini + Beta_(zone,1)*P_c_max(zone,t) + Delta_new_(t,opt_day,zone);
            s_phys_min = Eta_(zone,1)*tem_ini + Beta_(zone,1)*P_c_min(zone,t) + Delta_new_(t,opt_day,zone);
            s_valid_min = max(s_phys_min, T_c_ref_min_CC(zone,t));
            s_valid_max = min(s_phys_max, T_c_ref_max_CC(zone,t));
            if s_valid_min > s_valid_max
                warning('Infeasible candidate interval at t=%d, zone=%d. Relaxing to physical range.', t, zone);
                s_valid_min = s_phys_min; s_valid_max = s_phys_max;
            end
            [~, min_ind] = min(abs(s_sam - s_valid_min));
            [~, max_ind] = min(abs(s_sam - s_valid_max));
            if min_ind > max_ind
                tmp = min_ind; min_ind = max_ind; max_ind = tmp;
            end
            s = s_sam(min_ind:max_ind);
            p = (s - Eta_(zone,1)*tem_ini - Delta_new_(t,opt_day,zone)) / Beta_(zone,1);
            mask = (p >= P_c_min(zone,t)-1e-8) & (p <= P_c_max(zone,t)+1e-8);
            s = s(mask); p = p(mask);
            if isempty(s)
                p0 = min(max(0, P_c_min(zone,t)), P_c_max(zone,t));
                s0 = Eta_(zone,1)*tem_ini + Beta_(zone,1)*p0 + Delta_new_(t,opt_day,zone);
                s = s0; p = p0;
            end
            c = zeros(length(s), 4);
            for i = 1:length(s)
                c(i,1) = (1/12)*price_rt*p(i);
                c(i,2) = U_pen*(s(i)^2);
                c(i,3) = V_cal(v_dt1, s_sam, s(i));
                c(i,4) = sum(c(i,1:3));
            end
            [~, idx] = min(c(:,4));
            DP_v(t,zone) = p(idx);
            DP_tem(t,zone) = s(idx);
        end

        if t == 1
            DP_tem_ini = zeros(Z_num,1);
        else
            DP_tem_ini = reshape(DP_tem_recover(t-1,:), [], 1);
        end
        [DP_v_recover_t, DP_tem_recover_t] = Zone_recovery(Z_num,num_clusters,b_result,reshape(DP_v(t,:), [], 1),T_ref_max,T_ref_min,P_max,opt_day,eta_,beta_,delta_new,alpha_,DP_tem_ini,t);
        DP_v_recover(t,:) = DP_v_recover_t(:,1);
        DP_tem_recover(t,:) = DP_tem_recover_t(:,2:end);

        C_1 = C(1:20,1);
        for cc = 1:num_clusters
            if num_clusters == 1
                zone_list = 1:20;
            else
                zone_list = find(b_result(:, cc) >= 0.99);
            end
            for cci = 1:5
                cluster_idx = (cci-1)*num_clusters+cc;
                if isscalar(zone_list)
                    DP_tem_aggre(t,cluster_idx) = DP_tem_recover(t,(cci-1)*20+zone_list);
                else
                    zone_C = C_1(zone_list);
                    C_sum = sum(zone_C);
                    DP_tem_aggre(t,cluster_idx) = sum(zone_C.*reshape(DP_tem_recover(t,(cci-1)*20+zone_list), [], 1))/C_sum;
                end
            end
        end
    end

    PR_ele = (1/12) * sum(Pri_rt_Feb((opt_day-1)*T+1:opt_day*T, 1) .* sum(DP_v_recover, 2));
    PR_comfort = U_pen * sum(DP_tem_recover.^2, 'all');
    PR_cost = [PR_ele, PR_comfort, PR_ele + PR_comfort];
end

function [eta_, beta_, alpha_, delta_new, C, R, P_min, P_max, T_ref_min, T_ref_max, R_1, E, alpha_ij, COP] = build_base_model(T, N_days, Z_num)
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
    C = repmat(C0, 5, 1);
    R = 10000*ones(Z_num,Z_num);
    for ii = 1:5
        R((ii-1)*20+1:ii*20,(ii-1)*20+1:ii*20) = R0;
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
end

function [F_val, w_x1, w_x2, num_filtered] = F_value(s_sam, quantiles_x, quantiles_rt, s, eta_, beta_, u, v, Pmax, delta)
    x1 = eta_*s + beta_*Pmax + delta;
    x2 = eta_*s + delta;
    [w_x1, w_x2] = deal(v_u(s_sam, v, u, x1), v_u(s_sam, v, u, x2));
    P1 = eta_*w_x1*F_price(quantiles_x, quantiles_rt, -12*beta_*w_x1);
    P2 = eta_*w_x2*(1 - F_price(quantiles_x, quantiles_rt, -12*beta_*w_x2));
    weights0 = diff(quantiles_x);
    weights = [weights0(1), weights0];
    mask = (quantiles_rt >= -12*beta_*w_x1) & (quantiles_rt <= -12*beta_*w_x2);
    filtered_values = quantiles_rt(mask);
    filtered_weights = weights(mask);
    P3 = (eta_/(12*beta_))*sum(filtered_values .* filtered_weights);
    F_val = P1 + P2 - P3;
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
        p = quantiles_x(end);
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
    diff_vu = vx + ux;
end

function Vx = V_cal(v_dt1, s_sam, x)
    Vx = 0;
    for i = 1:length(s_sam)
        Vx = Vx + v_dt1(i)*max(x - s_sam(i), 0);
    end
end

function [PR_v_p_final_proposed, PR_v_tem_final_proposed] = Zone_recovery(Z_num,num_clusters,b_result,PR_v_c_p_opt,T_ref_max,T_ref_min,P_max,day_i,eta_,beta_,delta_new_in,alpha_,tem_ini,t_now)
    rng(42, 'twister');
    T_local = size(PR_v_c_p_opt,2);
    PR_v_p_final_proposed = zeros(Z_num, T_local);
    PR_v_tem_final_proposed = zeros(Z_num, T_local+1);
    for t = 1:T_local
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
                    if T_local >= 2
                        current_temps = PR_v_tem_final_proposed(zone_list, t);
                        current_temps_ = PR_v_tem_final_proposed(:, t);
                        K = eta_(zone_list).*current_temps + reshape(delta_new_in(t, day_i, zone_list), [], 1) + alpha_(zone_list,:)*current_temps_ + 0.1*randn(n_sub,1);
                        ub_t = P_max(zone_list,t);
                    else
                        current_temps = tem_ini(zone_list, 1);
                        current_temps_ = tem_ini(:, 1);
                        K = eta_(zone_list).*current_temps + reshape(delta_new_in(t_now, day_i, zone_list), [], 1) + alpha_(zone_list,:)*current_temps_ + 0.1*randn(n_sub,1);
                        ub_t = P_max(zone_list,t_now);
                    end
                    B = beta_(zone_list);
                    model.Q = sparse(diag(B.^2));
                    model.obj = 2 * (B .* K);
                    model.A = sparse(ones(1, n_sub));
                    model.rhs = p_cluster_target;
                    model.sense = '=';
                    model.lb = zeros(n_sub, 1);
                    model.ub = ub_t;
                    params.OutputFlag = 0;
                    results = gurobi(model, params);
                    if ~isfield(results,'x')
                        warning('Gurobi recovery failed; using proportional allocation.');
                        PR_v_p_final_proposed(zone_list,t) = min(ub_t, max(0, p_cluster_target/n_sub));
                    else
                        PR_v_p_final_proposed(zone_list, t) = results.x;
                    end
                end
            end
        end
        if T_local >= 2
            PR_v_tem_final_proposed(:, t+1) = eta_.*PR_v_tem_final_proposed(:, t) + beta_.*PR_v_p_final_proposed(:, t) + reshape(delta_new_in(t, day_i, :), [], 1) + alpha_*PR_v_tem_final_proposed(:, t);
        else
            PR_v_tem_final_proposed(:, t+1) = eta_.*tem_ini + beta_.*PR_v_p_final_proposed(:, t) + reshape(delta_new_in(t_now, day_i, :), [], 1) + alpha_*tem_ini;
            PR_v_tem_final_proposed(:, t+1) = max(PR_v_tem_final_proposed(:, t+1), T_ref_min(:, t_now));
            PR_v_tem_final_proposed(:, t+1) = min(PR_v_tem_final_proposed(:, t+1), T_ref_max(:, t_now));
        end
    end
end

function [b_result,e_result] = Zone_cluster(E,num_clusters,alpha_ij)
    num_nodes = 20;
    B = 1:num_nodes;
    num_edges = size(E, 1);
    Clu = 1:num_clusters;
    fprintf('Node number: %d\n', num_nodes);
    fprintf('Edge number: %d\n', num_edges);
    fprintf('Cluster number: %d\n', num_clusters);
    b = binvar(num_nodes, num_clusters, 'full');
    r = binvar(num_edges, num_clusters, 'full');
    e = binvar(num_edges, 1);
    Constraints = [];
    for i = B
        Constraints = [Constraints, sum(b(i, :)) == 1];
    end
    for c = Clu
        Constraints = [Constraints, sum(b(:, c)) >= 1];
    end
    for k = 1:num_edges
        i = E(k, 1); j = E(k, 2);
        Constraints = [Constraints, e(k) == sum(r(k, :))];
        for c = Clu
            Constraints = [Constraints, r(k, c) <= b(i, c), r(k, c) <= b(j, c), b(i, c) + b(j, c) <= r(k, c) + 1];
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
        Constraints = [Constraints, sum(b(:,c)) >= min_cluster_size, sum(b(:,c)) <= max_cluster_size];
    end
    omega = 0.5;
    Objective = sum((1-e).*alpha_ij) - omega*sum(sum(r.*repmat(alpha_ij,1,num_clusters)));
    ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
    sol = optimize(Constraints, Objective, ops);
    if sol.problem ~= 0
        warning('Zone clustering warning: %s', sol.info);
    end
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
            for cci = 1:5
                idx = (cci-1)*num_clusters+cc;
                z = (cci-1)*20+zone_list;
                Eta_(idx,1) = eta_(z,1);
                Beta_(idx,1) = beta_(z,1);
                Delta_new_(:,:,idx) = delta_new(:,:,z);
                P_c_max(idx,:) = P_max(z,:);
                T_c_ref_min(idx,:) = T_ref_min(z,:);
                T_c_ref_max(idx,:) = T_ref_max(z,:);
            end
        else
            for cci = 1:5
                idx = (cci-1)*num_clusters+cc;
                zone_C = C_1(zone_list);
                C_sum = sum(zone_C);
                Beta_(idx,1) = COP/C_sum;
                G_outdoor_sum = 0;
                for i = 1:length(zone_list)
                    z_idx = zone_list(i);
                    G_outdoor_sum = G_outdoor_sum + 1 / R_1(z_idx, z_idx);
                end
                Eta_(idx,1) = 1 - G_outdoor_sum/C_sum;
                weighted_delta = bsxfun(@times, delta_new(:,:,(cci-1)*20+zone_list), reshape(zone_C, 1, 1, []));
                Delta_new_(:,:,idx) = sum(weighted_delta, 3) / C_sum;
                P_c_max(idx,:) = sum(P_max((cci-1)*20+zone_list,:), 1);
                T_c_ref_min(idx,:) = sum(C((cci-1)*20+zone_list,1).*T_ref_min((cci-1)*20+zone_list,:))/sum(C((cci-1)*20+zone_list,1));
                T_c_ref_max(idx,:) = sum(C((cci-1)*20+zone_list,1).*T_ref_max((cci-1)*20+zone_list,:))/sum(C((cci-1)*20+zone_list,1));
            end
        end
    end
end
