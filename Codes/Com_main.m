clc; close all; clear all;

Total = zeros(2,5,3); %5: experiment times; 2: time and obj.;

for type_i = 1:2
    for run_i = 1:5
        if type_i == 1
            Price_sam = 20*run_i;
            c_time = 36;   %time length
        else
            Price_sam = 30;
            c_time = run_i*24;   %time length
        end
    
        rng(20);
        opt_day = 1; %5 = Feb 5
        
        % Parameters
        N_days = 10; %1-6, Feb, 2023
        U_pen = 1.5;
        P_min = 0.0; P_max = 2.0; % MW
        eta_ = 0.93; beta_ = 2.2;
        S0 = 20;
        Sref = [20, repmat(20, 1, 12 * 8), repmat(22, 1, 12 * 4), repmat(21, 1, 12 * 4), repmat(19, 1, 12 * 8)];
        data_in = load('Data.mat');
        delta = data_in.delta;
        Pri_rt_Feb = data_in.Pri_rt_Feb;
        Pri_da_Feb = data_in.Pri_da_Feb;
        
        quantiles_all = load('Data_quan_20.mat');
        s_sam_number = 501;
        quantiles_x = quantiles_all.(['x_' num2str(Price_sam)])(1,2:end-1);
        s_sam = -10 + 0.04*(0:s_sam_number-1)';
        u = zeros(s_sam_number, 2);
        u(:,1) = s_sam;
        u(:,2) = U_pen*sign(u(:,1));
        
        v = zeros(N_days, c_time, s_sam_number);
        v_3 = zeros(N_days, c_time, s_sam_number);
        tic;
        for d = opt_day:opt_day
            for t = c_time:-1:1
                if t == c_time
                    v(d,t,:) = zeros(1, 1, s_sam_number);
                else
                    quantile_y = quantiles_all.(['y_' num2str(d-1) '_' num2str(Price_sam)]);
                    quantiles_rt = quantile_y(t+1,:);
                    
                    for ii = 1:s_sam_number
                        [v(d,t,ii), w_x1, w_x2, v_3(d,t,ii)] = ...
                            F_value(s_sam, quantiles_x, quantiles_rt, s_sam(ii), ...
                                   eta_, beta_, u, squeeze(v(d,t+1,:)), P_max, ...
                                   delta(t,d) - Sref(t+1) + eta_*Sref(t));
                    end
                end
            end
        end
        time_elapsed1 = toc;
        disp(['new time: ' num2str(toc)]);
    
        aa = squeeze(v(d,:,:));
        % figure;
        % plot(aa');
        % xlabel('Time step');
        % ylabel('Value function');
        % figure;
        % plot(1:size(aa,1), max(aa,[],2), 'b-', 1:size(aa,1), min(aa,[],2), 'r-');
        % xlabel('Time step');
        % legend('Max', 'Min');
        
        aa1 = diff(aa, 1, 2);
        aa2 = aa1 > 0;
        disp(['Number of irnormal value: ', num2str(sum(aa2(:)))]);
        
        % Dynamic programming approach
        DP_obj_act = zeros(N_days, 1);
        DP_v = zeros(N_days, c_time);
        DP_tem = zeros(N_days, c_time);
        
        for d = opt_day:opt_day
            tic;
            for t = 1:c_time
                if t == 1
                    tem_ini = 0;
                else
                    tem_ini = DP_tem(d,t-1);
                end
                
                v_dt = squeeze(v(d,t,:));
                v_dt_0 = diff(v_dt);
                v_dt1 = [v_dt(1); v_dt_0];
                
                price_rt = Pri_rt_Feb((d-1)*288+t,1);
                p_2 = P_max;
                s_2 = eta_*tem_ini + beta_*P_max + (delta(t,d) - Sref(t+1) + eta_*Sref(t));
                p_3 = P_min;
                s_3 = eta_*tem_ini + beta_*P_min + (delta(t,d) - Sref(t+1) + eta_*Sref(t));
                
                [~, min_ind] = min(abs(s_sam - s_3));
                [~, max_ind] = min(abs(s_sam - s_2));
                disp([min_ind, max_ind]);
                
                s = s_sam(min_ind:max_ind);
                p = (s - eta_*tem_ini - (delta(t,d) - Sref(t+1) + eta_*Sref(t))) / beta_;
                
                if p(1) < P_min
                    s = s(2:end);
                    p = p(2:end);
                end
                if p(end) > P_max
                    s = s(1:end-1);
                    p = p(1:end-1);
                end
                
                c = zeros(length(s), 1);
                for i = 1:length(s)
                    c(i,1) = (1/12)*price_rt*p(i)+U_pen*abs(s(i))-V_cal(v_dt1, s_sam, s(i));
                end
                [~, idx] = min(c(:,1));
                DP_v(d,t) = p(idx);
                DP_tem(d,t) = s(idx);
            end
            time_elapsed2 = toc;
            
            if DP_v(d,t) > P_max || DP_v(d,t) < P_min
                disp('Wrong!');
            end
            
            DP_obj_act(d,1) = sum((1/12)*Pri_rt_Feb((d-1)*288+(1:c_time),1).*DP_v(d,1:c_time)')+sum(U_pen*abs(DP_tem(d,1:c_time)));
            
            disp(['time: ' num2str(time_elapsed2)]);
            disp(['solution: ' num2str(DP_obj_act(d,1))]);
            %disp('TCL: '); disp(DP_v(d,1:c_time)');
            %disp('Temperature: '); disp(DP_tem(d,1:c_time)');
            disp(min(DP_v(d,:)));
            disp(max(DP_v(d,:)));
            disp(min(DP_tem(d,:)));
            disp(max(DP_tem(d,:)));
            if min(DP_tem(d,:)) <= -10
                disp('Attention!!');
            end
        end
        Total(type_i,run_i,1) = time_elapsed1; Total(type_i,run_i,2) = time_elapsed2; Total(type_i,run_i,3) = DP_obj_act(d,1);
    end
end

Total1 = squeeze(Total(1,:,:));Total2 = squeeze(Total(2,:,:));

% Functions
function [F_val, w_x1, w_x2, num_filtered] = F_value(s_sam, quantiles_x, quantiles_rt, s, eta_, beta_, u, v, Pmax, delta)
    x1 = eta_*s + beta_*Pmax + delta;
    x2 = eta_*s + delta;
    
    [w_x1, w_x2] = deal(v_u(s_sam, v, u, x1), v_u(s_sam, v, u, x2));
    
    P1 = eta_*w_x1*F_price(quantiles_x, quantiles_rt, 12*beta_*w_x1);
    P2 = eta_*w_x2*(1 - F_price(quantiles_x, quantiles_rt, 12*beta_*w_x2));
    
    weights0 = diff(quantiles_x);
    weights = [weights0(1), weights0];
    
    mask = (quantiles_rt >= 12*beta_*w_x1) & (quantiles_rt <= 12*beta_*w_x2);
    filtered_values = quantiles_rt(mask);
    filtered_weights = weights(mask);
    
    P3 = (eta_/(12*beta_))*sum(filtered_values .* filtered_weights);
    F_val = P1 + P2 + P3;
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
    
    diff_vu = vx - ux;
end

function Vx = V_cal(v_dt1, s_sam, x)
    Vx = 0;
    for i = 1:length(s_sam)
        Vx = Vx + v_dt1(i)*max(x - s_sam(i), 0);
    end
end