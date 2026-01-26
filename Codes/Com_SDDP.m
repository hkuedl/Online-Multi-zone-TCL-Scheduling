clc; close all; clear all;

Total = zeros(2,5,3); %5: experiment times; 2: time and obj.;

for type_i = 1:2
    for run_i = 1:5
        if type_i == 1
            N_price = 20*run_i;
            H = 36;   %time length
        else
            N_price = 30;
            H = run_i*24;   %time length
        end
        rng(20);
        opt_day = 0; %4 = Feb 5
        
        Utili = "Abs";  %Fast only support Abs
        delta = load('Data.mat').delta;
        Sref = [repmat(20, 1, 12*8), repmat(22, 1, 12*4), repmat(21, 1, 12*4), repmat(19, 1, 12*8)];
        U_pen = 1.5;
        eta_ = 0.93;
        beta_ = 2.2;
        P_min = 0.0;
        P_max = 2.0; %MW
        
        price_rt = load('Data.mat').Pri_rt_Feb((12*24*opt_day+1):(12*24*opt_day+H),1);
        quantiles = load('Data_quan_20.mat').(['y_', num2str(opt_day), '_', num2str(N_price)]);
        % Creating a lattice
        lattice = Lattice.latticeEasy(H, N_price, @(t,i)Price(t,i,quantiles,price_rt(1,1)));
        lattice1 = Lattice.latticeEasy(H, N_price, @(t,i)Price1(t,i,price_rt));
        
        % Visualisation
        %figure;
        %lattice.plotLattice(@(data) num2str(data)) ;
        
        % Run SDDP
        params = sddpSettings('algo.McCount',100,...
                              'stop.iterationMax',20,...
                              'stop.pereiraCoef',2,...
                              'solver','gurobi',...
                              'stop.stdMcCoef',0.2,...
                              'stop.stopWhen','pereira and std');
        
        p = sddpVar(H,1);
        s = sddpVar(H,1);
        s_abs = sddpVar(H,1);
        lattice = compileLattice(lattice,  @(scenario)SDDP_nlds(scenario,p,s,s_abs,delta,Sref,U_pen,eta_,beta_,P_min,P_max,opt_day));
        lattice1 = compileLattice(lattice1,@(scenario)SDDP_nlds(scenario,p,s,s_abs,delta,Sref,U_pen,eta_,beta_,P_min,P_max,opt_day));
        tic;
        output = sddp(lattice,params);
        act_time1 = toc;
        %output1 = sddp(lattice1,params);
        
        % Visualise output
        %plotOutput(output);
        
        lattice_sol = output.lattice;
        %lattice_sol1 = output1.lattice;
        lattice_us_sol = lattice_sol;
        for us_i = 1:H
            if us_i == 1
                us_J = 1;
            else
                us_J = N_price;
            end
            for us_j = 1:us_J
                %lattice_us_sol.graph{us_i}{us_j}.model.cutCoeffs = lattice_sol.graph{us_i}{us_j}.model.cutCoeffs;
                %lattice_us_sol.graph{us_i}{us_j}.model.cutRHS = lattice_sol.graph{us_i}{us_j}.model.cutRHS;
                %lattice_us_sol.graph{us_i}{us_j}.model.modelCntrIdx = lattice_sol.graph{us_i}{us_j}.model.modelCntrIdx;
                lattice_us_sol.graph{us_i}{us_j}.model.c = lattice1.graph{us_i}{us_j}.model.c;
            end
        end
    
        nForward = 1 ;
        objVec_us = zeros(nForward,1);
        objVec_ori = zeros(nForward,1);
        objVec_truth = zeros(nForward,1);
        objVec_us_act = zeros(nForward,1);
        objVec_ori_act = zeros(nForward,1);
        objVec_truth_act = zeros(nForward,1);
        s_0_us = zeros(nForward,H);
        p_0_us = zeros(nForward,H);
        s_abs_us = zeros(nForward,H);
        s_0_ori = zeros(nForward,H);
        p_0_ori = zeros(nForward,H);
        s_0_truth = zeros(nForward,H);
        p_0_truth = zeros(nForward,H);
        
        dataForward = cell(nForward,1);
        for  i = 1:nForward
            tic;
            [objVec_us(i),~,~,solution_us] = forwardPass(lattice_us_sol,lattice_us_sol.randomPath(),params);
            act_time2 = toc;
            [objVec_ori(i),~,~,solution_ori] = forwardPass(lattice_sol,lattice_sol.randomPath(),params);
            %[objVec_truth(i),~,~,solution_truth] = forwardPass(lattice_sol1,lattice_sol1.randomPath(),params); 
            s_0_us(i,:) = lattice_us_sol.getPrimalSolution(s, solution_us);
            p_0_us(i,:) = lattice_us_sol.getPrimalSolution(p, solution_us);
            s_abs_us(i,:) = lattice_us_sol.getPrimalSolution(s_abs, solution_us);
            
            s_0_ori(i,:) = lattice_sol.getPrimalSolution(s, solution_ori);
            p_0_ori(i,:) = lattice_sol.getPrimalSolution(p, solution_ori);
            %s_0_truth(i,:) = lattice_sol1.getPrimalSolution(s, solution_truth);
            %p_0_truth(i,:) = lattice_sol1.getPrimalSolution(p, solution_truth);
            objVec_us_act(i) = (1/12)*sum(price_rt(1:H,1)'.*p_0_us(i,:)) + sum(U_pen*abs(s_0_us(i,:)));
            objVec_ori_act(i) = (1/12)*sum(price_rt(1:H,1)'.*p_0_ori(i,:)) + sum(U_pen*abs(s_0_ori(i,:)));
            %objVec_truth_act(i) = (1/12)*sum(price_rt(1:H,1)'.*p_0_truth(i,:)) + sum(U_pen*abs(s_0_truth(i,:)));
        end
        disp(min(p_0_us(1,:)));
        disp(max(p_0_us(1,:)));
        disp(min(s_0_us(1,:)));
        disp(max(s_0_us(1,:)));
        Total(type_i,run_i,1) = act_time1; Total(type_i,run_i,2) = act_time2; Total(type_i,run_i,3) = objVec_us(i);
    end
end

Total1 = squeeze(Total(1,:,:));Total2 = squeeze(Total(2,:,:));

function [cntr, obj] = SDDP_nlds(scenario, p,s,s_abs, delta,Sref,U_pen,eta_,beta_,P_min,P_max,opt_day)

    t = scenario.getTime() ;
    
    % Fuel cost
    obj = (1/12)*scenario.data*p(t) + U_pen*s_abs(t);
    % power limits
    power = [p(t) >= P_min, p(t) <= P_max];
    
    power_abs = [s_abs(t) - s(t) >= 0, s_abs(t) + s(t) >= 0];
    if t == 1
        model = s(1) - eta_*0 == beta_*p(t) + (delta(1,opt_day+1)-Sref(1)+eta_*20);
    else
        model = s(t) - eta_*s(t-1) == beta_*p(t) + (delta(t,opt_day+1)-Sref(t)+eta_*Sref(t-1));   
    end
    
    cntr = [model, power, power_abs];
end

function out = Price(t,i,quantiles,price_rt)
    if t == 1
        out = price_rt;
    else
        out = quantiles(t,i);
    end
end

function out = Price1(t,i,price_rt)
    out = price_rt(t);
end