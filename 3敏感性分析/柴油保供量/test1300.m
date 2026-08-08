%% 第二阶段： ADMM 优化
clc; close all;
fprintf('======================================================\n');
fprintf(' 启动 Phase Two: 纳什议价 ADMM 优化 \n');
fprintf('======================================================\n\n');

doPlot = true;
if exist('WS_BATCH_MODE', 'var') && WS_BATCH_MODE
    doPlot = false;
end

%% ==================== 1. 期望场景数据 ====================
T = 24; S = 4;
try
    load('Phase2_Clustered_Scenarios.mat', 'Expected_Ptot', 'Expected_Pg', 'Expected_Pd', 'Prob_Ptot', 'Prob_Pg', 'Prob_Pd');
catch
    error('未检测到 Phase2_Clustered_Scenarios.mat');
end

season_days = [92,92,91,90];
w_season = season_days / sum(season_days);   

P_demand = zeros(S, T); P_green  = zeros(S, T); P_diesel = zeros(S, T);
for s = 1:S
    P_demand(s,:) = max(0, Prob_Ptot(s,:) * squeeze(Expected_Ptot(s,:,:)));
    P_green(s,:)  = max(0, Prob_Pg(s,:)   * squeeze(Expected_Pg(s,:,:)));
    P_diesel(s,:) = max(0, Prob_Pd(s,:)   * squeeze(Expected_Pd(s,:,:)));
end

%% ==================== 2. 参数设定 ====================
price_grid = zeros(1,T); 
for t = 1:T
    h = t-1;
    if (h>=8 && h<12) || (h>=18 && h<22)
        price_grid(t) = 1.10; 
    elseif (h>=0 && h<8) || (h>=22 && h<24)
        price_grid(t) = 0.45; 
    else
        price_grid(t) = 0.75; 
    end
end
price_spot_green  = 0.40;   
price_spot_diesel = 0.90;   
price_shortage    = 2.10;   
price_prod_value  = 1.60;   
alpha_tr = 0.0005; beta_tr  = 0.004;           
I_out_prob = 0.06 * ones(1,T); I_out_prob([11:15, 20:24]) = 0.18; % 设定断电概率
rng(42); % 设定随机种子，保证每次运行生成的 0-1 停电状态一致，方便复现和对比

I_state = zeros(S, T); % 初始化真实的 0-1 断电状态矩阵
for s = 1:S
    I_state(s,:) = rand(1, T) < I_out_prob; 
end

P_contract_diesel =1300; % (kW)

%% ==================== 3. 非合作基准 ====================
P_grid0  = zeros(S, T); P_short0 = zeros(S, T);
T_p_0 = 0; T_g_0 = 0; T_d_0 = 0;
for s = 1:S
    P_grid0(s,:)  = P_demand(s,:) .* (1 - I_state(s,:)); 
    P_short0(s,:) = max(P_demand(s,:) - P_grid0(s,:), 0);
    
    T_p_0 = T_p_0 + w_season(s) * sum(price_prod_value .* P_grid0(s,:) - price_grid .* P_grid0(s,:) - price_shortage .* P_short0(s,:));
    T_g_0 = T_g_0 + w_season(s) * sum(price_spot_green .* P_green(s,:));
    T_d_0 = T_d_0 + w_season(s) * sum(price_spot_diesel .* P_diesel(s,:));
end
fprintf('非合作基准预期日利润 (T_p_0, T_g_0, T_d_0):\n  %.4f,  %.4f,  %.4f\n\n', T_p_0, T_g_0, T_d_0);

%% ==================== 4. P1： ====================
fprintf('--- 运行 P1 ---\n');
maxIter_P1 = 80; tol_primal = 0.1; tol_dual = 0.1; rho_p1 = 0.001; 

lambda_g = zeros(S, T); lambda_d = zeros(S, T);
val_Pg2p_gec = zeros(S, T); val_Pd2p_dgc = zeros(S, T);
Pg_gec_prev = zeros(S, T); Pd_dgc_prev = zeros(S, T);

v_Pg_pm = zeros(S,T); v_Pd_pm = zeros(S,T); v_Pgrid_pm = zeros(S,T); v_Pshort_pm = zeros(S,T);
v_Pspot_g = zeros(S,T); v_Pspot_d = zeros(S,T);
val_delta_p = zeros(1,S); val_delta_g = zeros(1,S); val_delta_d = zeros(1,S);
Obj_PM_store  = nan(1,maxIter_P1); Obj_GEC_store = nan(1,maxIter_P1); Obj_DGC_store = nan(1,maxIter_P1);
r_g_store = nan(1,maxIter_P1); r_d_store = nan(1,maxIter_P1);

for iter = 1:maxIter_P1
    r_g_max = 0; r_d_max = 0; s_g_max = 0; s_d_max = 0;
    
    for s = 1:S
        % PM 
        Pg_pm = sdpvar(1,T); Pd_pm = sdpvar(1,T); Pgrid_pm = sdpvar(1,T); Pshort_pm = sdpvar(1,T);
        C_pm = [ Pgrid_pm + Pg_pm + Pd_pm + Pshort_pm == P_demand(s,:), ...
                 Pgrid_pm >= 0, Pg_pm >= 0, Pd_pm >= 0, Pshort_pm >= 0, ...
                 Pgrid_pm <= P_demand(s,:) .* (1 - I_state(s,:)), ...
                 Pd_pm <= P_contract_diesel .* I_state(s,:) ]; 
                 
        delta_p_expr = sum(price_prod_value.*(Pgrid_pm+Pg_pm+Pd_pm) - price_grid.*Pgrid_pm - price_shortage.*Pshort_pm);
        Obj_pm = delta_p_expr - sum(lambda_g(s,:).*Pg_pm) - sum(lambda_d(s,:).*Pd_pm) - (rho_p1/2)*sum((Pg_pm - val_Pg2p_gec(s,:)).^2) - (rho_p1/2)*sum((Pd_pm - val_Pd2p_dgc(s,:)).^2);
        optimize(C_pm, -Obj_pm, sdpsettings('solver','mosek','verbose',0));
        v_Pg_pm(s,:) = value(Pg_pm); v_Pd_pm(s,:) = value(Pd_pm); v_Pgrid_pm(s,:) = value(Pgrid_pm); v_Pshort_pm(s,:) = value(Pshort_pm); val_delta_p(s) = value(delta_p_expr);

        % GEC 
        Pg_gec = sdpvar(1,T); Pspot_g = sdpvar(1,T);
        C_gec = [Pg_gec + Pspot_g <= P_green(s,:), Pg_gec >= 0, Pspot_g >= 0];
        delta_g_expr = sum(price_spot_green.*Pspot_g) - sum(alpha_tr*(Pg_gec.^2) + beta_tr*Pg_gec);
        Obj_gec = delta_g_expr + sum(lambda_g(s,:).*Pg_gec) - (rho_p1/2)*sum((v_Pg_pm(s,:) - Pg_gec).^2);
        optimize(C_gec, -Obj_gec, sdpsettings('solver','mosek','verbose',0));
        val_Pg2p_gec(s,:) = value(Pg_gec); v_Pspot_g(s,:) = value(Pspot_g); val_delta_g(s) = value(delta_g_expr);

        % DGC 
        Pd_dgc = sdpvar(1,T); Pspot_d = sdpvar(1,T);
        C_dgc = [Pd_dgc + Pspot_d <= P_diesel(s,:), Pd_dgc >= 0, Pspot_d >= 0, ...
                 Pd_dgc <= P_contract_diesel .* I_state(s,:)];
        delta_d_expr = sum(price_spot_diesel.*Pspot_d);
        Obj_dgc = delta_d_expr + sum(lambda_d(s,:).*Pd_dgc) - (rho_p1/2)*sum((v_Pd_pm(s,:) - Pd_dgc).^2);
        optimize(C_dgc, -Obj_dgc, sdpsettings('solver','mosek','verbose',0));
        val_Pd2p_dgc(s,:) = value(Pd_dgc); v_Pspot_d(s,:) = value(Pspot_d); val_delta_d(s) = value(delta_d_expr);

        lambda_g(s,:) = lambda_g(s,:) + rho_p1 * (v_Pg_pm(s,:) - val_Pg2p_gec(s,:));
        lambda_d(s,:) = lambda_d(s,:) + rho_p1 * (v_Pd_pm(s,:) - val_Pd2p_dgc(s,:));
        r_g_max = max(r_g_max, norm(v_Pg_pm(s,:) - val_Pg2p_gec(s,:), 2)); 
        r_d_max = max(r_d_max, norm(v_Pd_pm(s,:) - val_Pd2p_dgc(s,:), 2));
        s_g_max = max(s_g_max, norm(rho_p1*(val_Pg2p_gec(s,:) - Pg_gec_prev(s,:)), 2)); 
        s_d_max = max(s_d_max, norm(rho_p1*(val_Pd2p_dgc(s,:) - Pd_dgc_prev(s,:)), 2));
    end
    
    Obj_PM_store(iter)  = sum(w_season .* val_delta_p); Obj_GEC_store(iter) = sum(w_season .* val_delta_g); Obj_DGC_store(iter) = sum(w_season .* val_delta_d);
    r_g_store(iter) = r_g_max; r_d_store(iter) = r_d_max;
    Pg_gec_prev = val_Pg2p_gec; Pd_dgc_prev = val_Pd2p_dgc;

    if mod(iter,10)==0 || iter==1
        fprintf('  [P1] iter=%2d | r_g_max=%.2e r_d_max=%.2e \n', iter, r_g_max, r_d_max);
    end
    if (r_g_max<=tol_primal && r_d_max<=tol_primal && s_g_max<=tol_dual && s_d_max<=tol_dual)
        fprintf('>>> P1 在第 %d 次收敛！\n', iter);
        Obj_PM_store(iter:end)=Obj_PM_store(iter); Obj_GEC_store(iter:end)=Obj_GEC_store(iter); Obj_DGC_store(iter:end)=Obj_DGC_store(iter);
        r_g_store(iter:end)=r_g_max; r_d_store(iter:end)=r_d_max; break;
    end
end
lastP1 = find(~isnan(Obj_PM_store),1,'last');
delta_p_star = Obj_PM_store(lastP1); delta_g_star = Obj_GEC_store(lastP1); delta_d_star = Obj_DGC_store(lastP1);

%% ==================== 5. P2：纳什资金分配 ADMM ====================
fprintf('\n--- 运行 P2: 纳什议价 ADMM ---\n');
Scale = 800; 
dp_s = delta_p_star/Scale; dg_s = delta_g_star/Scale; dd_s = delta_d_star/Scale;
Tp0_s = T_p_0/Scale;       Tg0_s = T_g_0/Scale;       Td0_s = T_d_0/Scale;

maxIter_P2 = 80; tol_P2 = 0.0001; rho_p2 = 0.5;
gamma_upd = 0; gamma_cg = 0; gamma_cd = 0;
val_Cgp_gec = 0.5*(delta_g_star-T_g_0)/Scale;
val_Cdp_dgc = 0.5*(delta_d_star-T_d_0)/Scale;
val_Updp_dgc = 0.1*max(delta_d_star-T_d_0,0)/Scale;
Updp_hist = nan(1,maxIter_P2); Cgp_hist  = nan(1,maxIter_P2); Cdp_hist  = nan(1,maxIter_P2);

for iter = 1:maxIter_P2
    Updp_pm = sdpvar(1,1); Cgp_pm = sdpvar(1,1); Cdp_pm = sdpvar(1,1);
    surplus_p = dp_s + Updp_pm - Cgp_pm - Cdp_pm - Tp0_s;
    Obj2_pm = log(surplus_p) - gamma_upd*Updp_pm - gamma_cg*Cgp_pm - gamma_cd*Cdp_pm - (rho_p2/2)*(Updp_pm - val_Updp_dgc)^2 - (rho_p2/2)*(Cgp_pm - val_Cgp_gec)^2 - (rho_p2/2)*(Cdp_pm - val_Cdp_dgc)^2;
    optimize([surplus_p >= 1e-3, Updp_pm >= 0, Cgp_pm >= 0, Cdp_pm >= 0], -Obj2_pm, sdpsettings('solver','mosek','verbose',0));
    v_Updp_pm = value(Updp_pm); v_Cgp_pm = value(Cgp_pm); v_Cdp_pm = value(Cdp_pm);

    Cgp_gec = sdpvar(1,1); surplus_g = dg_s + Cgp_gec - Tg0_s;
    Obj2_gec = log(surplus_g) + gamma_cg*Cgp_gec - (rho_p2/2)*(v_Cgp_pm - Cgp_gec)^2;
    optimize([surplus_g >= 1e-3, Cgp_gec >= 0], -Obj2_gec, sdpsettings('solver','mosek','verbose',0)); val_Cgp_gec = value(Cgp_gec);

    Updp_dgc = sdpvar(1,1); Cdp_dgc = sdpvar(1,1); surplus_d = dd_s - Updp_dgc + Cdp_dgc - Td0_s;
    Obj2_dgc = log(surplus_d) + gamma_upd*Updp_dgc + gamma_cd*Cdp_dgc - (rho_p2/2)*(v_Updp_pm - Updp_dgc)^2 - (rho_p2/2)*(v_Cdp_pm - Cdp_dgc)^2;
    optimize([surplus_d >= 1e-3, Updp_dgc >= 0, Cdp_dgc >= 0], -Obj2_dgc, sdpsettings('solver','mosek','verbose',0)); val_Updp_dgc = value(Updp_dgc); val_Cdp_dgc = value(Cdp_dgc);

    gamma_upd = gamma_upd + rho_p2*(v_Updp_pm - val_Updp_dgc); gamma_cg  = gamma_cg  + rho_p2*(v_Cgp_pm - val_Cgp_gec); gamma_cd  = gamma_cd  + rho_p2*(v_Cdp_pm - val_Cdp_dgc);
    
    Updp_hist(iter) = v_Updp_pm * Scale / 10000; Cgp_hist(iter)  = v_Cgp_pm * Scale / 10000; Cdp_hist(iter)  = v_Cdp_pm * Scale / 10000;
    
    if max([abs(v_Updp_pm-val_Updp_dgc), abs(v_Cgp_pm-val_Cgp_gec), abs(v_Cdp_pm-val_Cdp_dgc)]) <= tol_P2
        fprintf('>>> P2 在第 %d 次收敛\n', iter);
        Updp_hist(iter:end)=Updp_hist(iter); Cgp_hist(iter:end)=Cgp_hist(iter); Cdp_hist(iter:end)=Cdp_hist(iter); break;
    end
end
Updp = v_Updp_pm * Scale; Cgp  = v_Cgp_pm * Scale; Cdp  = v_Cdp_pm * Scale;

%% ==================== 6. 终端对账单 ====================
days_year = 365;

S_machines = 8; 
theta_s = [120000, 85000, 150000, 60000, 90000, 110000, 75000, 130000]; 
Y_s     = [10, 12, 10, 15, 10, 10, 12, 15]; 
r       = 0.06; 
C_cap_p = 0;
for i = 1:S_machines
    CRF = (r * (1 + r)^Y_s(i)) / ((1 + r)^Y_s(i) - 1);
    C_cap_p = C_cap_p + CRF * theta_s(i);
end

% 计算年度能量统计 (kWh/year)
E_PM_grid_0  = sum(season_days .* sum(P_grid0, 2)');
E_PM_short_0 = sum(season_days .* sum(P_short0, 2)');
E_GEC_spot_0 = sum(season_days .* sum(P_green, 2)');
E_DGC_spot_0 = sum(season_days .* sum(P_diesel, 2)');

E_PM_grid_c  = sum(season_days .* sum(v_Pgrid_pm, 2)');
E_PM_gec_c   = sum(season_days .* sum(val_Pg2p_gec, 2)');
E_PM_dgc_c   = sum(season_days .* sum(val_Pd2p_dgc, 2)');
E_PM_short_c = sum(season_days .* sum(v_Pshort_pm, 2)');
E_GEC_spot_c = sum(season_days .* sum(v_Pspot_g, 2)');
E_DGC_spot_c = sum(season_days .* sum(v_Pspot_d, 2)');

% --- 6.2 计算各主体非合作总收益 (Baseline) ---
% PM: 产值 - 主网电费 - 缺电惩罚 - 资本成本
annual_grid_cost_non = 0;
for s = 1:4, annual_grid_cost_non = annual_grid_cost_non + season_days(s) * sum(price_grid .* P_grid0(s,:)); end
T_p_0_ann = price_prod_value * E_PM_grid_0 - annual_grid_cost_non - price_shortage * E_PM_short_0 - C_cap_p;

% GEC: 现货收益
T_g_0_ann = price_spot_green * E_GEC_spot_0;

% DGC: 现货收益
T_d_0_ann = price_spot_diesel * E_DGC_spot_0;

% --- 6.3 计算各主体合作总收益 (Cooperative) ---
% PM: 产值 + 柴油公司赔付 - 资本成本 - 主网费 - 支付给GEC - 支付给DGC - 缺电惩罚
annual_grid_cost_co = 0;
for s = 1:4, annual_grid_cost_co = annual_grid_cost_co + season_days(s) * sum(price_grid .* v_Pgrid_pm(s,:)); end
T_p_coop_ann = price_prod_value * (E_PM_grid_c + E_PM_gec_c + E_PM_dgc_c) + (Updp * days_year) ...
               - C_cap_p - annual_grid_cost_co - (Cgp * days_year) - (Cdp * days_year) - (price_shortage * E_PM_short_c);

% GEC: 现货收益 + PM支付 - 过网费
GEC_wheel_co = 0;
for s = 1:4
    daily_wheel = sum(alpha_tr * (val_Pg2p_gec(s,:).^2) + beta_tr * val_Pg2p_gec(s,:));
    GEC_wheel_co = GEC_wheel_co + daily_wheel * season_days(s);
end
T_g_coop_ann = price_spot_green * E_GEC_spot_c + (Cgp * days_year) - GEC_wheel_co;

% DGC: 现货收益 + PM支付 - 停电赔付
T_d_coop_ann = price_spot_diesel * E_DGC_spot_c + (Cdp * days_year) - (Updp * days_year);

% --- 6.4 打印对比表格 ---
T_total_0_ann = T_p_0_ann + T_g_0_ann + T_d_0_ann;
T_total_coop_ann = T_p_coop_ann + T_g_coop_ann + T_d_coop_ann;

fprintf('\n================ 全年能源交易量对比 (kWh/年) ================\n');
fprintf('[非合作模式]\n  PM: 购主网 = %10.1f | 停产缺口 = %10.1f\n  GEC: 售现货 = %10.1f\n  DGC: 售现货 = %10.1f\n', E_PM_grid_0, E_PM_short_0, E_GEC_spot_0, E_DGC_spot_0);
fprintf('[合作模式]\n  PM: 购主网 = %10.1f | 购绿电 = %10.1f | 购柴油 = %10.1f | 缺口 = %10.1f\n', E_PM_grid_c, E_PM_gec_c, E_PM_dgc_c, E_PM_short_c);

fprintf('\n================ 全年经济收益对比 ($/年) [明细对齐版] ================\n');
fprintf('                非合作利润       合作利润       纳什红利(净提升)     提升比例(%%)\n');
fprintf('  PM(制药厂): %12.2f  %12.2f  %12.2f      %8.2f%%\n', T_p_0_ann, T_p_coop_ann, T_p_coop_ann - T_p_0_ann, (T_p_coop_ann - T_p_0_ann)/abs(T_p_0_ann)*100);
fprintf('  GEC(绿电):  %12.2f  %12.2f  %12.2f      %8.2f%%\n', T_g_0_ann, T_g_coop_ann, T_g_coop_ann - T_g_0_ann, (T_g_coop_ann - T_g_0_ann)/abs(T_g_0_ann)*100);
fprintf('  DGC(柴油):  %12.2f  %12.2f  %12.2f      %8.2f%%\n', T_d_0_ann, T_d_coop_ann, T_d_coop_ann - T_d_0_ann, (T_d_coop_ann - T_d_0_ann)/abs(T_d_0_ann)*100);
fprintf('  ----------------------------------------------------------------------------------\n');
fprintf('  联盟总收益: %12.2f  %12.2f  %12.2f      %8.2f%%\n', T_total_0_ann, T_total_coop_ann, T_total_coop_ann - T_total_0_ann, (T_total_coop_ann - T_total_0_ann)/abs(T_total_0_ann)*100);
%% ==================== 7. 全套可视化 ====================
if doPlot
    fontN = 'Times New Roman'; 
    seasons_name = {'Spring', 'Summer', 'Autumn', 'Winter'};
    figure('Position',[50 50 900 750],'Color','w','Name','Convergence History'); 

    % --- (a) PM 目标函数收敛 ---
    subplot(2,2,3); 
    plot(Obj_PM_store, 'b.-', 'MarkerSize',8, 'LineWidth',1.2); 
    set(gca,'FontName',fontN,'FontSize',11,'XTick',1:20:maxIter_P1,'XLim',[1 maxIter_P1]); 
    ylabel('\it \delta^{p} \rm ($)'); 
    xlabel({'Iteration', '(a) PM objective function convergence'}, 'FontName', fontN, 'FontWeight', 'normal');
    box on; grid on;

    % --- (b) GEC 目标函数收敛 ---
    subplot(2,2,1); 
    plot(Obj_GEC_store, 'b.-', 'MarkerSize',8, 'LineWidth',1.2); 
    set(gca,'FontName',fontN,'FontSize',11,'XTick',1:20:maxIter_P1,'XLim',[1 maxIter_P1]); 
    ylabel('\it \delta^{g} \rm ($)'); 
    xlabel({'Iteration', '(b) GEC objective function convergence'}, 'FontName', fontN, 'FontWeight', 'normal'); 
    box on; grid on;
    
    % --- (c) DGC 目标函数收敛 ---
    subplot(2,2,2); 
    plot(Obj_DGC_store, 'b.-', 'MarkerSize',8, 'LineWidth',1.2); 
    set(gca,'FontName',fontN,'FontSize',11,'XTick',1:20:maxIter_P1,'XLim',[1 maxIter_P1]); 
    ylabel('\it \delta^{d} \rm ($)'); 
    xlabel({'Iteration', '(c) DGC objective function convergence'}, 'FontName', fontN, 'FontWeight', 'normal');
    box on; grid on;
    
    % --- (d) 交易支付项收敛 ---
    subplot(2,2,4); 
    plot(1:maxIter_P2, Cgp_hist, 'b-o', 1:maxIter_P2, Cdp_hist, 'm-d', 1:maxIter_P2, Updp_hist, 'k-*', 'MarkerSize',5, 'LineWidth',1.2); 
    set(gca,'FontName',fontN,'FontSize',11,'XLim',[1 maxIter_P2]); 
    ylabel('Payments 10^4 ($)'); 
    xlabel({'Iteration', '(d) Convergence of bargaining payments'}, 'FontName', fontN, 'FontWeight', 'normal');
    legend('\it C^p_{g}','\it C^p_{d}','\it U^p_{pd}','Location','best'); 
    box on; grid on;

    sub_labels = {'(a)', '(b)', '(c)', '(d)'}; 

    % ================= 1. PM 非合作状态能量获取 =================
    figure('Position',[100 100 1000 650],'Color','w','Name','PM Electricity Sourcing (Non-Coop)');
    for s = 1:4
        subplot(2,2,s); 
        bar(1:24, [P_grid0(s,:)', P_short0(s,:)'], 'stacked'); hold on; 
        plot(1:24, P_demand(s,:), 'k-o', 'LineWidth', 1.5, 'MarkerSize', 4); 
        xlim([0.5 24.5]); box on; set(gca,'FontName',fontN);
        
        ylabel('Electricity (kWh)', 'FontName', fontN);
        xlabel({'Time (h)', [sub_labels{s}, ' PM Non-Coop - ', seasons_name{s}]}, 'FontName', fontN, 'FontWeight', 'normal');
        
        if s == 1, legend('From Main Grid', 'Unmet Demand', 'Total Demand', 'Location','best'); end
    end

    % ================= 2. PM 合作状态能量获取 =================
    figure('Position',[150 150 1000 650],'Color','w','Name','PM Electricity Sourcing (Cooperative)');
    for s = 1:4
        subplot(2,2,s); 
        bar(1:24, [v_Pgrid_pm(s,:)', v_Pg_pm(s,:)', v_Pd_pm(s,:)', v_Pshort_pm(s,:)'], 'stacked'); hold on; 
        plot(1:24, P_demand(s,:), 'k-o', 'LineWidth', 1.5, 'MarkerSize', 4); 
        xlim([0.5 24.5]); box on; set(gca,'FontName',fontN);
        
        ylabel('Electricity (kWh)', 'FontName', fontN);
        xlabel({'Time (h)', [sub_labels{s}, ' PM Coop - ', seasons_name{s}]}, 'FontName', fontN, 'FontWeight', 'normal');
        
        if s == 1, legend('From Main Grid', 'From GEC', 'From DGC', 'Shortage', 'Total Demand', 'Location','best'); end
    end
    
    % ================= 3. GEC 合作状态绿电调度 =================
    figure('Position',[200 200 1000 650],'Color','w','Name','GEC Dispatch (Cooperative)');
    for s = 1:4
        subplot(2,2,s); 
        bar(1:24, [val_Pg2p_gec(s,:)', v_Pspot_g(s,:)'], 'stacked'); hold on; 
        plot(1:24, P_green(s,:), 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4); 
        xlim([0.5 24.5]); box on; set(gca,'FontName',fontN);
        
        ylabel('Electricity (kWh)', 'FontName', fontN);
        xlabel({'Time (h)', [sub_labels{s}, ' GEC Coop - ', seasons_name{s}]}, 'FontName', fontN, 'FontWeight', 'normal');
        
        if s == 1, legend('Sold to PM', 'Sold to Spot Market', 'Total Green Output', 'Location','best'); end
    end
    
    % ================= 4. DGC 合作状态柴油调度 =================
    figure('Position',[250 250 1000 650],'Color','w','Name','DGC Dispatch (Cooperative)');
    for s = 1:4
        subplot(2,2,s); 
        bar(1:24, [val_Pd2p_dgc(s,:)', v_Pspot_d(s,:)'], 'stacked'); hold on; 
        plot(1:24, P_diesel(s,:), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4); 
        xlim([0.5 24.5]); box on; set(gca,'FontName',fontN);
        
        ylabel('Electricity (kWh)', 'FontName', fontN);
        xlabel({'Time (h)', [sub_labels{s}, ' DGC Coop - ', seasons_name{s}]}, 'FontName', fontN, 'FontWeight', 'normal');
        
        if s == 1, legend('Guaranteed to PM', 'Sold to Spot Market', 'Total Diesel Capacity', 'Location','best'); end
    end

    figure('Position',[300 300 650 450],'Color','w','Name','Annual Profit Comparison');
    bar_data = [T_p_0_ann, T_p_coop_ann; T_g_0_ann, T_g_coop_ann; T_d_0_ann, T_d_coop_ann; (T_p_0_ann+T_g_0_ann+T_d_0_ann), (T_p_coop_ann+T_g_coop_ann+T_d_coop_ann)] / 10000; 
    b = bar(bar_data, 'FaceColor','flat', 'EdgeColor', 'k', 'LineWidth', 1);
    b(1).CData = [0.7 0.7 0.7]; b(2).CData = [0.2 0.6 0.8]; 
    set(gca, 'XTickLabel', {'PM', 'GEC', 'DGC', 'Total Coalition'}, 'FontName', fontN, 'FontSize', 12); ylabel('Annual Profit (10k $)', 'FontName', fontN, 'FontSize', 12);
    title('Comparison of Annual Profit', 'FontName', fontN, 'FontSize', 14); legend('Non-coop', 'Cooperative', 'Location', 'northwest', 'FontName', fontN); box on; grid on;

    %% ==================== 8.5 新增：对标顶级期刊的增强可视化 ====================
    
    % --- 新图 1: 对数坐标系下的残差收敛 ---
    figure('Position',[350 100 500 400],'Color','w','Name','ADMM Residual Decay (Log Scale)');
    semilogy(1:maxIter_P1, r_g_store, 'g-', 'LineWidth', 2, 'DisplayName', 'Residual: PM-GEC'); hold on; grid on; box on;
    semilogy(1:maxIter_P1, r_d_store, 'r-', 'LineWidth', 2, 'DisplayName', 'Residual: PM-DGC');
    xlabel('Iteration', 'FontName', fontN, 'FontSize', 11); ylabel('Primal Residuals (Log Scale)', 'FontName', fontN, 'FontSize', 11);
    title('ADMM Residual Decay', 'FontName', fontN, 'FontSize', 12);
    legend('Location','northeast', 'FontName', fontN); set(gca, 'FontName', fontN);

    % --- 新图 2: 四季协作能源构成堆叠图 ---
    figure('Position',[400 150 500 400],'Color','w','Name','Seasonal Energy Composition');
    energy_gec_seasonal = sum(val_Pg2p_gec, 2) .* season_days'; 
    energy_dgc_seasonal = sum(val_Pd2p_dgc, 2) .* season_days'; 
    b_comp = bar(1:4, [energy_gec_seasonal, energy_dgc_seasonal]/1000, 'stacked');
    b_comp(1).FaceColor = [0.9290 0.6940 0.1250]; % 太阳能(黄色)
    b_comp(2).FaceColor = [0.4 0.4 0.4];          % 柴油(灰色)
    set(gca, 'XTick', 1:4, 'XTickLabel', seasons_name, 'FontName', fontN, 'FontSize', 11);
    ylabel('Total Traded Energy (MWh)', 'FontName', fontN, 'FontSize', 11);
    title('Energy Supply by Source', 'FontName', fontN, 'FontSize', 12);
    legend('Solar Energy', 'Diesel Backup', 'Location', 'northwest', 'FontName', fontN);
    grid on; box on;

       % --- 新图 3: 时空全景热力图  ---
    figure('Position',[450 200 800 650],'Color','w','Name','Trading Heatmaps');
    
    % --- 子图 (a) 绿电交易热力图 ---
    subplot(2,1,1); 
    imagesc(1:24, 1:4, val_Pg2p_gec); colorbar;
    set(gca, 'YTick', 1:4, 'YTickLabel', seasons_name, 'XTick', 1:2:24, 'FontName', fontN, 'FontSize', 11);
    xlabel({'Time Slot (Hour)', '(a) Green energy traded to PM (kWh)'}, ...
        'FontName', fontN, 'FontSize', 11, 'FontWeight', 'normal');
    
    % --- 子图 (b) 柴油交易热力图 ---
    subplot(2,1,2); 
    imagesc(1:24, 1:4, val_Pd2p_dgc); colorbar;
    set(gca, 'YTick', 1:4, 'YTickLabel', seasons_name, 'XTick', 1:2:24, 'FontName', fontN, 'FontSize', 11);
    xlabel({'Time Slot (Hour)', '(b) Diesel energy traded to PM (kWh)'}, ...
        'FontName', fontN, 'FontSize', 11, 'FontWeight', 'normal');

    % --- 新图 4: 典型日(夏季)小时级精细化调度曲线 ---
    figure('Position',[500 250 800 400],'Color','w','Name','Hourly Dispatch Profile (Summer)');
    s_target = 2; 
    plot(1:24, v_Pgrid_pm(s_target,:), 'b-s', 'LineWidth', 1.5, 'MarkerFaceColor','w', 'DisplayName', 'Grid Supply'); hold on; grid on; box on;
    plot(1:24, val_Pg2p_gec(s_target,:), 'g-^', 'LineWidth', 1.5, 'MarkerFaceColor','w', 'DisplayName', 'Green Energy');
    plot(1:24, val_Pd2p_dgc(s_target,:), 'r-d', 'LineWidth', 1.5, 'MarkerFaceColor','w', 'DisplayName', 'Diesel Backup');
    plot(1:24, P_demand(s_target,:), 'k--o', 'LineWidth', 2, 'MarkerFaceColor','w', 'DisplayName', 'Total Demand');
    yline(P_contract_diesel, 'r:', 'LineWidth', 2, 'DisplayName', '\it {P_{d2p}^p}'); % 凸显契约上限
    title('Hourly Dispatch Profile during Summer Peak', 'FontName', fontN, 'FontSize', 13);
    xlabel('Time Slot (Hour)', 'FontName', fontN, 'FontSize', 11);
    ylabel('Power (kW)', 'FontName', fontN, 'FontSize', 11);
    legend('Location','best', 'FontName', fontN, 'FontSize', 11);
end

yalmip('clear');
%% ==================== 9. 明细分解输出 ====================
E_PM_grid_0   = sum(season_days .* sum(P_grid0, 2)');
E_PM_short_0  = sum(season_days .* sum(P_short0, 2)');
E_GEC_spot_0  = sum(season_days .* sum(P_green, 2)');
E_DGC_spot_0  = sum(season_days .* sum(P_diesel, 2)');

E_PM_grid_c   = sum(season_days .* sum(v_Pgrid_pm, 2)');
E_PM_gec_c    = sum(season_days .* sum(val_Pg2p_gec, 2)');
E_PM_dgc_c    = sum(season_days .* sum(val_Pd2p_dgc, 2)');
E_PM_short_c  = sum(season_days .* sum(v_Pshort_pm, 2)');
E_GEC_spot_c  = sum(season_days .* sum(v_Pspot_g, 2)');
E_DGC_spot_c  = sum(season_days .* sum(v_Pspot_d, 2)');

% ---------- PM 分解 ----------
PM_rev_non = price_prod_value * E_PM_grid_0;
PM_rev_co  = price_prod_value * (E_PM_grid_c + E_PM_gec_c + E_PM_dgc_c);

PM_penalty_non = 0;
PM_penalty_co  = Updp * days_year;  

S_machines = 8; 
theta_s = [120000, 85000, 150000, 60000, 90000, 110000, 75000, 130000]; 
Y_s     = [10, 12, 10, 15, 10, 10, 12, 15]; 
r       = 0.06; 
C_c_p = 0;
for i = 1:S_machines
    CRF = (r * (1 + r)^Y_s(i)) / ((1 + r)^Y_s(i) - 1);
    C_c_p = C_c_p + CRF * theta_s(i);
end
PM_cap_non = C_c_p;
PM_cap_co  = C_c_p;

annual_grid_cost_non = 0;
for s = 1:4
    annual_grid_cost_non = annual_grid_cost_non + season_days(s) * sum(price_grid .* P_grid0(s,:));
end
PM_cost_grid_non = annual_grid_cost_non;

annual_grid_cost_co = 0;
for s = 1:4
    annual_grid_cost_co = annual_grid_cost_co + season_days(s) * sum(price_grid .* v_Pgrid_pm(s,:));
end
PM_cost_grid_co = annual_grid_cost_co;

PM_cost_gec_non  = 0;
PM_cost_gec_co   = Cgp * days_year;   

PM_cost_dgc_non  = 0;
PM_cost_dgc_co   = Cdp * days_year;   

PM_short_cost_non = price_shortage * E_PM_short_0;
PM_short_cost_co  = price_shortage * E_PM_short_c;

PM_total_non = PM_rev_non + PM_penalty_non - PM_cap_non - PM_cost_grid_non - PM_cost_gec_non - PM_cost_dgc_non - PM_short_cost_non;
PM_total_co  = PM_rev_co  + PM_penalty_co  - PM_cap_co  - PM_cost_grid_co  - PM_cost_gec_co  - PM_cost_dgc_co  - PM_short_cost_co;

% ---------- GEC 分解 ----------
GEC_income_pm_non   = 0;
GEC_income_pm_co    = Cgp * days_year;                   
GEC_spot_rev_non    = price_spot_green * E_GEC_spot_0;
GEC_spot_rev_co     = price_spot_green * E_GEC_spot_c;

GEC_wheel_non = 0;
Annual_Wheel_Cost = 0;
for s = 1:4
    daily_wheel = sum(alpha_tr * (val_Pg2p_gec(s,:).^2) + beta_tr * val_Pg2p_gec(s,:));
    Annual_Wheel_Cost = Annual_Wheel_Cost + daily_wheel * season_days(s);
end
GEC_wheel_co = Annual_Wheel_Cost;

GEC_total_non = GEC_income_pm_non + GEC_spot_rev_non - GEC_wheel_non;
GEC_total_co  = GEC_income_pm_co  + GEC_spot_rev_co  - GEC_wheel_co;

% ---------- DGC 分解 ----------
DGC_income_pm_non   = 0;
DGC_income_pm_co    = Cdp * days_year;   
DGC_spot_rev_non    = price_spot_diesel * E_DGC_spot_0;
DGC_spot_rev_co     = price_spot_diesel * E_DGC_spot_c;

DGC_penalty_non     = 0;
DGC_penalty_co      = Updp * days_year;  

DGC_total_non = DGC_income_pm_non + DGC_spot_rev_non - DGC_penalty_non;
DGC_total_co  = DGC_income_pm_co  + DGC_spot_rev_co  - DGC_penalty_co;

%% ==================== 10. 控制台表格打印 ====================
fprintf('\n================= 制药厂（PM）收益分解 ($/year) =================\n');
fprintf('%-50s %-15s %-15s\n','Item','Cooperative','Noncooperative');
fprintf('%-50s %-15.2f %-15.2f\n','Annual revenue', PM_rev_co, PM_rev_non);
fprintf('%-50s %-15.2f %-15.2f\n','penalties from the diesel-guarantee company', PM_penalty_co, PM_penalty_non);
fprintf('%-50s %-15.2f %-15.2f\n','capital cost', PM_cap_co, PM_cap_non);
fprintf('%-50s %-15.2f %-15.2f\n','Electricity costs with power main-grid', PM_cost_grid_co, PM_cost_grid_non);
fprintf('%-50s %-15.2f %-15.2f\n','Electricity costs with green-electricity company', PM_cost_gec_co, PM_cost_gec_non);
fprintf('%-50s %-15.2f %-15.2f\n','Electricity costs with diesel-guarantee company', PM_cost_dgc_co, PM_cost_dgc_non);
fprintf('%-50s %-15.2f %-15.2f\n','Shortage penalty cost', PM_short_cost_co, PM_short_cost_non);
fprintf('%-50s %-15.2f %-15.2f\n','Total income', PM_total_co, PM_total_non);

fprintf('\n================= 绿电公司（GEC）收益分解 ($/year) =================\n');
fprintf('%-50s %-15s %-15s\n','Item','Cooperative','Noncooperative');
fprintf('%-50s %-15.2f %-15.2f\n','Electricity income with pharmaceutical manufacturer', GEC_income_pm_co, GEC_income_pm_non);
fprintf('%-50s %-15.2f %-15.2f\n','Annual revenue from spot market', GEC_spot_rev_co, GEC_spot_rev_non);
fprintf('%-50s %-15.2f %-15.2f\n','wheeling charge', GEC_wheel_co, GEC_wheel_non);
fprintf('%-50s %-15.2f %-15.2f\n','Total income', GEC_total_co, GEC_total_non);

fprintf('\n================= 柴油公司（DGC）收益分解 ($/year) =================\n');
fprintf('%-50s %-15s %-15s\n','Item','Cooperative','Noncooperative');
fprintf('%-50s %-15.2f %-15.2f\n','Electricity income with pharmaceutical manufacturer', DGC_income_pm_co, DGC_income_pm_non);
fprintf('%-50s %-15.2f %-15.2f\n','Annual revenue from spot market', DGC_spot_rev_co, DGC_spot_rev_non);
fprintf('%-50s %-15.2f %-15.2f\n','penalty to the pharmaceutical manufacturer', DGC_penalty_co, DGC_penalty_non);
fprintf('%-50s %-15.2f %-15.2f\n','Total income', DGC_total_co, DGC_total_non);

fprintf('\n================= 交易路径/主体行为 (kWh/year) =================\n');
fprintf('%-35s %-20s %-20s\n','Path','Baseline','Cooperative');
fprintf('%-35s %-20.1f %-20.1f\n','PM 从主网购电', E_PM_grid_0, E_PM_grid_c);
fprintf('%-35s %-20.1f %-20.1f\n','PM 从 GEC 购电', 0, E_PM_gec_c);
fprintf('%-35s %-20.1f %-20.1f\n','PM 从 DGC 购电', 0, E_PM_dgc_c);
fprintf('%-35s %-20.1f %-20.1f\n','PM 未满足缺电量', E_PM_short_0, E_PM_short_c);
fprintf('%-35s %-20.1f %-20.1f\n','GEC 售往现货市场', E_GEC_spot_0, E_GEC_spot_c);
fprintf('%-35s %-20.1f %-20.1f\n','DGC 售往现货市场', E_DGC_spot_0, E_DGC_spot_c);
%% =========================================================================
% [新增] 将所有运行结果导出为 Excel 和 MAT 文件
% =========================================================================
filename_excel = 'P_contract_diesel1300.xlsx';
filename_mat   = 'P_contract_diesel1300.mat';

% 1. 全年经济收益对比 (Profit Summary)
Profit_Summary = table(...
    {'PM (Pharm)'; 'GEC (Green)'; 'DGC (Diesel)'; 'Coalition Total'}, ...
    [T_p_0_ann; T_g_0_ann; T_d_0_ann; T_p_0_ann+T_g_0_ann+T_d_0_ann], ...
    [T_p_coop_ann; T_g_coop_ann; T_d_coop_ann; T_p_coop_ann+T_g_coop_ann+T_d_coop_ann], ...
    [T_p_coop_ann - T_p_0_ann; T_g_coop_ann - T_g_0_ann; T_d_coop_ann - T_d_0_ann; (T_p_coop_ann+T_g_coop_ann+T_d_coop_ann) - (T_p_0_ann+T_g_0_ann+T_d_0_ann)], ...
    'VariableNames', {'Agent_Name', 'Non_Coop_Profit', 'Coop_Profit', 'Nash_Surplus'});
writetable(Profit_Summary, filename_excel, 'Sheet', 'Profit_Summary');

% 2. 交易路径/主体行为 (Energy Volume)
Energy_Volume = table(...
    {'PM 从主网购电'; 'PM 从 GEC 购电'; 'PM 从 DGC 购电'; 'PM 未满足缺电量'; 'GEC 售往现货市场'; 'DGC 售往现货市场'}, ...
    [E_PM_grid_0; 0; 0; E_PM_short_0; E_GEC_spot_0; E_DGC_spot_0], ...
    [E_PM_grid_coop; E_PM_GEC_coop; E_PM_DGC_coop; E_PM_short_coop; E_GEC_spot_coop; E_DGC_spot_coop], ...
    'VariableNames', {'Path_Description', 'Baseline_kWh', 'Cooperative_kWh'});
writetable(Energy_Volume, filename_excel, 'Sheet', 'Energy_Volume');

% 3. 制药厂 (PM) 收益分解
PM_Breakdown = table(...
    {'Annual revenue'; 'Penalties from DGC'; 'Capital cost'; 'Grid electricity cost'; 'GEC electricity cost'; 'DGC electricity cost'; 'Shortage penalty cost'; 'Total income'}, ...
    [PM_rev_non; PM_penalty_non; PM_cap_non; PM_cost_grid_non; PM_cost_gec_non; PM_cost_dgc_non; PM_short_cost_non; PM_total_non], ...
    [PM_rev_co; PM_penalty_co; PM_cap_co; PM_cost_grid_co; PM_cost_gec_co; PM_cost_dgc_co; PM_short_cost_co; PM_total_co], ...
    'VariableNames', {'Item_Description', 'Non_Coop', 'Cooperative'});
writetable(PM_Breakdown, filename_excel, 'Sheet', 'PM_Breakdown');

% 4. 绿电公司 (GEC) 收益分解
GEC_Breakdown = table(...
    {'Income from PM'; 'Spot market revenue'; 'Wheeling charge'; 'Total income'}, ...
    [GEC_income_pm_non; GEC_spot_rev_non; GEC_wheel_non; GEC_total_non], ...
    [GEC_income_pm_co; GEC_spot_rev_co; GEC_wheel_co; GEC_total_co], ...
    'VariableNames', {'Item_Description', 'Non_Coop', 'Cooperative'});
writetable(GEC_Breakdown, filename_excel, 'Sheet', 'GEC_Breakdown');

% 5. 柴油公司 (DGC) 收益分解
DGC_Breakdown = table(...
    {'Income from PM'; 'Spot market revenue'; 'Penalty to PM'; 'Total income'}, ...
    [DGC_income_pm_non; DGC_spot_rev_non; DGC_penalty_non; DGC_total_non], ...
    [DGC_income_pm_co; DGC_spot_rev_co; DGC_penalty_co; DGC_total_co], ...
    'VariableNames', {'Item_Description', 'Non_Coop', 'Cooperative'});
writetable(DGC_Breakdown, filename_excel, 'Sheet', 'DGC_Breakdown');

% 6. 同时保存所有的工作区变量到 mat 文件
save(filename_mat);

fprintf('\n>>>数据保存完毕！\n');
fprintf('>>> 1. 结构化明细已保存至 Excel 文件: %s (内含5个Sheet)\n', filename_excel);
fprintf('>>> 2. 完整工作区变量已保存至 MAT 文件: %s\n\n', filename_mat);