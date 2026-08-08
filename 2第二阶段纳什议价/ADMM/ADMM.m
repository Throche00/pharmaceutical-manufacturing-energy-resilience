%% 第二阶段： ADMM 优化
fprintf('======================================================\n');
fprintf(' 启动 Phase Two: 纳什议价 ADMM 优化\n');
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
price_spot_green  = 0.40;   % $/kWh, GEC spot-market price
price_spot_diesel = 0.90;   % $/kWh, DGC spot-market opportunity price
price_shortage    = 2.10;   % $/kWh, unmet-load penalty
price_prod_value  = 1.60;   % $/kWh-equivalent, fulfilled production value
diesel_reserve_cost_per_kw_day = 0.08; % $/(kW day), reserve-capacity cost
alpha_tr = 0.0005; beta_tr  = 0.004;           

hours = 0:T-1;
I_out_prob = 0.06 * ones(1,T);
outage_peak_mask = (hours>=11 & hours<15) | (hours>=20 & hours<24);
I_out_prob(outage_peak_mask) = 0.18; % outage-risk peak periods
grid_availability = 1 - I_out_prob;
outage_probability = I_out_prob;
green_delivery_factor = 1 - 0.20 * outage_probability;
P_green_deliverable = P_green .* repmat(green_delivery_factor, S, 1);

P_contract_diesel = 900; % (kW) 
solver_opts = sdpsettings('verbose',0);

%% ==================== 3. 非合作基准 ====================
P_grid0  = zeros(S, T); P_short0 = zeros(S, T);
T_p_0 = 0; T_g_0 = 0; T_d_0 = 0;
for s = 1:S
    P_grid0(s,:)  = P_demand(s,:) .* grid_availability; 
    P_short0(s,:) = max(P_demand(s,:) - P_grid0(s,:), 0);
    
    T_p_0 = T_p_0 + w_season(s) * sum(price_prod_value .* P_grid0(s,:) - price_grid .* P_grid0(s,:) - price_shortage .* P_short0(s,:));
    T_g_0 = T_g_0 + w_season(s) * sum(price_spot_green .* P_green_deliverable(s,:));
    T_d_0 = T_d_0 + w_season(s) * sum(price_spot_diesel .* P_diesel(s,:));
end
fprintf('非合作基准预期日利润 (T_p_0, T_g_0, T_d_0):\n  %.4f,  %.4f,  %.4f\n\n', T_p_0, T_g_0, T_d_0);

%% ==================== 4. P1：物理调度分季节解耦极速求解 ====================
fprintf('--- 运行 P1 (四季完全解耦) ---\n');
maxIter_P1 = 80; tol_primal = 0.1; tol_dual = 0.1; rho_p1 = 0.001; % 修正电量级参数避免震荡

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
                 Pgrid_pm <= P_demand(s,:) .* grid_availability, ...
                 Pd_pm <= P_contract_diesel .* outage_probability ]; 
                 
        delta_p_expr = sum(price_prod_value.*(Pgrid_pm+Pg_pm+Pd_pm) - price_grid.*Pgrid_pm - price_shortage.*Pshort_pm);
        Obj_pm = delta_p_expr - sum(lambda_g(s,:).*Pg_pm) - sum(lambda_d(s,:).*Pd_pm) - (rho_p1/2)*sum((Pg_pm - val_Pg2p_gec(s,:)).^2) - (rho_p1/2)*sum((Pd_pm - val_Pd2p_dgc(s,:)).^2);
        optimize(C_pm, -Obj_pm, solver_opts);
        v_Pg_pm(s,:) = value(Pg_pm); v_Pd_pm(s,:) = value(Pd_pm); v_Pgrid_pm(s,:) = value(Pgrid_pm); v_Pshort_pm(s,:) = value(Pshort_pm); val_delta_p(s) = value(delta_p_expr);

        % GEC 
        Pg_gec = sdpvar(1,T); Pspot_g = sdpvar(1,T);
        C_gec = [Pg_gec + Pspot_g <= P_green_deliverable(s,:), Pg_gec >= 0, Pspot_g >= 0];
        delta_g_expr = sum(price_spot_green.*Pspot_g) - sum(alpha_tr*(Pg_gec.^2) + beta_tr*Pg_gec);
        Obj_gec = delta_g_expr + sum(lambda_g(s,:).*Pg_gec) - (rho_p1/2)*sum((v_Pg_pm(s,:) - Pg_gec).^2);
        optimize(C_gec, -Obj_gec, solver_opts);
        val_Pg2p_gec(s,:) = value(Pg_gec); v_Pspot_g(s,:) = value(Pspot_g); val_delta_g(s) = value(delta_g_expr);

        % DGC 
        Pd_dgc = sdpvar(1,T); Pspot_d = sdpvar(1,T);
        C_dgc = [Pd_dgc + Pspot_d <= P_diesel(s,:), Pd_dgc >= 0, Pspot_d >= 0, ...
                 Pd_dgc <= P_contract_diesel .* outage_probability];
        delta_d_expr = sum(price_spot_diesel.*Pspot_d) - diesel_reserve_cost_per_kw_day * P_contract_diesel;
        Obj_dgc = delta_d_expr + sum(lambda_d(s,:).*Pd_dgc) - (rho_p1/2)*sum((v_Pd_pm(s,:) - Pd_dgc).^2);
        optimize(C_dgc, -Obj_dgc, solver_opts);
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
    optimize([surplus_p >= 1e-3, Updp_pm >= 0, Cgp_pm >= 0, Cdp_pm >= 0], -Obj2_pm, solver_opts);
    v_Updp_pm = value(Updp_pm); v_Cgp_pm = value(Cgp_pm); v_Cdp_pm = value(Cdp_pm);

    Cgp_gec = sdpvar(1,1); surplus_g = dg_s + Cgp_gec - Tg0_s;
    Obj2_gec = log(surplus_g) + gamma_cg*Cgp_gec - (rho_p2/2)*(v_Cgp_pm - Cgp_gec)^2;
    optimize([surplus_g >= 1e-3, Cgp_gec >= 0], -Obj2_gec, solver_opts); val_Cgp_gec = value(Cgp_gec);

    Updp_dgc = sdpvar(1,1); Cdp_dgc = sdpvar(1,1); surplus_d = dd_s - Updp_dgc + Cdp_dgc - Td0_s;
    Obj2_dgc = log(surplus_d) + gamma_upd*Updp_dgc + gamma_cd*Cdp_dgc - (rho_p2/2)*(v_Updp_pm - Updp_dgc)^2 - (rho_p2/2)*(v_Cdp_pm - Cdp_dgc)^2;
    optimize([surplus_d >= 1e-3, Updp_dgc >= 0, Cdp_dgc >= 0], -Obj2_dgc, solver_opts); val_Updp_dgc = value(Updp_dgc); val_Cdp_dgc = value(Cdp_dgc);

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

% --- 6.1  ---
% 计算 PM 资本成本 (CRF)
S_machines = 12; 
theta_s = [120000, 85000, 150000, 60000, 90000, 110000, 75000, 130000, 95000, 100000, 140000, 105000]; 
Y_s     = [10, 12, 10, 15, 10, 10, 12, 15, 12, 12, 10, 12]; 
r       = 0.06; 
C_cap_p = 0;
for i = 1:S_machines
    CRF = (r * (1 + r)^Y_s(i)) / ((1 + r)^Y_s(i) - 1);
    C_cap_p = C_cap_p + CRF * theta_s(i);
end

% 计算年度能量统计 (kWh/year)
E_PM_grid_0  = sum(season_days .* sum(P_grid0, 2)');
E_PM_short_0 = sum(season_days .* sum(P_short0, 2)');
E_GEC_spot_0 = sum(season_days .* sum(P_green_deliverable, 2)');
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
DGC_reserve_cost_non = 0;
DGC_reserve_cost_co = diesel_reserve_cost_per_kw_day * P_contract_diesel * days_year;
T_d_coop_ann = price_spot_diesel * E_DGC_spot_c + (Cdp * days_year) - (Updp * days_year) - DGC_reserve_cost_co;

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
%% ==================== 7. 全套可视化====================
if doPlot
    fontN = 'Times New Roman'; 
    seasons_name = {'Spring', 'Summer', 'Autumn', 'Winter'};
    figure('Position',[50 50 900 750],'Color','w','Name','Convergence History'); % 

    % --- (a) PM 目标函数收敛 ---
    subplot(2,2,1); 
    plot(Obj_PM_store, 'b.-', 'MarkerSize',8, 'LineWidth',1.2); 
    set(gca,'FontName',fontN,'FontSize',11,'XTick',1:20:maxIter_P1,'XLim',[1 maxIter_P1]); 
    ylabel('\it \delta^{p} \rm ($)'); 
    xlabel({'Iteration', '(a) PM objective function convergence'}, 'FontName', fontN, 'FontWeight', 'normal');
    box on; grid on;

    % --- (b) GEC 目标函数收敛 ---
    subplot(2,2,2); 
    plot(Obj_GEC_store, 'b.-', 'MarkerSize',8, 'LineWidth',1.2); 
    set(gca,'FontName',fontN,'FontSize',11,'XTick',1:20:maxIter_P1,'XLim',[1 maxIter_P1]); 
    ylabel('\it \delta^{g} \rm ($)');
    xlabel({'Iteration', '(b) GEC objective function convergence'}, 'FontName', fontN, 'FontWeight', 'normal'); 
    box on; grid on;
    
    % --- (c) DGC 目标函数收敛 ---
    subplot(2,2,3); 
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
        
        ylabel('Electricity (kW·h)', 'FontName', fontN);
        xlabel({'Time (h)', [sub_labels{s}, ' PM Non-Coop - ', seasons_name{s}]}, 'FontName', fontN, 'FontWeight', 'normal');
        
        if s == 1
            lgd = legend('From Main Grid', 'Unmet Demand', 'Total Demand', ...
                'NumColumns', 3, 'FontSize', 6, 'Box', 'on');
            lgd.ItemTokenSize = [7, 4];
            set(lgd, 'Units', 'normalized', 'Position', [0.34, 0.470, 0.32, 0.03]);
        end
    end

    % ================= 2. PM 合作状态能量获取 =================
    figure('Position',[150 150 1000 650],'Color','w','Name','PM Electricity Sourcing (Cooperative)');
    for s = 1:4
        subplot(2,2,s); 
        bar(1:24, [v_Pgrid_pm(s,:)', v_Pg_pm(s,:)', v_Pd_pm(s,:)', v_Pshort_pm(s,:)'], 'stacked'); hold on; 
        plot(1:24, P_demand(s,:), 'k-o', 'LineWidth', 1.5, 'MarkerSize', 4); 
        xlim([0.5 24.5]); box on; set(gca,'FontName',fontN);
        
        ylabel('Electricity (kW·h)', 'FontName', fontN);
        xlabel({'Time (h)', [sub_labels{s}, ' PM Coop - ', seasons_name{s}]}, 'FontName', fontN, 'FontWeight', 'normal');
        
        if s == 1
            lgd = legend('From Main Grid', 'From GEC', 'From DGC', 'Shortage', 'Total Demand', ...
                'NumColumns', 5, 'FontSize', 6, 'Box', 'on');
            lgd.ItemTokenSize = [7, 4];
            set(lgd, 'Units', 'normalized', 'Position', [0.24, 0.470, 0.52, 0.03]);
        end
    end
    
    % ================= 3. GEC 合作状态绿电调度 =================
    figure('Position',[200 200 1000 650],'Color','w','Name','GEC Dispatch (Cooperative)');
    for s = 1:4
        subplot(2,2,s); 
        bar(1:24, [val_Pg2p_gec(s,:)', v_Pspot_g(s,:)'], 'stacked'); hold on; 
        plot(1:24, P_green_deliverable(s,:), 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4); 
        xlim([0.5 24.5]); box on; set(gca,'FontName',fontN);
        
        ylabel('Electricity (kW·h)', 'FontName', fontN);
        xlabel({'Time (h)', [sub_labels{s}, ' GEC Coop - ', seasons_name{s}]}, 'FontName', fontN, 'FontWeight', 'normal');
        
        if s == 1
            lgd = legend('Sold to PM', 'Sold to Spot Market', 'Total Green Output', ...
                'NumColumns', 3, 'FontSize', 6, 'Box', 'on');
            lgd.ItemTokenSize = [7, 4];
            set(lgd, 'Units', 'normalized', 'Position', [0.32, 0.470, 0.36, 0.03]);
        end
    end
    
    % ================= 4. DGC 合作状态柴油调度 =================
    figure('Position',[250 250 1000 650],'Color','w','Name','DGC Dispatch (Cooperative)');
    for s = 1:4
        subplot(2,2,s); 
        bar(1:24, [val_Pd2p_dgc(s,:)', v_Pspot_d(s,:)'], 'stacked'); hold on; 
        plot(1:24, P_diesel(s,:), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4); 
        xlim([0.5 24.5]); box on; set(gca,'FontName',fontN);
        
        ylabel('Electricity (kW·h)', 'FontName', fontN);
        xlabel({'Time (h)', [sub_labels{s}, ' DGC Coop - ', seasons_name{s}]}, 'FontName', fontN, 'FontWeight', 'normal');
        
        if s == 1
            lgd = legend('Guaranteed to PM', 'Sold to Spot Market', 'Total Diesel Capacity', ...
                'NumColumns', 3, 'FontSize', 6, 'Box', 'on');
            lgd.ItemTokenSize = [7, 4];
            set(lgd, 'Units', 'normalized', 'Position', [0.30, 0.470, 0.40, 0.03]);
        end
    end

    figure('Position',[300 300 650 450],'Color','w','Name','Annual Profit Comparison');
    bar_data = [T_p_0_ann, T_p_coop_ann; T_g_0_ann, T_g_coop_ann; T_d_0_ann, T_d_coop_ann; (T_p_0_ann+T_g_0_ann+T_d_0_ann), (T_p_coop_ann+T_g_coop_ann+T_d_coop_ann)] / 10000; 
    b = bar(bar_data, 'FaceColor','flat', 'EdgeColor', 'k', 'LineWidth', 1);
    b(1).CData = [0.7 0.7 0.7]; b(2).CData = [0.2 0.6 0.8]; 
    set(gca, 'XTickLabel', {'PM', 'GEC', 'DGC', 'Total Coalition'}, 'FontName', fontN, 'FontSize', 12); ylabel('Annual Profit (10k $)', 'FontName', fontN, 'FontSize', 12);
    title(''); legend('Non-coop', 'Cooperative', 'Location', 'northwest', 'FontName', fontN); box on; grid on;

    %% ==================== 8.5 新增：对标顶级期刊的增强可视化 ====================
    
    % --- 新图 1: 对数坐标系下的残差收敛 ---
    figure('Position',[350 100 500 400],'Color','w','Name','ADMM Residual Decay (Log Scale)');
    semilogy(1:maxIter_P1, r_g_store, 'g-', 'LineWidth', 2, 'DisplayName', 'Residual: PM-GEC'); hold on; grid on; box on;
    semilogy(1:maxIter_P1, r_d_store, 'r-', 'LineWidth', 2, 'DisplayName', 'Residual: PM-DGC');
    xlabel('Iteration', 'FontName', fontN, 'FontSize', 11); ylabel('Primal Residuals (Log Scale)', 'FontName', fontN, 'FontSize', 11);
    title('ADMM Residual Decay', 'FontName', fontN, 'FontSize', 12);
    legend('Location','northeast', 'FontName', fontN); set(gca, 'FontName', fontN);

    % --- 新图 2: 四季协作能源构成堆叠图  ---
    figure('Position',[400 150 500 400],'Color','w','Name','Seasonal Energy Composition');
    energy_gec_seasonal = sum(val_Pg2p_gec, 2) .* season_days'; 
    energy_dgc_seasonal = sum(val_Pd2p_dgc, 2) .* season_days'; 
    b_comp = bar(1:4, [energy_gec_seasonal, energy_dgc_seasonal]/1000, 'stacked');
    b_comp(1).FaceColor = [0.9290 0.6940 0.1250]; % 太阳能(黄色)
    b_comp(2).FaceColor = [0.4 0.4 0.4];          % 柴油(灰色)
    set(gca, 'XTick', 1:4, 'XTickLabel', seasons_name, 'FontName', fontN, 'FontSize', 11);
    ylabel('Total Traded Energy (MW·h)', 'FontName', fontN, 'FontSize', 11);
    title('Energy Supply by Source', 'FontName', fontN, 'FontSize', 12);
    legend('Solar Energy', 'Diesel Backup', 'Location', 'northwest', 'FontName', fontN);
    grid on; box on;

       % --- 新图 3: 时空全景热力图  ---
    figure('Position',[450 200 800 650],'Color','w','Name','Trading Heatmaps');
    
    % --- 子图 (a) 绿电交易热力图 ---
    subplot(2,1,1); 
    imagesc(1:24, 1:4, val_Pg2p_gec); colorbar;
    set(gca, 'YTick', 1:4, 'YTickLabel', seasons_name, 'XTick', 1:2:24, 'FontName', fontN, 'FontSize', 11);
    xlabel({'Time (h)', '(a) Green energy traded to PM (kW·h)'}, ...
        'FontName', fontN, 'FontSize', 11, 'FontWeight', 'normal');
    
    % --- 子图 (b) 柴油交易热力图 ---
    subplot(2,1,2); 
    imagesc(1:24, 1:4, val_Pd2p_dgc); colorbar;
    set(gca, 'YTick', 1:4, 'YTickLabel', seasons_name, 'XTick', 1:2:24, 'FontName', fontN, 'FontSize', 11);
    xlabel({'Time (h)', '(b) Diesel energy traded to PM (kW·h)'}, ...
        'FontName', fontN, 'FontSize', 11, 'FontWeight', 'normal');

    % --- 新图 4: 典型日(夏季)小时级精细化调度曲线 (严格对标参考图15) ---
    figure('Position',[500 250 800 400],'Color','w','Name','Hourly Dispatch Profile (Summer)');
    s_target = 2; % 选取夏季展示
    plot(1:24, v_Pgrid_pm(s_target,:), 'b-s', 'LineWidth', 1.5, 'MarkerFaceColor','w', 'DisplayName', 'Grid Supply'); hold on; grid on; box on;
    plot(1:24, val_Pg2p_gec(s_target,:), 'g-^', 'LineWidth', 1.5, 'MarkerFaceColor','w', 'DisplayName', 'Green Energy');
    plot(1:24, val_Pd2p_dgc(s_target,:), 'r-d', 'LineWidth', 1.5, 'MarkerFaceColor','w', 'DisplayName', 'Diesel Backup');
    plot(1:24, P_demand(s_target,:), 'k--o', 'LineWidth', 2, 'MarkerFaceColor','w', 'DisplayName', 'Total Demand');
    yline(P_contract_diesel, 'r:', 'LineWidth', 2, 'DisplayName', '\it {P_{d2p}^p}'); % 凸显契约上限
    title('');
    xlabel('Time (h)', 'FontName', fontN, 'FontSize', 11);
    ylabel('Power (kW)', 'FontName', fontN, 'FontSize', 11);
    legend('Location','best', 'FontName', fontN, 'FontSize', 11);
end

yalmip('clear');
%% ==================== 9. 明细分解输出（按论文表格口径） ====================
% ---------- 年度能量统计 ----------
E_PM_grid_0   = sum(season_days .* sum(P_grid0, 2)');
E_PM_short_0  = sum(season_days .* sum(P_short0, 2)');
E_GEC_spot_0  = sum(season_days .* sum(P_green_deliverable, 2)');
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

S_machines = 12; 
theta_s = [120000, 85000, 150000, 60000, 90000, 110000, 75000, 130000, 95000, 100000, 140000, 105000]; 
Y_s     = [10, 12, 10, 15, 10, 10, 12, 15, 12, 12, 10, 12]; 
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
DGC_reserve_cost_non = 0;
DGC_reserve_cost_co  = diesel_reserve_cost_per_kw_day * P_contract_diesel * days_year;

DGC_total_non = DGC_income_pm_non + DGC_spot_rev_non - DGC_penalty_non - DGC_reserve_cost_non;
DGC_total_co  = DGC_income_pm_co  + DGC_spot_rev_co  - DGC_penalty_co  - DGC_reserve_cost_co;

%% ==================== 10. 控制台表格打印 ====================
fprintf('\n================= 制药厂（PM）收益分解 ($/year) =================\n');
fprintf('%-50s %-15s %-15s\n','Item','Cooperative','Noncooperative');
fprintf('%-50s %-15.2f %-15.2f\n','Annual revenue', PM_rev_co, PM_rev_non);
fprintf('%-50s %-15.2f %-15.2f\n','penalties from the diesel-guarantee company', PM_penalty_co, PM_penalty_non);
fprintf('%-50s %-15.2f %-15.2f\n','capital cost', PM_cap_co, PM_cap_non);
fprintf('%-50s %-15.2f %-15.2f\n','Electricity costs with power main-grid', PM_cost_grid_co, PM_cost_grid_non);
fprintf('%-50s %-15.2f %-15.2f\n','Electricity costs with green-electricity company', PM_cost_gec_co, PM_cost_gec_non);
fprintf('%-50s %-15.2f %-15.2f\n','Nash payment to diesel-guarantee company', PM_cost_dgc_co, PM_cost_dgc_non);
fprintf('%-50s %-15.2f %-15.2f\n','Shortage penalty cost', PM_short_cost_co, PM_short_cost_non);
fprintf('%-50s %-15.2f %-15.2f\n','Net profit', PM_total_co, PM_total_non);

fprintf('\n================= 绿电公司（GEC）收益分解 ($/year) =================\n');
fprintf('%-50s %-15s %-15s\n','Item','Cooperative','Noncooperative');
fprintf('%-50s %-15.2f %-15.2f\n','Electricity income with pharmaceutical manufacturer', GEC_income_pm_co, GEC_income_pm_non);
fprintf('%-50s %-15.2f %-15.2f\n','Annual revenue from spot market', GEC_spot_rev_co, GEC_spot_rev_non);
fprintf('%-50s %-15.2f %-15.2f\n','wheeling charge', GEC_wheel_co, GEC_wheel_non);
fprintf('%-50s %-15.2f %-15.2f\n','Net profit', GEC_total_co, GEC_total_non);

fprintf('\n================= 柴油公司（DGC）收益分解 ($/year) =================\n');
fprintf('%-50s %-15s %-15s\n','Item','Cooperative','Noncooperative');
fprintf('%-50s %-15.2f %-15.2f\n','Electricity income with pharmaceutical manufacturer', DGC_income_pm_co, DGC_income_pm_non);
fprintf('%-50s %-15.2f %-15.2f\n','Annual revenue from spot market', DGC_spot_rev_co, DGC_spot_rev_non);
fprintf('%-50s %-15.2f %-15.2f\n','penalty to the pharmaceutical manufacturer', DGC_penalty_co, DGC_penalty_non);
fprintf('%-50s %-15.2f %-15.2f\n','reserve capacity cost', DGC_reserve_cost_co, DGC_reserve_cost_non);
fprintf('%-50s %-15.2f %-15.2f\n','Net profit', DGC_total_co, DGC_total_non);

fprintf('\n================= 交易路径/主体行为 (kWh/year) =================\n');
fprintf('%-35s %-20s %-20s\n','Path','Baseline','Cooperative');
fprintf('%-35s %-20.1f %-20.1f\n','PM 从主网购电', E_PM_grid_0, E_PM_grid_c);
fprintf('%-35s %-20.1f %-20.1f\n','PM 从 GEC 购电', 0, E_PM_gec_c);
fprintf('%-35s %-20.1f %-20.1f\n','PM 从 DGC 购电', 0, E_PM_dgc_c);
fprintf('%-35s %-20.1f %-20.1f\n','PM 未满足缺电量', E_PM_short_0, E_PM_short_c);
fprintf('%-35s %-20.1f %-20.1f\n','GEC 售往现货市场', E_GEC_spot_0, E_GEC_spot_c);
fprintf('%-35s %-20.1f %-20.1f\n','DGC 售往现货市场', E_DGC_spot_0, E_DGC_spot_c);
%% ==================== 追加：各主体细分收益/成本结构对比图 ====================
filename_excel = 'Baseline_ADMM.xlsx';
filename_mat = 'Baseline_ADMM.mat';

Profit_Summary = table(...
    {'PM (Pharm)'; 'GEC (Green)'; 'DGC (Diesel)'; 'Coalition Total'}, ...
    [T_p_0_ann; T_g_0_ann; T_d_0_ann; T_total_0_ann], ...
    [T_p_coop_ann; T_g_coop_ann; T_d_coop_ann; T_total_coop_ann], ...
    [T_p_coop_ann - T_p_0_ann; T_g_coop_ann - T_g_0_ann; T_d_coop_ann - T_d_0_ann; T_total_coop_ann - T_total_0_ann], ...
    'VariableNames', {'Agent_Name', 'Non_Coop_Profit', 'Coop_Profit', 'Nash_Surplus'});
writetable(Profit_Summary, filename_excel, 'Sheet', 'Profit_Summary');

Energy_Volume = table(...
    {'PM grid purchase'; 'PM purchase from GEC'; 'PM purchase from DGC'; 'PM unmet electricity'; 'GEC spot-market sale'; 'DGC spot-market sale'}, ...
    [E_PM_grid_0; 0; 0; E_PM_short_0; E_GEC_spot_0; E_DGC_spot_0], ...
    [E_PM_grid_c; E_PM_gec_c; E_PM_dgc_c; E_PM_short_c; E_GEC_spot_c; E_DGC_spot_c], ...
    'VariableNames', {'Path_Description', 'Baseline_kWh', 'Cooperative_kWh'});
writetable(Energy_Volume, filename_excel, 'Sheet', 'Energy_Volume');

PM_Breakdown = table(...
    {'Annual revenue'; 'Penalties from DGC'; 'Capital cost'; 'Grid electricity cost'; 'GEC electricity cost'; 'Nash payment to DGC'; 'Shortage penalty cost'; 'Net profit'}, ...
    [PM_rev_non; PM_penalty_non; PM_cap_non; PM_cost_grid_non; PM_cost_gec_non; PM_cost_dgc_non; PM_short_cost_non; PM_total_non], ...
    [PM_rev_co; PM_penalty_co; PM_cap_co; PM_cost_grid_co; PM_cost_gec_co; PM_cost_dgc_co; PM_short_cost_co; PM_total_co], ...
    'VariableNames', {'Item_Description', 'Non_Coop', 'Cooperative'});
writetable(PM_Breakdown, filename_excel, 'Sheet', 'PM_Breakdown');

GEC_Breakdown = table(...
    {'Income from PM'; 'Spot market revenue'; 'Wheeling charge'; 'Net profit'}, ...
    [GEC_income_pm_non; GEC_spot_rev_non; GEC_wheel_non; GEC_total_non], ...
    [GEC_income_pm_co; GEC_spot_rev_co; GEC_wheel_co; GEC_total_co], ...
    'VariableNames', {'Item_Description', 'Non_Coop', 'Cooperative'});
writetable(GEC_Breakdown, filename_excel, 'Sheet', 'GEC_Breakdown');

DGC_Breakdown = table(...
    {'Income from PM'; 'Spot market revenue'; 'Penalty to PM'; 'Reserve capacity cost'; 'Net profit'}, ...
    [DGC_income_pm_non; DGC_spot_rev_non; DGC_penalty_non; DGC_reserve_cost_non; DGC_total_non], ...
    [DGC_income_pm_co; DGC_spot_rev_co; DGC_penalty_co; DGC_reserve_cost_co; DGC_total_co], ...
    'VariableNames', {'Item_Description', 'Non_Coop', 'Cooperative'});
writetable(DGC_Breakdown, filename_excel, 'Sheet', 'DGC_Breakdown');

save(filename_mat);
fprintf('\n>>> Baseline tables saved to %s and %s.\n', filename_excel, filename_mat);

scale = 1e6;

% 1. PM (制药厂) 数据提取
% 结构：[生产收入, 主网电费, 绿电费, 柴油费, 缺电惩罚, 柴油商赔偿金, 资本成本]
pm_non_data  = [PM_rev_non, PM_cost_grid_non, PM_cost_gec_non, PM_cost_dgc_non, PM_short_cost_non, PM_penalty_non, PM_cap_non] / scale;
pm_coop_data = [PM_rev_co,  PM_cost_grid_co,  PM_cost_gec_co,  PM_cost_dgc_co,  PM_short_cost_co,  PM_penalty_co,  PM_cap_co] / scale;

% 2. GEC (绿电商) 数据提取
% 结构：[现货市场收入, 售予PM收入, 过网费]
gec_non_data  = [GEC_spot_rev_non, GEC_income_pm_non, GEC_wheel_non] / scale;
gec_coop_data = [GEC_spot_rev_co,  GEC_income_pm_co,  GEC_wheel_co] / scale;

% 3. DGC (柴油商) 数据提取
% 结构：[现货市场收入, 售予PM收入, 支付给PM的违约金]
dgc_non_data  = [DGC_spot_rev_non, DGC_income_pm_non, DGC_penalty_non, DGC_reserve_cost_non] / scale;
dgc_coop_data = [DGC_spot_rev_co,  DGC_income_pm_co,  DGC_penalty_co, DGC_reserve_cost_co] / scale;

% 创建宽幅画布，包含1行3列的子图
figure('Name', 'Financial Breakdown Comparison', 'Color', 'w', 'Position', [100, 200, 1400, 500]);
fontN = 'Times New Roman';
color_non = [0.8 0.4 0.4]; % 非合作 - 浅红色
color_coop = [0.2 0.6 0.8]; % 合作 - 浅蓝色

% --- 子图 1: 制药厂 (PM) 的账本 ---
subplot(1, 3, 1);
pm_bar = bar([pm_non_data', pm_coop_data'], 'grouped', 'EdgeColor', 'k');
pm_bar(1).FaceColor = color_non; pm_bar(2).FaceColor = color_coop;
set(gca, 'XTickLabel', {'\it U^p_{fd}', '\it C^p_{p}', '\it C^p_{g}', '\it C^p_{d}', '\it C^p_{sp}', '\it U^p_{pd}', '\it C^p_{c}'}, ...
    'FontName', fontN, 'FontSize', 9);
ylabel('Income and Cost (Million)', 'FontName', fontN, 'FontSize', 12);
xlabel({'(a) PM financial breakdown structure'}, 'FontName', fontN, 'FontSize', 12);
legend('Non-Cooperative', 'Cooperative', 'Location', 'northwest', 'FontName', fontN, 'FontSize', 9);
grid on; box on;

% --- 子图 2: 绿电 (GEC) 的账本 ---
subplot(1, 3, 2);
gec_bar = bar([gec_non_data', gec_coop_data'], 'grouped', 'EdgeColor', 'k');
gec_bar(1).FaceColor = color_non; gec_bar(2).FaceColor = color_coop;
set(gca, 'XTickLabel', {'\it U^{g}_{m}', '\it U^{g}_{s}', '\it C^{g}_{tr}'}, 'FontName', fontN, 'FontSize', 10);
ylabel('Income and Cost (Million)', 'FontName', fontN, 'FontSize', 12);
xlabel({'(b) GEC financial breakdown structure'}, 'FontName', fontN, 'FontSize', 12);
legend('Non-Cooperative', 'Cooperative', 'Location', 'northeast', 'FontName', fontN, 'FontSize', 9);
grid on; box on;

% --- 子图 3: 柴油 (DGC) 的账本 ---
subplot(1, 3, 3);
dgc_bar = bar([dgc_non_data', dgc_coop_data'], 'grouped', 'EdgeColor', 'k');
dgc_bar(1).FaceColor = color_non; dgc_bar(2).FaceColor = color_coop;
set(gca, 'XTickLabel', {'\it U^{d}_{m}', '\it C^{p}_{d}', '\it U^{p}_{pd}', '\it C^{d}_{res}'}, 'FontName', fontN, 'FontSize', 10);
ylabel('Income and Cost (Million)', 'FontName', fontN, 'FontSize', 12);
xlabel({'(c) DGC financial breakdown structure'}, 'FontName', fontN, 'FontSize', 12);
legend('Non-Cooperative', 'Cooperative', 'Location', 'northeast', 'FontName', fontN, 'FontSize', 9);
grid on; box on;
