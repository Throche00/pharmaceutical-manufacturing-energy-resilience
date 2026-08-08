clc; close all;
fprintf('======================================================\n');
doPlot = true;
if exist('WS_BATCH_MODE', 'var') && WS_BATCH_MODE
    doPlot = false;
end

T = 24; S = 4;
try
    load('Phase2_Clustered_Scenarios.mat', 'Expected_Ptot', 'Expected_Pg', 'Expected_Pd', 'Prob_Ptot', 'Prob_Pg', 'Prob_Pd');
catch
    error('No Phase2_Clustered_Scenarios.mat');
end

season_days = [92,92,91,90];
w_season = season_days / sum(season_days);   

P_demand = zeros(S, T); P_green  = zeros(S, T); P_diesel = zeros(S, T);
for s = 1:S
    P_demand(s,:) = max(0, Prob_Ptot(s,:) * squeeze(Expected_Ptot(s,:,:)));
    P_green(s,:)  = max(0, Prob_Pg(s,:)   * squeeze(Expected_Pg(s,:,:)));
    P_diesel(s,:) = max(0, Prob_Pd(s,:)   * squeeze(Expected_Pd(s,:,:)));
end

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
I_out_prob = 0.04 * ones(1,T); I_out_prob([11:15, 20:24]) = 0.12; 
rng(42); 

I_state = zeros(S, T); 
for s = 1:S
    I_state(s,:) = rand(1, T) < I_out_prob; 
end

P_contract_diesel =900;

P_grid0  = zeros(S, T); P_short0 = zeros(S, T);
T_p_0 = 0; T_g_0 = 0; T_d_0 = 0;
for s = 1:S
    P_grid0(s,:)  = P_demand(s,:) .* (1 - I_state(s,:)); 
    P_short0(s,:) = max(P_demand(s,:) - P_grid0(s,:), 0);
    
    T_p_0 = T_p_0 + w_season(s) * sum(price_prod_value .* P_grid0(s,:) - price_grid .* P_grid0(s,:) - price_shortage .* P_short0(s,:));
    T_g_0 = T_g_0 + w_season(s) * sum(price_spot_green .* P_green(s,:));
    T_d_0 = T_d_0 + w_season(s) * sum(price_spot_diesel .* P_diesel(s,:));
end

%% ==================== 4. P1： ====================
fprintf('--- run P1 ---\n');
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
        fprintf('>>> P1 convergence on the %d time\n', iter);
        Obj_PM_store(iter:end)=Obj_PM_store(iter); Obj_GEC_store(iter:end)=Obj_GEC_store(iter); Obj_DGC_store(iter:end)=Obj_DGC_store(iter);
        r_g_store(iter:end)=r_g_max; r_d_store(iter:end)=r_d_max; break;
    end
end
lastP1 = find(~isnan(Obj_PM_store),1,'last');
delta_p_star = Obj_PM_store(lastP1); delta_g_star = Obj_GEC_store(lastP1); delta_d_star = Obj_DGC_store(lastP1);

fprintf('\n--- run P2 ---\n');
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
        fprintf('>>> P2 convergence on the %d time\n', iter);
        Updp_hist(iter:end)=Updp_hist(iter); Cgp_hist(iter:end)=Cgp_hist(iter); Cdp_hist(iter:end)=Cdp_hist(iter); break;
    end
end
Updp = v_Updp_pm * Scale; Cgp  = v_Cgp_pm * Scale; Cdp  = v_Cdp_pm * Scale;

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

annual_grid_cost_non = 0;
for s = 1:4, annual_grid_cost_non = annual_grid_cost_non + season_days(s) * sum(price_grid .* P_grid0(s,:)); end
T_p_0_ann = price_prod_value * E_PM_grid_0 - annual_grid_cost_non - price_shortage * E_PM_short_0 - C_cap_p;

T_g_0_ann = price_spot_green * E_GEC_spot_0;

T_d_0_ann = price_spot_diesel * E_DGC_spot_0;

annual_grid_cost_co = 0;
for s = 1:4, annual_grid_cost_co = annual_grid_cost_co + season_days(s) * sum(price_grid .* v_Pgrid_pm(s,:)); end
T_p_coop_ann = price_prod_value * (E_PM_grid_c + E_PM_gec_c + E_PM_dgc_c) + (Updp * days_year) ...
               - C_cap_p - annual_grid_cost_co - (Cgp * days_year) - (Cdp * days_year) - (price_shortage * E_PM_short_c);

GEC_wheel_co = 0;
for s = 1:4
    daily_wheel = sum(alpha_tr * (val_Pg2p_gec(s,:).^2) + beta_tr * val_Pg2p_gec(s,:));
    GEC_wheel_co = GEC_wheel_co + daily_wheel * season_days(s);
end
T_g_coop_ann = price_spot_green * E_GEC_spot_c + (Cgp * days_year) - GEC_wheel_co;

T_d_coop_ann = price_spot_diesel * E_DGC_spot_c + (Cdp * days_year) - (Updp * days_year);

T_total_0_ann = T_p_0_ann + T_g_0_ann + T_d_0_ann;
T_total_coop_ann = T_p_coop_ann + T_g_coop_ann + T_d_coop_ann;

if doPlot
    fontN = 'Times New Roman'; 
    seasons_name = {'Spring', 'Summer', 'Autumn', 'Winter'};
    figure('Position',[50 50 900 750],'Color','w','Name','Convergence History'); 

    subplot(2,2,3); 
    plot(Obj_PM_store, 'b.-', 'MarkerSize',8, 'LineWidth',1.2); 
    set(gca,'FontName',fontN,'FontSize',11,'XTick',1:20:maxIter_P1,'XLim',[1 maxIter_P1]); 
    ylabel('\it \delta^{p} \rm ($)'); 
    xlabel({'Iteration', '(a) PM objective function convergence'}, 'FontName', fontN, 'FontWeight', 'normal');
    box on; grid on;

    subplot(2,2,1); 
    plot(Obj_GEC_store, 'b.-', 'MarkerSize',8, 'LineWidth',1.2); 
    set(gca,'FontName',fontN,'FontSize',11,'XTick',1:20:maxIter_P1,'XLim',[1 maxIter_P1]); 
    ylabel('\it \delta^{g} \rm ($)'); 
    xlabel({'Iteration', '(b) GEC objective function convergence'}, 'FontName', fontN, 'FontWeight', 'normal'); 
    box on; grid on;
    
    subplot(2,2,2); 
    plot(Obj_DGC_store, 'b.-', 'MarkerSize',8, 'LineWidth',1.2); 
    set(gca,'FontName',fontN,'FontSize',11,'XTick',1:20:maxIter_P1,'XLim',[1 maxIter_P1]); 
    ylabel('\it \delta^{d} \rm ($)'); 
    xlabel({'Iteration', '(c) DGC objective function convergence'}, 'FontName', fontN, 'FontWeight', 'normal');
    box on; grid on;
    
    subplot(2,2,4); 
    plot(1:maxIter_P2, Cgp_hist, 'b-o', 1:maxIter_P2, Cdp_hist, 'm-d', 1:maxIter_P2, Updp_hist, 'k-*', 'MarkerSize',5, 'LineWidth',1.2); 
    set(gca,'FontName',fontN,'FontSize',11,'XLim',[1 maxIter_P2]); 
    ylabel('Payments 10^4 ($)'); 
    xlabel({'Iteration', '(d) Convergence of bargaining payments'}, 'FontName', fontN, 'FontWeight', 'normal');
    legend('\it C^p_{g}','\it C^p_{d}','\it U^p_{pd}','Location','best'); 
    box on; grid on;

    sub_labels = {'(a)', '(b)', '(c)', '(d)'}; 

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

    
    figure('Position',[350 100 500 400],'Color','w','Name','ADMM Residual Decay (Log Scale)');
    semilogy(1:maxIter_P1, r_g_store, 'g-', 'LineWidth', 2, 'DisplayName', 'Residual: PM-GEC'); hold on; grid on; box on;
    semilogy(1:maxIter_P1, r_d_store, 'r-', 'LineWidth', 2, 'DisplayName', 'Residual: PM-DGC');
    xlabel('Iteration', 'FontName', fontN, 'FontSize', 11); ylabel('Primal Residuals (Log Scale)', 'FontName', fontN, 'FontSize', 11);
    title('ADMM Residual Decay', 'FontName', fontN, 'FontSize', 12);
    legend('Location','northeast', 'FontName', fontN); set(gca, 'FontName', fontN);

    figure('Position',[400 150 500 400],'Color','w','Name','Seasonal Energy Composition');
    energy_gec_seasonal = sum(val_Pg2p_gec, 2) .* season_days'; 
    energy_dgc_seasonal = sum(val_Pd2p_dgc, 2) .* season_days'; 
    b_comp = bar(1:4, [energy_gec_seasonal, energy_dgc_seasonal]/1000, 'stacked');
    b_comp(1).FaceColor = [0.9290 0.6940 0.1250]; 
    b_comp(2).FaceColor = [0.4 0.4 0.4];          
    set(gca, 'XTick', 1:4, 'XTickLabel', seasons_name, 'FontName', fontN, 'FontSize', 11);
    ylabel('Total Traded Energy (MWh)', 'FontName', fontN, 'FontSize', 11);
    title('Energy Supply by Source', 'FontName', fontN, 'FontSize', 12);
    legend('Solar Energy', 'Diesel Backup', 'Location', 'northwest', 'FontName', fontN);
    grid on; box on;

    figure('Position',[450 200 800 650],'Color','w','Name','Trading Heatmaps');
    
    subplot(2,1,1); 
    imagesc(1:24, 1:4, val_Pg2p_gec); colorbar;
    set(gca, 'YTick', 1:4, 'YTickLabel', seasons_name, 'XTick', 1:2:24, 'FontName', fontN, 'FontSize', 11);
    xlabel({'Time Slot (Hour)', '(a) Green energy traded to PM (kWh)'}, ...
        'FontName', fontN, 'FontSize', 11, 'FontWeight', 'normal');
    
    subplot(2,1,2); 
    imagesc(1:24, 1:4, val_Pd2p_dgc); colorbar;
    set(gca, 'YTick', 1:4, 'YTickLabel', seasons_name, 'XTick', 1:2:24, 'FontName', fontN, 'FontSize', 11);
    xlabel({'Time Slot (Hour)', '(b) Diesel energy traded to PM (kWh)'}, ...
        'FontName', fontN, 'FontSize', 11, 'FontWeight', 'normal');

    figure('Position',[500 250 800 400],'Color','w','Name','Hourly Dispatch Profile (Summer)');
    s_target = 2; 
    plot(1:24, v_Pgrid_pm(s_target,:), 'b-s', 'LineWidth', 1.5, 'MarkerFaceColor','w', 'DisplayName', 'Grid Supply'); hold on; grid on; box on;
    plot(1:24, val_Pg2p_gec(s_target,:), 'g-^', 'LineWidth', 1.5, 'MarkerFaceColor','w', 'DisplayName', 'Green Energy');
    plot(1:24, val_Pd2p_dgc(s_target,:), 'r-d', 'LineWidth', 1.5, 'MarkerFaceColor','w', 'DisplayName', 'Diesel Backup');
    plot(1:24, P_demand(s_target,:), 'k--o', 'LineWidth', 2, 'MarkerFaceColor','w', 'DisplayName', 'Total Demand');
    yline(P_contract_diesel, 'r:', 'LineWidth', 2, 'DisplayName', '\it {P_{d2p}^p}'); 
    title('Hourly Dispatch Profile during Summer Peak', 'FontName', fontN, 'FontSize', 13);
    xlabel('Time Slot (Hour)', 'FontName', fontN, 'FontSize', 11);
    ylabel('Power (kW)', 'FontName', fontN, 'FontSize', 11);
    legend('Location','best', 'FontName', fontN, 'FontSize', 11);
end

yalmip('clear');
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

DGC_income_pm_non   = 0;
DGC_income_pm_co    = Cdp * days_year;   
DGC_spot_rev_non    = price_spot_diesel * E_DGC_spot_0;
DGC_spot_rev_co     = price_spot_diesel * E_DGC_spot_c;

DGC_penalty_non     = 0;
DGC_penalty_co      = Updp * days_year;  

DGC_total_non = DGC_income_pm_non + DGC_spot_rev_non - DGC_penalty_non;
DGC_total_co  = DGC_income_pm_co  + DGC_spot_rev_co  - DGC_penalty_co;

fprintf('\n=================（PM）=================\n');
fprintf('%-50s %-15s %-15s\n','Item','Cooperative','Noncooperative');
fprintf('%-50s %-15.2f %-15.2f\n','Annual revenue', PM_rev_co, PM_rev_non);
fprintf('%-50s %-15.2f %-15.2f\n','penalties from the diesel-guarantee company', PM_penalty_co, PM_penalty_non);
fprintf('%-50s %-15.2f %-15.2f\n','capital cost', PM_cap_co, PM_cap_non);
fprintf('%-50s %-15.2f %-15.2f\n','Electricity costs with power main-grid', PM_cost_grid_co, PM_cost_grid_non);
fprintf('%-50s %-15.2f %-15.2f\n','Electricity costs with green-electricity company', PM_cost_gec_co, PM_cost_gec_non);
fprintf('%-50s %-15.2f %-15.2f\n','Electricity costs with diesel-guarantee company', PM_cost_dgc_co, PM_cost_dgc_non);
fprintf('%-50s %-15.2f %-15.2f\n','Shortage penalty cost', PM_short_cost_co, PM_short_cost_non);
fprintf('%-50s %-15.2f %-15.2f\n','Total income', PM_total_co, PM_total_non);

fprintf('\n=================（GEC） =================\n');
fprintf('%-50s %-15s %-15s\n','Item','Cooperative','Noncooperative');
fprintf('%-50s %-15.2f %-15.2f\n','Electricity income with pharmaceutical manufacturer', GEC_income_pm_co, GEC_income_pm_non);
fprintf('%-50s %-15.2f %-15.2f\n','Annual revenue from spot market', GEC_spot_rev_co, GEC_spot_rev_non);
fprintf('%-50s %-15.2f %-15.2f\n','wheeling charge', GEC_wheel_co, GEC_wheel_non);
fprintf('%-50s %-15.2f %-15.2f\n','Total income', GEC_total_co, GEC_total_non);

fprintf('\n=================（DGC） =================\n');
fprintf('%-50s %-15s %-15s\n','Item','Cooperative','Noncooperative');
fprintf('%-50s %-15.2f %-15.2f\n','Electricity income with pharmaceutical manufacturer', DGC_income_pm_co, DGC_income_pm_non);
fprintf('%-50s %-15.2f %-15.2f\n','Annual revenue from spot market', DGC_spot_rev_co, DGC_spot_rev_non);
fprintf('%-50s %-15.2f %-15.2f\n','penalty to the pharmaceutical manufacturer', DGC_penalty_co, DGC_penalty_non);
fprintf('%-50s %-15.2f %-15.2f\n','Total income', DGC_total_co, DGC_total_non);

fprintf('%-35s %-20s %-20s\n','Path','Baseline','Cooperative');
fprintf('%-35s %-20.1f %-20.1f\n','PM purchases electricity from the main grid', E_PM_grid_0, E_PM_grid_c);
fprintf('%-35s %-20.1f %-20.1f\n','PM purchases electricity from GEC', 0, E_PM_gec_c);
fprintf('%-35s %-20.1f %-20.1f\n','PM purchases electricity from DGC', 0, E_PM_dgc_c);
fprintf('%-35s %-20.1f %-20.1f\n','PM electricity shortage', E_PM_short_0, E_PM_short_c);
fprintf('%-35s %-20.1f %-20.1f\n','GEC sells to the spot market', E_GEC_spot_0, E_GEC_spot_c);
fprintf('%-35s %-20.1f %-20.1f\n','DGC sells to the spot market', E_DGC_spot_0, E_DGC_spot_c);
filename_excel = 'P_contract_diesel1300.xlsx';
filename_mat   = 'P_contract_diesel1300.mat';

Profit_Summary = table(...
    {'PM (Pharm)'; 'GEC (Green)'; 'DGC (Diesel)'; 'Coalition Total'}, ...
    [T_p_0_ann; T_g_0_ann; T_d_0_ann; T_p_0_ann+T_g_0_ann+T_d_0_ann], ...
    [T_p_coop_ann; T_g_coop_ann; T_d_coop_ann; T_p_coop_ann+T_g_coop_ann+T_d_coop_ann], ...
    [T_p_coop_ann - T_p_0_ann; T_g_coop_ann - T_g_0_ann; T_d_coop_ann - T_d_0_ann; (T_p_coop_ann+T_g_coop_ann+T_d_coop_ann) - (T_p_0_ann+T_g_0_ann+T_d_0_ann)], ...
    'VariableNames', {'Agent_Name', 'Non_Coop_Profit', 'Coop_Profit', 'Nash_Surplus'});
writetable(Profit_Summary, filename_excel, 'Sheet', 'Profit_Summary');

Energy_Volume = table(...
    {'PM grid purchase'; 'PM purchase from GEC'; 'PM purchase from DGC'; 'PM unmet electricity'; 'GEC spot-market sale'; 'DGC spot-market sale'}, ...
    [E_PM_grid_0; 0; 0; E_PM_short_0; E_GEC_spot_0; E_DGC_spot_0], ...
    [E_PM_grid_c; E_PM_gec_c; E_PM_dgc_c; E_PM_short_c; E_GEC_spot_c; E_DGC_spot_c], ...
    'VariableNames', {'Path_Description', 'Baseline_kWh', 'Cooperative_kWh'});
writetable(Energy_Volume, filename_excel, 'Sheet', 'Energy_Volume');

PM_Breakdown = table(...
    {'Annual revenue'; 'Penalties from DGC'; 'Capital cost'; 'Grid electricity cost'; 'GEC electricity cost'; 'DGC electricity cost'; 'Shortage penalty cost'; 'Total income'}, ...
    [PM_rev_non; PM_penalty_non; PM_cap_non; PM_cost_grid_non; PM_cost_gec_non; PM_cost_dgc_non; PM_short_cost_non; PM_total_non], ...
    [PM_rev_co; PM_penalty_co; PM_cap_co; PM_cost_grid_co; PM_cost_gec_co; PM_cost_dgc_co; PM_short_cost_co; PM_total_co], ...
    'VariableNames', {'Item_Description', 'Non_Coop', 'Cooperative'});
writetable(PM_Breakdown, filename_excel, 'Sheet', 'PM_Breakdown');

GEC_Breakdown = table(...
    {'Income from PM'; 'Spot market revenue'; 'Wheeling charge'; 'Total income'}, ...
    [GEC_income_pm_non; GEC_spot_rev_non; GEC_wheel_non; GEC_total_non], ...
    [GEC_income_pm_co; GEC_spot_rev_co; GEC_wheel_co; GEC_total_co], ...
    'VariableNames', {'Item_Description', 'Non_Coop', 'Cooperative'});
writetable(GEC_Breakdown, filename_excel, 'Sheet', 'GEC_Breakdown');

DGC_Breakdown = table(...
    {'Income from PM'; 'Spot market revenue'; 'Penalty to PM'; 'Total income'}, ...
    [DGC_income_pm_non; DGC_spot_rev_non; DGC_penalty_non; DGC_total_non], ...
    [DGC_income_pm_co; DGC_spot_rev_co; DGC_penalty_co; DGC_total_co], ...
    'VariableNames', {'Item_Description', 'Non_Coop', 'Cooperative'});
writetable(DGC_Breakdown, filename_excel, 'Sheet', 'DGC_Breakdown');

save(filename_mat);

fprintf('\n>>>Data saved！\n');