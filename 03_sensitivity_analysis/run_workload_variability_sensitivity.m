%% Sens_03_Fast.m 
if ~exist('RUN_FROM_MASTER', 'var') || ~RUN_FROM_MASTER
    clc; clear; close all;
else
    close all;
end
fprintf('=== start Phase 3 ===\n');

load('SensData_PH1.mat', 'SensData_PH1');
load('SensData_Monte.mat');

if ~exist('mu_list','var') && exist('mu_scale_list','var')
    mu_list = mu_scale_list;
end
if ~exist('sigma_list','var') && exist('sigma_scale_list','var')
    sigma_list = sigma_scale_list;
end

Mat_Lambda  = zeros(Nm, Ns); Mat_Demand  = zeros(Nm, Ns); Mat_Surplus = zeros(Nm, Ns);
Slice_PM_non = zeros(1, Nm); Slice_PM_coop = zeros(1, Nm);
Slice_GEC_non = zeros(1, Nm); Slice_GEC_coop = zeros(1, Nm);
Slice_DGC_non = zeros(1, Nm); Slice_DGC_coop = zeros(1, Nm);
Slice_Total_non = zeros(1, Nm); Slice_Total_coop = zeros(1, Nm);
slice_idx = 3; 
Slice_E_Grid_coop = zeros(1, Nm); Slice_E_Green_coop = zeros(1, Nm); Slice_E_Diesel_coop = zeros(1, Nm);
Slice_Cgp = zeros(1, Nm); Slice_Cdp = zeros(1, Nm); Slice_Updp = zeros(1, Nm);
Slice_Shortage_non = zeros(1, Nm); Slice_Shortage_coop = zeros(1, Nm);

T = 24; S_energy = 4; days_year = 365; season_days = [92,92,91,90]; w_season = season_days/days_year;
price_grid = [repmat(0.45,1,8), repmat(1.10,1,4), repmat(0.75,1,6), repmat(1.10,1,4), repmat(0.45,1,2)];
price_spot_green = 0.40; price_spot_diesel = 0.90; price_shortage = 2.10; price_prod_value = 1.60;
diesel_reserve_cost_per_kw_day = 0.08;
alpha_tr = 0.0005; beta_tr = 0.004; P_contract_diesel = 900;
Value_per_Mu = 48.5;
mu_base = [5.0, 5.0, 8.0, 8.0, 3.0, 3.0];

hours = 0:T-1;
I_out_prob = 0.06 * ones(1,T);
outage_peak_mask = (hours>=11 & hours<15) | (hours>=20 & hours<24);
I_out_prob(outage_peak_mask) = 0.18;
grid_availability = 1 - I_out_prob;
outage_probability = I_out_prob;
green_delivery_factor = 1 - 0.20 * outage_probability;
ops = sdpsettings('verbose',0, 'warning', 0);

for i = 1:Nm
    for j = 1:Ns
        yalmip('clear');
        
        Mat_Lambda(i,j) = SensData_PH1{i,j}.lambda_star;
        
        P_demand = zeros(S_energy,T); P_green = zeros(S_energy,T); P_diesel = zeros(S_energy,T);
        for s=1:S_energy
            P_demand(s,:) = max(0, SensData_Monte{i,j}.Prob_Joint(s,:) * squeeze(SensData_Monte{i,j}.Expected_Ptot(s,:,:)));
            P_green(s,:)  = max(0, SensData_Monte{i,j}.Prob_Joint(s,:) * squeeze(SensData_Monte{i,j}.Expected_Pg(s,:,:)));
            P_diesel(s,:) = max(0, SensData_Monte{i,j}.Prob_Joint(s,:) * squeeze(SensData_Monte{i,j}.Expected_Pd(s,:,:)));
        end
        P_green_deliverable = P_green .* repmat(green_delivery_factor, S_energy, 1);
        Mat_Demand(i,j) = sum(sum(P_demand, 2) .* season_days');
        
        T_p_0 = 0; T_g_0 = 0; T_d_0 = 0;
        P_short0_all = zeros(S_energy, T); 
        for s = 1:S_energy
            P_grid0 = P_demand(s,:) .* grid_availability;
            P_short0 = max(P_demand(s,:) - P_grid0, 0);
            P_short0_all(s,:) = P_short0;

           
            current_total_mu = sum(mu_base) * mu_list(i); 
            Hourly_Revenue_Cap = Mat_Lambda(i,j) * current_total_mu * Value_per_Mu;

            Nominal_Revenue0 = price_prod_value * P_grid0; 
            Real_Revenue0 = min(Nominal_Revenue0, Hourly_Revenue_Cap); 

            T_p_0 = T_p_0 + w_season(s)*sum(Real_Revenue0 - price_grid.*P_grid0 - price_shortage*P_short0);
            T_g_0 = T_g_0 + w_season(s)*sum(price_spot_green*P_green_deliverable(s,:));
            T_d_0 = T_d_0 + w_season(s)*sum(price_spot_diesel*P_diesel(s,:));
        end
        
        T_coop_total = 0;
        delta_g_star = 0; 
        delta_d_star = 0; 

        Egrid_day = 0; Egreen_day = 0; Ediesel_day = 0;
        Eshort_non_day = 0; Eshort_coop_day = 0;
        
        for s = 1:S_energy
            Pgrid_pm = sdpvar(1,T); Pg_pm = sdpvar(1,T); Pd_pm = sdpvar(1,T); Pshort_pm = sdpvar(1,T);
            Pg_gec = sdpvar(1,T); Pspot_g = sdpvar(1,T);
            Pd_dgc = sdpvar(1,T); Pspot_d = sdpvar(1,T);
            
            Constraints = [
                Pgrid_pm + Pg_pm + Pd_pm + Pshort_pm == P_demand(s,:),
                Pgrid_pm <= P_demand(s,:) .* grid_availability,
                Pg_gec + Pspot_g <= P_green_deliverable(s,:),
                Pd_dgc + Pspot_d <= P_diesel(s,:),
                Pd_dgc <= P_contract_diesel .* outage_probability,
                Pd_dgc <= P_demand(s,:) .* outage_probability,   
                Pg_pm == Pg_gec, 
                Pd_pm == Pd_dgc, 
                Pgrid_pm >= 0, Pg_pm >= 0, Pd_pm >= 0, Pshort_pm >= 0,
                Pg_gec >= 0, Pspot_g >= 0, Pd_dgc >= 0, Pspot_d >= 0
            ];
            
            Nominal_Revenue_Coop = price_prod_value * (Pgrid_pm + Pg_pm + Pd_pm);
            Real_Revenue_Coop = min(Nominal_Revenue_Coop, Hourly_Revenue_Cap); 

            Obj_Total = sum(Real_Revenue_Coop - price_grid.*Pgrid_pm - price_shortage*Pshort_pm) ...
                      + sum(price_spot_green*Pspot_g) - sum(alpha_tr*(Pg_gec.^2) + beta_tr*Pg_gec) ...
                      + sum(price_spot_diesel*Pspot_d) - diesel_reserve_cost_per_kw_day * P_contract_diesel;
                  
            optimize(Constraints, -Obj_Total, ops);
            T_coop_total = T_coop_total + w_season(s) * value(Obj_Total);
            
            val_delta_g = sum(price_spot_green*value(Pspot_g)) - sum(alpha_tr*(value(Pg_gec).^2) + beta_tr*value(Pg_gec));
            val_delta_d = sum(price_spot_diesel*value(Pspot_d)) - diesel_reserve_cost_per_kw_day * P_contract_diesel;
            
            delta_g_star = delta_g_star + w_season(s) * val_delta_g;
            delta_d_star = delta_d_star + w_season(s) * val_delta_d;

            Egrid_day = Egrid_day + w_season(s) * sum(value(Pgrid_pm));
            Egreen_day = Egreen_day + w_season(s) * sum(value(Pg_pm));
            Ediesel_day = Ediesel_day + w_season(s) * sum(value(Pd_pm));
            Eshort_non_day = Eshort_non_day + w_season(s) * sum(P_short0_all(s,:));
            Eshort_coop_day = Eshort_coop_day + w_season(s) * sum(value(Pshort_pm));
        end
        
        Total_Surplus = T_coop_total - (T_p_0 + T_g_0 + T_d_0);
        Mat_Surplus(i,j) = Total_Surplus * days_year;
        
        T_p_coop = T_p_0 + Total_Surplus / 3;
        T_g_coop = T_g_0 + Total_Surplus / 3;
        T_d_coop = T_d_0 + Total_Surplus / 3;
        
        if j == slice_idx
            Slice_PM_non(i) = T_p_0*days_year; Slice_PM_coop(i) = T_p_coop*days_year;
            Slice_GEC_non(i) = T_g_0*days_year; Slice_GEC_coop(i) = T_g_coop*days_year;
            Slice_DGC_non(i) = T_d_0*days_year; Slice_DGC_coop(i) = T_d_coop*days_year;
            Slice_Total_non(i) = Slice_PM_non(i)+Slice_GEC_non(i)+Slice_DGC_non(i);
            Slice_Total_coop(i) = Slice_PM_coop(i)+Slice_GEC_coop(i)+Slice_DGC_coop(i);
            
            Slice_E_Grid_coop(i) = Egrid_day * days_year;
            Slice_E_Green_coop(i) = Egreen_day * days_year;
            Slice_E_Diesel_coop(i) = Ediesel_day * days_year;
            
            Slice_Shortage_non(i) = Eshort_non_day * days_year;
            Slice_Shortage_coop(i) = Eshort_coop_day * days_year;
            
            Slice_Cgp(i) = (Total_Surplus / 3 + T_g_0 - delta_g_star) * days_year;
            Slice_Cdp(i) = (Total_Surplus / 3 + T_d_0 - delta_d_star) * days_year;
            Slice_Updp(i) = 0; 
        end
    end
end


fontN = 'Times New Roman';
figure('Name','Combined Fig 1: Manufacturing-Energy Constraints','Color','w','Position',[50 400 1300 400]);

subplot(1,3,1);
imagesc(sigma_list, mu_list, Mat_Lambda); set(gca, 'YDir', 'normal'); colorbar; colormap('parula');
ylabel('Workload Mean (\mu scale)', 'FontName', fontN, 'FontSize', 12);
xlabel({'Workload Volatility (\sigma scale)', '', '(a) Throughput \lambda^* Sensitivity'}, 'FontName', fontN, 'FontSize', 13);
set(gca,'FontName',fontN); box on;

subplot(1,3,2);
imagesc(sigma_list, mu_list, Mat_Demand); set(gca, 'YDir', 'normal'); colorbar; colormap('parula');
ylabel('Workload Mean (\mu scale)', 'FontName', fontN, 'FontSize', 12);
xlabel({'Workload Volatility (\sigma scale)', '', '(b) Annual Electricity Demand (kW·h)'}, 'FontName', fontN, 'FontSize', 13);
set(gca,'FontName',fontN); box on;

subplot(1,3,3);
yyaxis left; plot(mu_list, Mat_Lambda(:, slice_idx), 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor','w');
ylabel('Theoretical Max Throughput (\lambda^*)', 'FontName', fontN, 'FontSize', 12); set(gca, 'ycolor', 'b');
yyaxis right; plot(mu_list, Mat_Demand(:, slice_idx), 'r-s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor','w');
ylabel('Annual Electricity Demand (kW·h)', 'FontName', fontN, 'FontSize', 12); set(gca, 'ycolor', 'r');
xlabel({['Workload Mean (\mu scale)' num2str(sigma_list(slice_idx)) ']'], '', '(c) Manufacturing-Energy Coupling'}, 'FontName', fontN, 'FontSize', 13);
grid on; box on; set(gca,'FontName',fontN);

figure('Name','Combined Fig 2: Economic & Coalition Analysis','Color','w','Position',[50 200 1300 400]);

subplot(1,3,1);
imagesc(sigma_list, mu_list, Mat_Surplus); set(gca, 'YDir', 'normal'); colorbar; colormap('parula');
ylabel('Workload Mean (\mu scale)', 'FontName', fontN, 'FontSize', 12);
xlabel({'Workload Volatility (\sigma scale)', '', '(a) Nash Coalition Surplus ($)'}, 'FontName', fontN, 'FontSize', 13);
set(gca,'FontName',fontN); box on;

subplot(1,3,2);
plot(mu_list, Slice_PM_coop/1e6, 'b-o', 'LineWidth', 2, 'MarkerFaceColor','w', 'DisplayName', 'PM (Pharm)'); hold on;
plot(mu_list, Slice_GEC_coop/1e6, 'r-s', 'LineWidth', 2, 'MarkerFaceColor','w', 'DisplayName', 'GEC (Green)');
plot(mu_list, Slice_DGC_coop/1e6, 'y-d', 'LineWidth', 2, 'MarkerFaceColor','w', 'DisplayName', 'DGC (Diesel)');
plot(mu_list, Slice_Total_coop/1e6, 'k-^', 'LineWidth', 2.5, 'MarkerFaceColor','k', 'DisplayName', 'Coalition Total');
ylabel('Annual Coop Profit (Million $)', 'FontName', fontN, 'FontSize', 12);
xlabel({'Workload Mean (\mu scale)', '', '(b) Profit Sensitivity to Workload Mean'}, 'FontName', fontN, 'FontSize', 13);
lgd_profit = legend('Location', 'west', 'FontName', fontN, 'FontSize', 10, 'Box', 'on');
lgd_profit.ItemTokenSize = [14, 8];
grid on; box on; set(gca,'FontName',fontN);

subplot(1,3,3);
plot(mu_list, Slice_Cgp/1e6, 'g-o', 'LineWidth', 2, 'MarkerFaceColor','w', 'DisplayName', 'Net Payment to GEC'); hold on;
plot(mu_list, Slice_Cdp/1e6, 'k-s', 'LineWidth', 2, 'MarkerFaceColor','w', 'DisplayName', 'Net Payment to DGC');
ylabel('Net Transfer Amount (Million $)', 'FontName', fontN, 'FontSize', 12);
xlabel({'Workload Mean (\mu scale)', '', '(c) Nash Bargaining Net Transfers'}, 'FontName', fontN, 'FontSize', 13);
lgd_transfer = legend('Location', 'west', 'FontName', fontN, 'FontSize', 10, 'Box', 'on');
lgd_transfer.ItemTokenSize = [14, 8];
grid on; box on; set(gca,'FontName',fontN);

figure('Name','Combined Fig 3: Resilience & Energy Mix','Color','w','Position',[50 50 1000 450]);

subplot(1,2,1);
b_short = bar(mu_list, [Slice_Shortage_non', Slice_Shortage_coop']/1e3, 'grouped', 'EdgeColor', 'k');
b_short(1).FaceColor = [0.8 0.4 0.4];
b_short(2).FaceColor = [0.2 0.6 0.8];
ylabel('Annual Power Shortage (MW·h)', 'FontName', fontN, 'FontSize', 12);
xlabel({'Workload Mean (\mu scale)', '', '(a) Power Shortage'}, 'FontName', fontN, 'FontSize', 13);
legend('Non-cooperative Base', 'Nash Cooperative Mode', 'Location', 'northwest', 'FontName', fontN);
grid on; box on; set(gca,'FontName',fontN);

subplot(1,2,2);
Total_Energy = Slice_E_Grid_coop + Slice_E_Green_coop + Slice_E_Diesel_coop + Slice_Shortage_coop;
Total_Energy = max(Total_Energy, 1e-9);
Mix_Grid = Slice_E_Grid_coop ./ Total_Energy * 100;
Mix_Green = Slice_E_Green_coop ./ Total_Energy * 100;
Mix_Diesel = Slice_E_Diesel_coop ./ Total_Energy * 100;
Mix_Short = Slice_Shortage_coop ./ Total_Energy * 100;

b_mix = bar(mu_list, [Mix_Grid', Mix_Green', Mix_Diesel', Mix_Short'], 'stacked', 'EdgeColor', 'k');
b_mix(1).FaceColor = [0.2 0.5 0.8];
b_mix(2).FaceColor = [0.3 0.7 0.3];
b_mix(3).FaceColor = [0.4 0.4 0.4];
b_mix(4).FaceColor = [0.8 0.2 0.2];
ylabel('Energy Supply Proportion (%)', 'FontName', fontN, 'FontSize', 12);
xlabel({'Workload Mean (\mu scale)', '', '(b) Energy Supply Structure'}, 'FontName', fontN, 'FontSize', 13);
lgd_mix = legend('Main Grid', 'Green Energy', 'Diesel Backup', 'Unmet Shortage', ...
    'Location', 'south', 'NumColumns', 2, 'FontName', fontN, 'FontSize', 10, 'Box', 'on');
lgd_mix.ItemTokenSize = [14, 8];
grid on; box on; set(gca,'FontName',fontN); ylim([0 100]);
