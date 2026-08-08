clc; clear; close all;

x_labels = {'2%, 8%', '4%, 12%', '6%, 18%', '8%, 22%', '10%, 25%'};
x_vals = 1:5;

nash_surplus = [2410077.15, 3193936.69, 3991052.13, 5452434.87, 6806131.89];

pm_0 = [7595363.35, 6800332.16, 5941299.81, 4423760.09, 3007885.24];
pm_c = [8398846.55, 7865339.38, 7273019.965, 6244574.32, 5283263.06];
gec_0 = [3312014.33, 3312014.33, 3312014.33, 3312014.33, 3312014.33];
gec_c = [4114995.37, 4375623.74, 4640386.45, 5126608.38, 5577193.40];
dgc_0 = [8683577.22, 8683577.22, 8683577.221, 8683577.22, 8683577.22];
dgc_c = [9487168.97, 9748890.28, 10014491.56, 10500501.06, 10949117.93];

shortage_0 = [573556.48, 857534.84, 1138113.777, 1652619.72, 2130434.69];
shortage_c = [0.01, 0.03, 0.03, 0.03, 0.02]; 

energy_grid_c = [7597277.44, 7409156.93, 7200919.17, 6822500.29, 6461574.97];
energy_gec_c = [2910878.16, 2939808.72, 3004389.09, 3114162.68, 3169004.24];
energy_dgc_c = [322388.20, 481569.27, 625250.76, 893931.80, 1199972.79];

fontN = 'Times New Roman';
set(groot,'defaultAxesFontName',fontN,'defaultAxesFontSize',12);

coalition_0 = pm_0 + gec_0 + dgc_0; 
coalition_c = pm_c + gec_c + dgc_c; 
pay_to_gec = [2778956.08, 3068861.07, 3371304.531, 3936102.07, 4421021.74];
pay_to_dgc = [1095371.02, 1499130.32, 1893673.178, 2621492.39, 3345543.93];

figure('Name', 'Combined Outage Probability Sensitivity', 'Color', 'w', 'Position', [100, 100, 1600, 500]);

%% ==================== (a) ====================
subplot(1, 3, 1);
yyaxis left;
plot(x_vals, coalition_c/1e6, 'k-^', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor','k'); hold on;
plot(x_vals, coalition_0/1e6, 'k--', 'LineWidth', 1.5);
plot(x_vals, pm_c/1e6, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor','w');
plot(x_vals, pm_0/1e6, 'b--', 'LineWidth', 1.5);
ylabel('Profit (Million $)'); 
set(gca, 'ycolor', 'k'); ylim([0, 30]);

yyaxis right;
b = bar(x_vals, [shortage_0'/1e6, shortage_c'/1e6], 0.6, 'EdgeColor', 'k');
b(1).FaceColor = [0.8 0.4 0.4]; b(1).FaceAlpha = 0.7; 
b(2).FaceColor = [0.2 0.6 0.8]; b(2).FaceAlpha = 0.7; 
ylabel('Unmet Power Shortage (GWh)'); 
set(gca, 'ycolor', [0.8 0.3 0.3]); ylim([0, 2.5]);

xticks(x_vals); xticklabels(x_labels);
xlabel({'Grid Outage Probability', '(a) System Resilience against Grid Outage'}, ...
    'FontName', fontN, 'FontSize', 12, 'FontWeight', 'normal');
legend({'Coalition Profit (Coop)', 'Coalition Profit (Non-Coop)', 'PM Profit (Coop)', 'PM Profit (Non-Coop)', ...
        'Shortage (Non-Coop)', 'Shortage (Coop)'}, 'Location', 'southwest', 'NumColumns', 2);
grid on; box on; xlim([0.5, 5.5]);

%% ==================== (b) ====================
subplot(1, 3, 2);
energy_matrix = [energy_dgc_c', energy_gec_c', energy_grid_c'] / 1e6; 
a = area(x_vals, energy_matrix);
a(1).FaceColor = [0.5 0.5 0.5]; 
a(2).FaceColor = [0.3 0.7 0.4]; 
a(3).FaceColor = [0.2 0.5 0.8]; 
a(1).FaceAlpha = 0.8; a(2).FaceAlpha = 0.8; a(3).FaceAlpha = 0.8;

xticks(x_vals); xticklabels(x_labels);
ylabel('Annual Procured Energy (GWh)');
xlabel({'Grid Outage Probability', '(b) Energy Supply Structure'}, ...
    'FontName', fontN, 'FontSize', 12, 'FontWeight', 'normal');
legend({'From DGC (Diesel)', 'From GEC (Green)', 'From Main Grid'}, 'Location', 'eastoutside');
grid on; box on; xlim([1, 5]); ylim([0, 14]);

%% ====================  (c) ====================
subplot(1, 3, 3);
plot(x_vals, pay_to_dgc/1e6, 'k-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor','w'); hold on;
plot(x_vals, pay_to_gec/1e6, 'g-d', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor','w');

fill([x_vals, fliplr(x_vals)], [nash_surplus/1e6, zeros(1,5)], [0.9 0.8 0.2], 'FaceAlpha', 0.4, 'EdgeColor', 'none');

xticks(x_vals); xticklabels(x_labels);
ylabel('Payment $)');
xlabel({'Grid Outage Probability', '(c) Nash Bargaining Financial Transfers'}, ...
    'FontName', fontN, 'FontSize', 12, 'FontWeight', 'normal');
legend({'Payment to DGC (C^p_{d})', 'Payment to GEC (C^p_{g})', 'Nash Surplus'}, 'Location', 'northwest');
grid on; box on; xlim([0.8, 5.2]);