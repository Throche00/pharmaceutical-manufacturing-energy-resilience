clc; clear; close all;
fontN = 'Times New Roman';

cap = [500, 700, 900, 1100, 1300]; 

coalition_coop = [21.003, 21.406, 21.825, 21.972, 22.124];
coalition_non = repmat(17.937, 1, 5);

shortage_coop = [29.306, 0.040, 0.0120, 0.007, 0.015];

grid_energy = [7.201, 7.201, 7.201, 7.201, 7.201];
green_energy = [3.175, 3.090, 3.004, 2.977, 2.962];
diesel_energy = [0.428, 0.540, 0.625, 0.653, 0.668];

pay_to_gec = [3.532, 3.461, 3.371, 3.342, 3.328];
pay_to_dgc = [1.676, 1.811, 1.894, 1.920, 1.935];

figure('Name', 'Combined Capacity Sensitivity', 'Color', 'w', 'Position', [100, 100, 1600, 500]);

% ==================== (a): ====================
subplot(1, 3, 1);
yyaxis left;
plot(cap, coalition_coop, 'k-^', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor','k'); hold on;
plot(cap, coalition_non, 'k--', 'LineWidth', 1.5);
plot(cap, coalition_coop-coalition_non, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor','w');
ylabel('Profit (Million $)', 'FontName', fontN, 'FontSize', 13);
set(gca, 'ycolor', 'k'); ylim([0, 30]);

yyaxis right;
b = bar(cap, shortage_coop, 0.4, 'FaceColor', [0.8 0.3 0.3], 'EdgeColor', 'k', 'FaceAlpha', 0.6);
ylabel('Unmet Power Shortage (MWh)', 'FontName', fontN, 'FontSize', 13);
set(gca, 'ycolor', [0.8 0.3 0.3]); ylim([0, 40]);

xlabel({'Diesel Contract Capacity P^p_{{d2p}}   (kW)', '(a) Impact of Contract Capacity on Profit and System Resilience'}, ...
    'FontName', fontN, 'FontSize', 13, 'FontWeight', 'normal');
legend('Coalition Profit (Coop)', 'Coalition Profit (Non-Coop)', 'PM Profit (Coop)', 'Unmet Shortage', 'Location', 'west', 'FontName', fontN);
grid on; box on; set(gca,'FontName',fontN,'FontSize',11);


% ==================== (b) ====================
subplot(1, 3, 2);
energy_matrix = [diesel_energy', green_energy', grid_energy'];
a = area(cap, energy_matrix);
a(1).FaceColor = [0.5 0.5 0.5]; 
a(2).FaceColor = [0.3 0.7 0.4]; 
a(3).FaceColor = [0.2 0.5 0.8];
a(1).FaceAlpha = 0.8; a(2).FaceAlpha = 0.8; a(3).FaceAlpha = 0.8;

ylabel('Annual Procured Energy (GWh)', 'FontName', fontN, 'FontSize', 13);
xlabel({'Diesel Contract Capacity P^p_{{d2p}}   (kW)', '(b) Energy Procurement Structure'}, ...
    'FontName', fontN, 'FontSize', 13, 'FontWeight', 'normal');
legend('From DGC (Diesel)', 'From GEC (Green)', 'From Main Grid', 'Location', 'eastoutside', 'FontName', fontN);
grid on; box on; set(gca,'FontName',fontN,'FontSize',11);
xlim([500, 1300]); ylim([0, 13]);


% ==================== (c) ====================
subplot(1, 3, 3);
plot(cap, pay_to_dgc, 'k-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor','w'); hold on;
plot(cap, pay_to_gec, 'g-d', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor','w');
plot(cap, pay_to_dgc + pay_to_gec, 'b-.', 'LineWidth', 2);

ylabel('Internal Payment (Million $)', 'FontName', fontN, 'FontSize', 13);
xlabel({'Diesel Contract Capacity P^p_{{d2p}}   (kW)', '(c) Nash Bargaining Financial Transfers'}, ...
    'FontName', fontN, 'FontSize', 13, 'FontWeight', 'normal');
legend('Payment to DGC (C^p_{d})', 'Payment to GEC (C^p_{g})', 'Total Energy Expenditure', 'Location', 'east', 'FontName', fontN);
grid on; box on; set(gca,'FontName',fontN,'FontSize',11);