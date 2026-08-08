%% 停电概率敏感性分析可视化
clc; clear; close all;

% X轴标签：
x_labels = {'2%, 8%', '4%, 12%', '6%, 18%', '8%, 22%', '10%, 25%'};
x_vals = 1:5;

% ==================== 1. 从 Excel 文件中提取的数据 ====================
% 1. 纳什红利 (Nash Surplus)
nash_surplus = [2410077.15, 3193936.69, 3991052.13, 5452434.87, 6806131.89];

% 2. 利润数据 (PM, GEC, DGC - 分别为非合作与合作)
pm_0 = [7595363.35, 6800332.16, 5941299.81, 4423760.09, 3007885.24];
pm_c = [8398846.55, 7865339.38, 7273019.965, 6244574.32, 5283263.06];
gec_0 = [3312014.33, 3312014.33, 3312014.33, 3312014.33, 3312014.33];
gec_c = [4114995.37, 4375623.74, 4640386.45, 5126608.38, 5577193.40];
dgc_0 = [8683577.22, 8683577.22, 8683577.221, 8683577.22, 8683577.22];
dgc_c = [9487168.97, 9748890.28, 10014491.56, 10500501.06, 10949117.93];

% 3. 缺电量 (Shortage - Baseline vs Cooperative)
shortage_0 = [573556.48, 857534.84, 1138113.777, 1652619.72, 2130434.69];
shortage_c = [0.01, 0.03, 0.03, 0.03, 0.02]; 

% 4. 合作模式下的制药厂能源构成
energy_grid_c = [7597277.44, 7409156.93, 7200919.17, 6822500.29, 6461574.97];
energy_gec_c = [2910878.16, 2939808.72, 3004389.09, 3114162.68, 3169004.24];
energy_dgc_c = [322388.20, 481569.27, 625250.76, 893931.80, 1199972.79];

% ==================== 2. 全套学术图表绘制 ====================
fontN = 'Times New Roman';
set(groot,'defaultAxesFontName',fontN,'defaultAxesFontSize',12);

% --- 数据预处理 ---
coalition_0 = pm_0 + gec_0 + dgc_0; % 非合作联盟总利润
coalition_c = pm_c + gec_c + dgc_c; % 合作联盟总利润
pay_to_gec = [2778956.08, 3068861.07, 3371304.531, 3936102.07, 4421021.74];
pay_to_dgc = [1095371.02, 1499130.32, 1893673.178, 2621492.39, 3345543.93];

figure('Name', 'Combined Outage Probability Sensitivity', 'Color', 'w', 'Position', [100, 100, 1600, 500]);

%% ==================== 子图 (a): 利润韧性与短缺惩罚对比 ====================
subplot(1, 3, 1);
% 左侧Y轴：利润
yyaxis left;
plot(x_vals, coalition_c/1e6, 'k-^', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor','k'); hold on;
plot(x_vals, coalition_0/1e6, 'k--', 'LineWidth', 1.5);
plot(x_vals, pm_c/1e6, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor','w');
plot(x_vals, pm_0/1e6, 'b--', 'LineWidth', 1.5);
ylabel('Profit (Million $)'); 
set(gca, 'ycolor', 'k'); ylim([0, 30]);

% 右侧Y轴：短缺电量
yyaxis right;
b = bar(x_vals, [shortage_0'/1e6, shortage_c'/1e6], 0.6, 'EdgeColor', 'k');
b(1).FaceColor = [0.8 0.4 0.4]; b(1).FaceAlpha = 0.7; % 红色
b(2).FaceColor = [0.2 0.6 0.8]; b(2).FaceAlpha = 0.7; % 蓝色
ylabel('Unmet Power Shortage (GWh)'); 
set(gca, 'ycolor', [0.8 0.3 0.3]); ylim([0, 2.5]);

% 坐标轴与图例设置
xticks(x_vals); xticklabels(x_labels);
xlabel({'Grid Outage Probability', '(a) System Resilience against Grid Outage'}, ...
    'FontName', fontN, 'FontSize', 12, 'FontWeight', 'normal');
legend({'Coalition Profit (Coop)', 'Coalition Profit (Non-Coop)', 'PM Profit (Coop)', 'PM Profit (Non-Coop)', ...
        'Shortage (Non-Coop)', 'Shortage (Coop)'}, 'Location', 'southwest', 'NumColumns', 2);
grid on; box on; xlim([0.5, 5.5]);

%% ==================== 子图 (b): 能源供应结构的动态演变 ====================
subplot(1, 3, 2);
energy_matrix = [energy_dgc_c', energy_gec_c', energy_grid_c'] / 1e6; % 转换为 GWh
a = area(x_vals, energy_matrix);
a(1).FaceColor = [0.5 0.5 0.5]; % 柴油 灰色
a(2).FaceColor = [0.3 0.7 0.4]; % 绿电 绿色
a(3).FaceColor = [0.2 0.5 0.8]; % 主网 蓝色
a(1).FaceAlpha = 0.8; a(2).FaceAlpha = 0.8; a(3).FaceAlpha = 0.8;

% 坐标轴与图例设置
xticks(x_vals); xticklabels(x_labels);
ylabel('Annual Procured Energy (GWh)');
xlabel({'Grid Outage Probability', '(b) Energy Supply Structure'}, ...
    'FontName', fontN, 'FontSize', 12, 'FontWeight', 'normal');
legend({'From DGC (Diesel)', 'From GEC (Green)', 'From Main Grid'}, 'Location', 'eastoutside');
grid on; box on; xlim([1, 5]); ylim([0, 14]);

%% ==================== 子图 (c): 风险溢价下的内部转移支付 ====================
subplot(1, 3, 3);
plot(x_vals, pay_to_dgc/1e6, 'k-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor','w'); hold on;
plot(x_vals, pay_to_gec/1e6, 'g-d', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor','w');

fill([x_vals, fliplr(x_vals)], [nash_surplus/1e6, zeros(1,5)], [0.9 0.8 0.2], 'FaceAlpha', 0.4, 'EdgeColor', 'none');

% 坐标轴与图例设置
xticks(x_vals); xticklabels(x_labels);
ylabel('Payment $)');
xlabel({'Grid Outage Probability', '(c) Nash Bargaining Financial Transfers'}, ...
    'FontName', fontN, 'FontSize', 12, 'FontWeight', 'normal');
legend({'Payment to DGC (C^p_{d})', 'Payment to GEC (C^p_{g})', 'Nash Surplus'}, 'Location', 'northwest');
grid on; box on; xlim([0.8, 5.2]);