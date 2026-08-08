%% =========================================================================
clc; close all;

doPlot = true;
if exist('WS_BATCH_MODE', 'var') && WS_BATCH_MODE
    doPlot = false;
end

fprintf('======================================================\n');
fprintf(' 蒙特卡洛模拟 (2000次) 与聚类\n');
fprintf('======================================================\n\n');

%% ==================== 1. 加载基准数据 ====================
try
    if exist('Phase1_Stochastic_Loads.mat', 'file')
        load('Phase1_Stochastic_Loads.mat', 'Ptot');
        fprintf('成功加载排产数据 (Phase1_Stochastic_Loads.mat)！\n');
    else
        load('Phase1_Styled_Seasonal_Loads.mat', 'Ptot');
        fprintf('成功加载基础排产数据 (Phase1_Styled_Seasonal_Loads.mat)！\n');
    end
    
    load('Phase2_GreenDiesel_Profiles.mat', 'Pg', 'Pd', 'seasons'); 
catch
    error('未找到 .mat 数据文件，请确保 Phase1 和 GEC/DGC 数据在当前目录。');
end

%% ==================== 2. 参数设定 ====================
n_scenarios = 2000; % 蒙特卡洛生成 2000 个场景
n_clusters = 3;     % 聚类为 3 个典型场景

% 波动率设定 (绿电 15%, 负荷 12%, 柴油 10%)
err_Pg = 0.15; err_Ptot = 0.12; err_Pd = 0.10;

% 预分配联合聚类结果存储变量
Expected_Pg   = zeros(4, n_clusters, 24);
Expected_Ptot = zeros(4, n_clusters, 24);
Expected_Pd   = zeros(4, n_clusters, 24);

% 联合概率（统一概率，保证物理时空一致性）
Prob_Joint = zeros(4, n_clusters);

seed_mc = 2026;
if exist('WS_SEED', 'var')
    seed_mc = WS_SEED + 30000;
end
rng(seed_mc, 'twister'); % 固定随机种子保证复现性

%% ==================== 3. 联合蒙特卡洛模拟与聚类 ====================
for s = 1:4
    % 获取该季节的 24h 基准曲线
    base_Pg   = Pg(s, :);
    base_Ptot = Ptot(s, :);
    base_Pd   = Pd(s, :);
    
    % 1) 绿电场景生成与期望对齐修正 (修复截断偏差)
    scen_Pg = repmat(base_Pg, n_scenarios, 1) + randn(n_scenarios, 24) .* (base_Pg * err_Pg);
    scen_Pg = max(0, scen_Pg);
    scen_Pg = scen_Pg .* (base_Pg ./ max(mean(scen_Pg, 1), 1e-6)); % 无损期望对齐
    
    % 2) 负荷场景生成与期望对齐修正
    scen_Ptot = repmat(base_Ptot, n_scenarios, 1) + randn(n_scenarios, 24) .* (base_Ptot * err_Ptot);
    scen_Ptot = max(0, scen_Ptot);
    scen_Ptot = scen_Ptot .* (base_Ptot ./ max(mean(scen_Ptot, 1), 1e-6));
    
    % 3) 柴油场景生成与期望对齐修正
    scen_Pd = repmat(base_Pd, n_scenarios, 1) + randn(n_scenarios, 24) .* (base_Pd * err_Pd);
    scen_Pd = max(0, scen_Pd);
    scen_Pd = scen_Pd .* (base_Pd ./ max(mean(scen_Pd, 1), 1e-6));
    
    % 4) 核心突破：将三者拼接为 72 维向量进行【联合聚类】
    % 这样保证了同一个场景下，风光、负荷、柴油是发生在一个物理平行宇宙里的
    joint_scenarios = [scen_Pg, scen_Ptot, scen_Pd]; 
    
    [idx, joint_centers] = kmeans(joint_scenarios, n_clusters, 'Replicates', 15, 'MaxIter', 1500);
    
    % 计算联合发生概率
    probs = [sum(idx==1), sum(idx==2), sum(idx==3)] / n_scenarios;
    
    % 为了图例颜色的一致性，按 Ptot（负荷）的总大小降序排列场景 S1, S2, S3
    Ptot_centers = joint_centers(:, 25:48);
    [~, sort_idx] = sort(sum(Ptot_centers, 2), 'descend');
    
    joint_centers = joint_centers(sort_idx, :);
    probs = probs(sort_idx);
    
    % 5) 将聚类后的 72 维中心解耦回 3 个变量
    Expected_Pg(s,:,:)   = joint_centers(:, 1:24);
    Expected_Ptot(s,:,:) = joint_centers(:, 25:48);
    Expected_Pd(s,:,:)   = joint_centers(:, 49:72);
    Prob_Joint(s,:)      = probs;
end

%% ==================== 4. 数据导出 ====================
Prob_Pg   = Prob_Joint;
Prob_Ptot = Prob_Joint;
Prob_Pd   = Prob_Joint;

save('Phase2_Clustered_Scenarios.mat', 'Expected_Pg', 'Expected_Ptot', 'Expected_Pd', ...
     'Prob_Pg', 'Prob_Ptot', 'Prob_Pd', 'seasons');
fprintf('联合聚类完成！典型场景数据已保存至 Phase2_Clustered_Scenarios.mat\n\n');

%% ==================== 5. 3D 可视化渲染  ====================
if doPlot
    z_labels = {'Green-electricity Gen. (kW)', 'Electricity Demand (kW)', 'Diesel Supp. (kW)'};
    fig_titles = {'Typical Scenarios of Green-electricity Generation', ...
                  'Typical Scenarios of Electricity Demand', ...
                  'Typical Scenarios of Diesel Supplementation'};
    colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
    plot_data = {Expected_Pg, Expected_Ptot, Expected_Pd};
    sub_titles = {'(a) Spring', '(b) Summer', '(c) Autumn', '(d) Winter'};
    
    for var_idx = 1:3
        current_centers = plot_data{var_idx};
        
        figure('Name', fig_titles{var_idx}, 'Position', [100, 100, 950, 700], 'Color', 'w');
        for s = 1:4
            subplot(2, 2, s);
            hold on; grid on;
            
            for k = 1:n_clusters
                plot3(1:24, ones(1,24)*k, squeeze(current_centers(s, k, :))', ...
                    'Color', colors(k,:), 'LineWidth', 2, ...
                    'DisplayName', sprintf('S%d (p=%.2f)', k, Prob_Joint(s,k)));
            end
            
            title('');
            xlabel('Time (h)', 'FontName', 'Times New Roman', 'FontSize', 12);
            ylabel('Typical Scenario', 'FontName', 'Times New Roman', 'FontSize', 12);
            zlabel(z_labels{var_idx}, 'FontName', 'Times New Roman', 'FontSize', 12);
            xticks([1 6 12 18 24]); yticks([1 2 3]);
            view(-35, 30); 
            
            ax = gca; ax.Box = 'on'; ax.LineWidth = 0.8; ax.GridColor = [0.6 0.6 0.6];
            legend('Location', 'northeast', 'FontSize', 8, 'FontName', 'Times New Roman');
            ax_pos = ax.Position;
            annotation(gcf, 'textbox', [ax_pos(1), max(ax_pos(2)-0.075, 0.01), ax_pos(3), 0.035], ...
                'String', sub_titles{s}, 'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontName', 'Times New Roman', 'FontSize', 12);
        end
    end
end
