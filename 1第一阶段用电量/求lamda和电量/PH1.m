%% =========================================================================
if ~exist('RUN_FROM_MASTER', 'var') || ~RUN_FROM_MASTER
    clc; clear; close all;
else
    close all;
end

seed_phase1 = 2026;
seed_env = getenv('WS_SEED');
if ~isempty(seed_env)
    seed_candidate = str2double(seed_env);
    if ~isnan(seed_candidate)
        seed_phase1 = seed_candidate;
    end
elseif exist('WS_SEED', 'var')
    seed_phase1 = WS_SEED;
end
rng(seed_phase1, 'twister');

fprintf('======================================================\n');
fprintf(' 启动 Phase One: 吞吐量容量规划与负荷生成\n');
fprintf('======================================================\n\n');

%% -------------------- [A] 计算全年最大吞吐量 lambda* 与分配矩阵 delta* --------------------
% 1. 生产网络设定（问题规模）
K = 6; % 任务类型: 6类药品工序
S = 12; % 机器数量: 12台异构设备

% 2. 核心参数生成方式
% T1:阿司匹林制粒, T2:阿司匹林压片, T3:二甲双胍制粒, T4:二甲双胍压片, T5:奥美拉唑制粒, T6:奥美拉唑包衣
mu = [5.0, 5.0, 8.0, 8.0, 3.0, 3.0]; 

v = [
    8.5, 7.2, 3.1, 4.5, 2.5, 3.0;  % M1
    8.0, 7.5, 3.5, 4.0, 2.8, 3.2;  % M2
    4.0, 3.5, 8.8, 7.6, 4.5, 3.8;  % M3
    4.5, 3.8, 8.2, 7.1, 4.0, 3.5;  % M4
    3.0, 4.0, 4.5, 5.0, 8.5, 7.0;  % M5
    3.5, 4.2, 4.8, 5.5, 8.0, 6.8;  % M6
    6.0, 6.0, 6.0, 6.0, 6.0, 6.0;  % M7  
    5.5, 5.5, 5.5, 5.5, 5.5, 5.5;  % M8  
    7.0, 4.0, 7.0, 4.0, 7.0, 4.0;  % M9  
    4.0, 7.0, 4.0, 7.0, 4.0, 7.0;  % M10 
    9.5, 2.0, 9.5, 2.0, 9.5, 2.0;  % M11 
    2.0, 9.5, 2.0, 9.5, 2.0, 9.5   % M12 
];
% -----------------------------------------------------------

% 任务部署/准备时间系数 \alpha_{sk} ~ U(0.1, 0.3)
alpha_sk = unifrnd(0.1, 0.3, [S, K]);

% 能耗参数 
% 工作状态电功率 we_s ~ U(50, 100) kW
we = unifrnd(50, 100, [S, 1]); 
% 待机空闲状态电功率 ie_s ~ U(5, 10) kW
ie = unifrnd(5, 10, [S, 1]);  

% 3. 构建并求解线性规划 (LP)
n_delta = S*K; n_var = n_delta + 1; idxLam = n_var;
f = zeros(n_var,1); f(idxLam) = 1;

A = zeros(S,n_var); b = zeros(S,1);
for s = 1:S
    for k = 1:K
        idx = (s-1)*K + k;
        A(s,idx) = mu(k)/v(s,k) + alpha_sk(s,k)*mu(k);
    end
    A(s,idxLam) = -1;
end

Aeq = zeros(K,n_var); beq = ones(K,1);
for k = 1:K
    for s = 1:S
        idx = (s-1)*K + k;
        Aeq(k,idx) = 1;
    end
end
lb = zeros(n_var,1); ub = [ones(n_delta,1); inf];

opts = optimoptions('linprog','Display','none','Algorithm','dual-simplex');
[x,fval,ef] = linprog(f,A,b,Aeq,beq,lb,ub,opts);
lambda_star = 1/fval;
delta_opt = reshape(x(1:n_delta),K,S)';

fprintf('【系统上限】\n');
fprintf('任务数量 K=%d, 机器数量 S=%d\n', K, S);
fprintf('全年理论最大吞吐量 lambda* = %.6f\n\n', lambda_star);
disp('最优任务分配矩阵 delta* (S x K):'); 
disp(delta_opt);

%% -------------------- [B] 改进版：依据公式映射生成四季典型日负荷 --------------------
season_names = {'Spring','Summer','Autumn','Winter'};

lambda_season = lambda_star * [0.60, 0.96, 0.69, 0.85]; 

T = 24; t_axis = 1:T;
Pmach = zeros(4,T); Pe = zeros(4,T); Ptot = zeros(4,T);
base_env = [380, 650, 460, 570]; 

fprintf('【生成四季典型日用电需求】\n');

for sidx = 1:4
    % 1. 计算单台机器处理一个订单的平均耗时 T_s
    T_s = zeros(S, 1);
    for s = 1:S
        for k = 1:K
            T_s(s) = T_s(s) + delta_opt(s,k) * mu(k) * (1/v(s,k) + alpha_sk(s,k));
        end
    end
    
    % 2. 依据到达率计算机器的平均利用率
    utilization = lambda_season(sidx) * T_s; 
    
    % 3. 生成基准机器负荷曲线
    P_mach_expected = sum(utilization .* we + (1 - utilization) .* ie);
    P_mach_base_t = zeros(1, T);
    for i = 1:T
        if i >= 8 && i <= 18
            P_mach_base_t(i) = P_mach_expected * 1.08; 
        else
            P_mach_base_t(i) = P_mach_expected * 0.92;
        end
    end
    
    % 4. 生成基准环境负荷曲线 (强化白天气温极值)
    P_e_base_t = zeros(1,T);
    for i=1:T
        if i >= 10 && i <= 16
            P_e_base_t(i) = base_env(sidx) * 1.35; % 正午及下午阳光直射，空调冷负荷极大
        elseif i >= 8 && i <= 18
            P_e_base_t(i) = base_env(sidx) * 1.15;
        else
            P_e_base_t(i) = base_env(sidx) * 0.85; % 夜间气温下降
        end
    end
    
    % 5. 加入正态分布扰动
    noise_mach = randn(1, T) .* (0.08 * P_mach_base_t); 
    noise_env  = randn(1, T) .* (0.11 * P_e_base_t);  
    
    % 加上扰动并防止出现负数
    Pmach(sidx,:) = max(P_mach_base_t + noise_mach, 0);
    Pe(sidx,:)    = max(P_e_base_t + noise_env, 0);
end

Ptot = Pmach + Pe;
fprintf('强季节性随机负荷生成完毕。\n\n');

%% -------------------- [C] 绘图展示 (组合图) --------------------
colors = [0.30 0.65 0.35; 0.90 0.25 0.15; 0.80 0.45 0.20; 0.25 0.50 0.80];
mk = {'o','^','d','s'};
set(groot,'defaultAxesFontName','Times New Roman','defaultAxesFontSize',11);

% 创建一个宽幅的 Figure 容纳三个子图
figure('Name','Combined Load Profiles','Color','w','Position',[50 200 1600 400]);

% --- (a) 环境用电 ---
subplot(1,3,1);
hold on; grid on; box on;
for sidx=1:4
    plot(t_axis, Pe(sidx,:), '-','Color',colors(sidx,:),'LineWidth',1.8,...
        'Marker',mk{sidx},'MarkerSize',4.5,'DisplayName',season_names{sidx});
end
ylabel('Power (kW)');
% 使用 xlabel 放置横坐标文本和子标题
xlabel({'Time (h)', '(a) Stochastic Environmental Demand ($P^e_t$)'}, 'Interpreter', 'latex', 'FontSize', 12);
legend('Location','best'); xlim([1 24]);

% --- (b) 机器用电 ---
subplot(1,3,2);
hold on; grid on; box on;
for sidx=1:4
    plot(t_axis, Pmach(sidx,:), '-','Color',colors(sidx,:),'LineWidth',1.8,...
        'Marker',mk{sidx},'MarkerSize',4.5,'DisplayName',season_names{sidx});
end
ylabel('Power (kW)');
xlabel({'Time (h)', '(b) Stochastic Production Demand ($P^p_{\bar{c},t}$)'}, 'Interpreter', 'latex', 'FontSize', 12);
legend('Location','best'); xlim([1 24]); ylim([min(Pmach(:))*0.8, max(Pmach(:))*1.2]);

% --- (c) 总用电 ---
subplot(1,3,3);
hold on; grid on; box on;
for sidx=1:4
    plot(t_axis, Ptot(sidx,:), '-','Color',colors(sidx,:),'LineWidth',2.0,...
        'Marker',mk{sidx},'MarkerSize',5.0,'DisplayName',season_names{sidx});
end
ylabel('Total Power (kW)');
xlabel({'Time (h)', '(c) Total Electricity Demand'}, 'Interpreter', 'latex', 'FontSize', 12);
legend('Location','best'); xlim([1 24]);

%% -------------------- [D] 保存数据 --------------------
save('Phase1_Stochastic_Loads.mat','lambda_star','delta_opt','Pmach','Pe','Ptot','alpha_sk','we','ie');
disp('数据已保存至 Phase1_Stochastic_Loads.mat，准备输入至 Phase Two。');
