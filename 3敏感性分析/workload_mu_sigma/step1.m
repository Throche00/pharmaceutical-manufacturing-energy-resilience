%% Sens_01_PH1.m 
if ~exist('RUN_FROM_MASTER', 'var') || ~RUN_FROM_MASTER
    clc; clear; close all;
else
    close all;
end
fprintf('=== 启动 Phase 1: Workload 二维敏感性分析 (带物理上限) ===\n');

% 参数网格 (5x5)
mu_scale_list    = [0.6, 0.8, 1.0, 1.2, 1.4];   % 所有任务均值同步缩放
sigma_scale_list = [0.2, 0.6, 1.0, 1.4, 1.8];   % 所有任务标准差同步缩放
Nm = length(mu_scale_list); Ns = length(sigma_scale_list);

K = 6; S = 12; 

mu_base = [5.0, 5.0, 8.0, 8.0, 3.0, 3.0];
sigma_base = 0.5 * ones(1,K);   % 基准：所有任务同一标准差

v_base = [
    8.5, 7.2, 3.1, 4.5, 2.5, 3.0;
    8.0, 7.5, 3.5, 4.0, 2.8, 3.2;
    4.0, 3.5, 8.8, 7.6, 4.5, 3.8;
    4.5, 3.8, 8.2, 7.1, 4.0, 3.5;
    3.0, 4.0, 4.5, 5.0, 8.5, 7.0;
    3.5, 4.2, 4.8, 5.5, 8.0, 6.8;
    6.0, 6.0, 6.0, 6.0, 6.0, 6.0;
    5.5, 5.5, 5.5, 5.5, 5.5, 5.5;
    7.0, 4.0, 7.0, 4.0, 7.0, 4.0;
    4.0, 7.0, 4.0, 7.0, 4.0, 7.0;
    9.5, 2.0, 9.5, 2.0, 9.5, 2.0;
    2.0, 9.5, 2.0, 9.5, 2.0, 9.5
];
if exist('Phase1_Stochastic_Loads.mat', 'file')
    base_data = load('Phase1_Stochastic_Loads.mat');
    alpha_base = base_data.alpha_sk; % 将基准的 alpha_sk 赋给 step1 使用的 alpha_base
    we_base    = base_data.we;       % 将基准的 we 赋给 step1 使用的 we_base
    ie_base    = base_data.ie;       % 将基准的 ie 赋给 step1 使用的 ie_base
else
    error('未找到 Phase1_Stochastic_Loads.mat，请先运行更新后的 PH1.m');
end
base_env = [380, 650, 460, 570]; 

mu_baseline = mu_base;
n_delta = S*K; n_var = n_delta + 1; idxLam = n_var;
f = zeros(n_var,1); f(idxLam) = 1;
A = zeros(S,n_var); b = zeros(S,1);
for s = 1:S
    for k = 1:K, A(s,(s-1)*K+k) = mu_baseline(k)/v_base(s,k) + alpha_base(s,k)*mu_baseline(k); end
    A(s,idxLam) = -1;
end
Aeq = zeros(K,n_var); beq = ones(K,1);
for k = 1:K, for s = 1:S, Aeq(k,(s-1)*K+k) = 1; end; end
opts = optimoptions('linprog','Display','none');
[~,fval,~] = linprog(f,A,b,Aeq,beq,zeros(n_var,1),[ones(n_delta,1); inf],opts);
lambda_baseline = 1/fval;
fixed_lambda = lambda_baseline * [0.75, 0.96, 0.82, 0.89];
fprintf('基准设定完毕: mu~N(5, 0.5), 锁定夏季基准到达率 = %.4f\n\n', fixed_lambda(2));

SensData_PH1 = cell(Nm, Ns);
T = 24;

for i = 1:Nm
    for j = 1:Ns
        mu_k_mean  = mu_base * mu_scale_list(i);             % 全任务均值同比例变化
        mu_k_sigma = sigma_base * sigma_scale_list(j);       % 全任务标准差同比例变化
        
        if mu_scale_list(i) == 1.0 && sigma_scale_list(j) == 1.0
            mu_k = mu_base; 
        else
            rng(1000 * i + j, 'twister'); 
            mu_k = abs(normrnd(mu_k_mean, mu_k_sigma, [1, K]));
        end
        
        A = zeros(S,n_var);
        for s = 1:S
            for k = 1:K, A(s,(s-1)*K+k) = mu_k(k)/v_base(s,k) + alpha_base(s,k)*mu_k(k); end
            A(s,idxLam) = -1;
        end
        [x,fval,~] = linprog(f,A,b,Aeq,beq,zeros(n_var,1),[ones(n_delta,1); inf],opts);
        lambda_star = 1/fval;
        delta_opt = reshape(x(1:n_delta),K,S)';
        
        Pmach = zeros(4,T); Pe = zeros(4,T);
        for sidx = 1:4
            T_s = zeros(S, 1);
            for s = 1:S
                for k = 1:K, T_s(s) = T_s(s) + delta_opt(s,k)*mu_k(k)*(1/v_base(s,k)+alpha_base(s,k)); end
            end
            
            % 【核心修正：物理上限截断！】利用率最高100%，机器24小时无休
            utilization = min(1.0, fixed_lambda(sidx) * T_s); 
            
            P_mach_expected = sum(utilization .* we_base + (1 - utilization) .* ie_base);
            P_mach_base_t = repmat(P_mach_expected, 1, T); 
            
            P_e_base_t = zeros(1,T);
            for h=1:T
                 if h >= 10 && h <= 16
                    P_e_base_t(h) = base_env(sidx) * 1.35;
                  elseif h >= 8 && h <= 18
                    P_e_base_t(h) = base_env(sidx) * 1.15;
                  else
                    P_e_base_t(h) = base_env(sidx) * 0.85; 
                  end
            end
            
            rng(2026 + i + j, 'twister');
            noise_mach = randn(1, T) .* (0.08 * P_mach_base_t);
            noise_env  = randn(1, T) .* (0.11 * P_e_base_t); 
            Pmach(sidx,:) = max(P_mach_base_t + noise_mach, 0);
            Pe(sidx,:)    = max(P_e_base_t + noise_env, 0);
        end
        
        res.lambda_star = lambda_star;
        res.util_max = max(min(1.0, fixed_lambda(2) * T_s)); 
        res.Ptot = Pmach + Pe;
        SensData_PH1{i, j} = res;
    end
end
save('SensData_PH1.mat', 'SensData_PH1', 'mu_scale_list', 'sigma_scale_list', 'Nm', 'Ns', ...
     'mu_base', 'sigma_base', 'fixed_lambda');
disp('Phase 1 数据生成完毕。');
