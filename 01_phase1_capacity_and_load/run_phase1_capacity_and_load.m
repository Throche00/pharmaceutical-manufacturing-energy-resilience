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
fprintf(' Phase One: Throughput Capacity Planning and Load\n');
fprintf('======================================================\n\n');

%% --------------------lambda*  delta* --------------------
K = 6; 
S = 12; 

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

alpha_sk = unifrnd(0.1, 0.3, [S, K]);

we = unifrnd(50, 100, [S, 1]); 
ie = unifrnd(5, 10, [S, 1]);  

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

fprintf('System Limit\n');
fprintf('Tasks K=%d, Machines S=%d\n', K, S);
fprintf('lambda* = %.6f\n\n', lambda_star);
disp('delta* (S x K):'); 
disp(delta_opt);

season_names = {'Spring','Summer','Autumn','Winter'};

lambda_season = lambda_star * [0.60, 0.96, 0.69, 0.85]; 

T = 24; t_axis = 1:T;
Pmach = zeros(4,T); Pe = zeros(4,T); Ptot = zeros(4,T);
base_env = [380, 650, 460, 570]; 

fprintf('Typical daily electricity demand in four seasons\n');

for sidx = 1:4
    T_s = zeros(S, 1);
    for s = 1:S
        for k = 1:K
            T_s(s) = T_s(s) + delta_opt(s,k) * mu(k) * (1/v(s,k) + alpha_sk(s,k));
        end
    end
    
    utilization = lambda_season(sidx) * T_s; 
    
    P_mach_expected = sum(utilization .* we + (1 - utilization) .* ie);
    P_mach_base_t = zeros(1, T);
    for i = 1:T
        if i >= 8 && i <= 18
            P_mach_base_t(i) = P_mach_expected * 1.08; 
        else
            P_mach_base_t(i) = P_mach_expected * 0.92;
        end
    end
    
    P_e_base_t = zeros(1,T);
    for i=1:T
        if i >= 10 && i <= 16
            P_e_base_t(i) = base_env(sidx) * 1.35; 
        elseif i >= 8 && i <= 18
            P_e_base_t(i) = base_env(sidx) * 1.15;
        else
            P_e_base_t(i) = base_env(sidx) * 0.85;
        end
    end
    noise_mach = randn(1, T) .* (0.08 * P_mach_base_t); 
    noise_env  = randn(1, T) .* (0.11 * P_e_base_t);  
    
    Pmach(sidx,:) = max(P_mach_base_t + noise_mach, 0);
    Pe(sidx,:)    = max(P_e_base_t + noise_env, 0);
end

Ptot = Pmach + Pe;

%% -------------------- plots --------------------
colors = [0.30 0.65 0.35; 0.90 0.25 0.15; 0.80 0.45 0.20; 0.25 0.50 0.80];
mk = {'o','^','d','s'};
set(groot,'defaultAxesFontName','Times New Roman','defaultAxesFontSize',11);
figure('Name','Combined Load Profiles','Color','w','Position',[50 200 1600 400]);

% --- (a) ---
subplot(1,3,1);
hold on; grid on; box on;
for sidx=1:4
    plot(t_axis, Pe(sidx,:), '-','Color',colors(sidx,:),'LineWidth',1.8,...
        'Marker',mk{sidx},'MarkerSize',4.5,'DisplayName',season_names{sidx});
end
ylabel('Power (kW)');
xlabel({'Time (h)', '(a) Stochastic Environmental Demand ($P^e_t$)'}, 'Interpreter', 'latex', 'FontSize', 12);
legend('Location','best'); xlim([1 24]);

% --- (b) ---
subplot(1,3,2);
hold on; grid on; box on;
for sidx=1:4
    plot(t_axis, Pmach(sidx,:), '-','Color',colors(sidx,:),'LineWidth',1.8,...
        'Marker',mk{sidx},'MarkerSize',4.5,'DisplayName',season_names{sidx});
end
ylabel('Power (kW)');
xlabel({'Time (h)', '(b) Stochastic Production Demand ($P^p_{\bar{c},t}$)'}, 'Interpreter', 'latex', 'FontSize', 12);
legend('Location','best'); xlim([1 24]); ylim([min(Pmach(:))*0.8, max(Pmach(:))*1.2]);

% --- (c) ---
subplot(1,3,3);
hold on; grid on; box on;
for sidx=1:4
    plot(t_axis, Ptot(sidx,:), '-','Color',colors(sidx,:),'LineWidth',2.0,...
        'Marker',mk{sidx},'MarkerSize',5.0,'DisplayName',season_names{sidx});
end
ylabel('Total Power (kW)');
xlabel({'Time (h)', '(c) Total Electricity Demand'}, 'Interpreter', 'latex', 'FontSize', 12);
legend('Location','best'); xlim([1 24]);

save('Phase1_Stochastic_Loads.mat','lambda_star','delta_opt','Pmach','Pe','Ptot','alpha_sk','we','ie');
disp('saved to Phase1_Stochastic_Loads.mat');
