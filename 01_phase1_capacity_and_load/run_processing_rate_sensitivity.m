%% =========================================================================
if exist('RUN_FROM_MASTER', 'var') && RUN_FROM_MASTER
    clc; close all;
else
    clc; clear; close all;
end
set(groot,'defaultAxesFontName','Times New Roman','defaultAxesFontSize',11);

if exist('WS_SEED', 'var')
    base_seed = WS_SEED + 202;
else
    base_seed = 42;
end
rng(base_seed);

%% 1. Load Pre-calculated Baseline Data
if isfile('Phase1_Stochastic_Loads.mat')
    load('Phase1_Stochastic_Loads.mat', 'lambda_star');
    fprintf('Successfully loaded baseline lambda_star = %.4f\n', lambda_star);
else
    error('File Phase1_Stochastic_Loads.mat not found! Please run PH1.m first.');
end
fixed_arrival_rate = lambda_star; 

%% 2. Re-initialize Physical Machine Parameters
K = 6; S = 12;  

mu_base = [5.0, 5.0, 8.0, 8.0, 3.0, 3.0]; 
v_base = [
    8.5, 7.2, 3.1, 4.5, 2.5, 3.0;  8.0, 7.5, 3.5, 4.0, 2.8, 3.2;  
    4.0, 3.5, 8.8, 7.6, 4.5, 3.8;  4.5, 3.8, 8.2, 7.1, 4.0, 3.5;  
    3.0, 4.0, 4.5, 5.0, 8.5, 7.0;  3.5, 4.2, 4.8, 5.5, 8.0, 6.8;  
    6.0, 6.0, 6.0, 6.0, 6.0, 6.0;  5.5, 5.5, 5.5, 5.5, 5.5, 5.5;  
    7.0, 4.0, 7.0, 4.0, 7.0, 4.0;  4.0, 7.0, 4.0, 7.0, 4.0, 7.0;  
    9.5, 2.0, 9.5, 2.0, 9.5, 2.0;  2.0, 9.5, 2.0, 9.5, 2.0, 9.5   
];
alpha_sk = unifrnd(0.1, 0.3, [S, K]);
we = unifrnd(50, 100, [S, 1]); ie = unifrnd(5, 10, [S, 1]);  

% LP Setup
n_delta = S*K; n_var = n_delta + 1; idxLam = n_var;
f = zeros(n_var,1); f(idxLam) = 1;
Aeq = zeros(K,n_var); beq = ones(K,1);
for k = 1:K
    for s = 1:S
        idx = (s-1)*K + k; Aeq(k,idx) = 1;
    end
end
lb = zeros(n_var,1); ub = [ones(n_delta,1); inf];
opts = optimoptions('linprog','Display','none','Algorithm','dual-simplex');

%% 3. Experiment Settings
scale_ratios = 0.5 : 0.05 : 1.5;
N_tests = length(scale_ratios);
num_trials = 100; 
T_sim = 168; 
rec_lam_star  = zeros(1, N_tests);
energy_mean   = zeros(1, N_tests);
energy_ci     = zeros(1, N_tests);
delay_mean    = zeros(1, N_tests);
delay_ci      = zeros(1, N_tests);

%% 4. Execute Monte Carlo Simulation Loop
for i = 1:N_tests
    rng(base_seed + 1000); 
    
    current_v = v_base * scale_ratios(i); 
    
    A_test = zeros(S,n_var);
    for s = 1:S
        for k = 1:K
            idx = (s-1)*K + k;
            A_test(s,idx) = mu_base(k)/current_v(s,k) + alpha_sk(s,k)*mu_base(k);
        end
        A_test(s,idxLam) = -1;
    end
    [x, fval] = linprog(f,A_test,zeros(S,1),Aeq,beq,lb,ub,opts);
    rec_lam_star(i) = 1/fval;
    delta_opt = reshape(x(1:n_delta),K,S)';
    
    temp_energy = zeros(1, num_trials);
    temp_delay  = zeros(1, num_trials);
    
    for trial = 1:num_trials
        AT_current = 0; arrival_times = [];
        while AT_current < T_sim
            AT_current = AT_current + exprnd(1/fixed_arrival_rate);
            arrival_times = [arrival_times; AT_current];
        end
        N_orders = length(arrival_times);
        
        w_n = abs(normrnd(repmat(mu_base, N_orders, 1), 0.5));
        EP_prev = zeros(S, 1);       
        PT = zeros(S, T_sim);        
        Delay_list = zeros(N_orders, 1); 
        
        for n = 1:N_orders
            AT = arrival_times(n); 
            max_EP_this_order = AT; 
            
            for s = 1:S
                S_sn = sum(alpha_sk(s,:) .* w_n(n,:) .* delta_opt(s,:));
                P_sn = sum((w_n(n,:) .* delta_opt(s,:)) ./ current_v(s,:));
                
                if (S_sn + P_sn) > 0
                    SP_sn = max(AT, EP_prev(s)); 
                    EP_sn = SP_sn + S_sn + P_sn;
                    EP_prev(s) = EP_sn;
                    max_EP_this_order = max(max_EP_this_order, EP_sn);
                    
                    t_start_idx = floor(SP_sn) + 1;
                    t_end_idx   = ceil(EP_sn);
                    for t_hour = t_start_idx : min(t_end_idx, T_sim)
                        overlap_start = max(SP_sn, t_hour - 1);
                        overlap_end   = min(EP_sn, t_hour);
                        if overlap_start < overlap_end
                            PT(s, t_hour) = PT(s, t_hour) + (overlap_end - overlap_start);
                        end
                    end
                end
            end
            Delay_list(n) = max_EP_this_order - AT; 
        end
        
        PT = min(PT, 1.0); 
        
        P_mach_t = zeros(1, T_sim);
        for t_hour = 1:T_sim
            P_mach_t(t_hour) = sum(PT(:, t_hour) .* we + (1 - PT(:, t_hour)) .* ie);
        end
        
        temp_energy(trial) = mean(P_mach_t(48:end)) * 24; 
        temp_delay(trial) = mean(Delay_list);
    end
    
    energy_mean(i) = mean(temp_energy);
    energy_ci(i)   = 1.96 * std(temp_energy) / sqrt(num_trials);
    delay_mean(i)  = mean(temp_delay);
    delay_ci(i)    = 1.96 * std(temp_delay) / sqrt(num_trials);
    
    fprintf('Completed speed ratio %.2f (Capacity Limit: %.2f)\n', scale_ratios(i), rec_lam_star(i));
end

%% 5. Plot
figure('Position',[100 150 1200 450],'Color','w','Name','Machine Speed v_sk'); 

cross_idx = find(rec_lam_star > fixed_arrival_rate, 1);
if ~isempty(cross_idx)
    cross_x = scale_ratios(cross_idx) * 100;
end

x_axis = scale_ratios * 100;
fill_patch = @(x, y_mean, y_ci, color) fill([x, fliplr(x)], [y_mean + y_ci, fliplr(max(0, y_mean - y_ci))], color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% --- Plot (a) ---
subplot(1,3,1); hold on; grid on; box on;
plot(x_axis, rec_lam_star, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 5, 'DisplayName', 'Capacity Limit (\lambda^*)');
yline(fixed_arrival_rate, 'k-', 'Fixed Market Demand', 'LineWidth', 2, ...
    'LabelHorizontalAlignment', 'right', 'DisplayName', 'Market Demand');
if ~isempty(cross_idx)
    xline(cross_x, 'r--', 'Safety Threshold', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
end
ylabel('Throughput', 'FontWeight', 'normal');
xlabel({'Machine Speed v_{sk} Scaling (%)', '(a) Maximum Throughput'}, 'FontWeight', 'normal');
lgd_a = legend('Location', 'northwest', 'NumColumns', 2, 'FontSize', 8, 'Box', 'on');
lgd_a.ItemTokenSize = [10, 6];

% --- Plot (b) ---
subplot(1,3,2); hold on; grid on; box on;
fill_patch(x_axis, energy_mean, energy_ci, 'k'); 
plot(x_axis, energy_mean, 'k-s', 'LineWidth', 2.5, 'MarkerSize', 5, 'DisplayName', 'Mean Daily Energy');
if ~isempty(cross_idx)
    xline(cross_x, 'r--', 'Safety Threshold', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
end
ylabel('Energy Consumption (kW·h)', 'FontWeight', 'normal');
xlabel({'Machine Speed v_{sk} Scaling (%)', '(b) Average Electricity Consumption'}, 'FontWeight', 'normal');
legend('Location', 'northeast');
ylim([min(energy_mean-energy_ci)*0.95, max(energy_mean+energy_ci)*1.05]);

% --- Plot (c) ---
subplot(1,3,3); hold on; grid on; box on;
fill_patch(x_axis, delay_mean, delay_ci, 'm'); 
plot(x_axis, delay_mean, 'm-^', 'LineWidth', 2.5, 'MarkerSize', 5, 'MarkerFaceColor', 'm', 'DisplayName', 'Mean Order Lead Time');
if ~isempty(cross_idx)
    xline(cross_x, 'r--', 'Safety Threshold', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
end
ylabel('Order Delay Time (h)', 'FontWeight', 'normal');
xlabel({'Machine Speed v_{sk} Scaling (%)', '(c) Average Order Delay'}, 'FontWeight', 'normal');
legend('Location', 'northeast');
if ~isempty(cross_idx)
    ylim([0 max(delay_mean(x_axis >= cross_x - 10))*1.5]);
else
    ylim([0 max(delay_mean)*1.5]);
end

