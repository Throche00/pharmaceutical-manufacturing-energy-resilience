%% Sens_02_Monte.m 
if ~exist('RUN_FROM_MASTER', 'var') || ~RUN_FROM_MASTER
    clc; clear; close all;
else
    close all;
end

load('SensData_PH1.mat');

if ~exist('mu_list','var') && exist('mu_scale_list','var')
    mu_list = mu_scale_list;
end
if ~exist('sigma_list','var') && exist('sigma_scale_list','var')
    sigma_list = sigma_scale_list;
end

try
    load('Phase2_GreenDiesel_Profiles.mat', 'Pg', 'Pd'); 
catch
    error('No Phase2_GreenDiesel_Profiles.mat');
end

SensData_Monte = cell(Nm, Ns);
n_scenarios = 2000; 
n_clusters = 3;     
rng(2026, 'twister'); 

for i = 1:Nm
    for j = 1:Ns
        Ptot = SensData_PH1{i,j}.Ptot;
        Exp_Pg = zeros(4, n_clusters, 24);
        Exp_Ptot = zeros(4, n_clusters, 24);
        Exp_Pd = zeros(4, n_clusters, 24);
        Prob_Joint = zeros(4, n_clusters);
        
        for s = 1:4
            scen_Pg = max(0, repmat(Pg(s,:), n_scenarios, 1) + randn(n_scenarios, 24).*(0.15*Pg(s,:)));
            scen_Pg = scen_Pg .* (Pg(s,:) ./ max(mean(scen_Pg, 1), 1e-6));
            
            scen_Ptot = max(0, repmat(Ptot(s,:), n_scenarios, 1) + randn(n_scenarios, 24).*(0.12*Ptot(s,:)));
            scen_Ptot = scen_Ptot .* (Ptot(s,:) ./ max(mean(scen_Ptot, 1), 1e-6));
            
            scen_Pd = max(0, repmat(Pd(s,:), n_scenarios, 1) + randn(n_scenarios, 24).*(0.10*Pd(s,:)));
            scen_Pd = scen_Pd .* (Pd(s,:) ./ max(mean(scen_Pd, 1), 1e-6));
            
            joint_scenarios = [scen_Pg, scen_Ptot, scen_Pd]; 
            [idx_c, centers] = kmeans(joint_scenarios, n_clusters, 'MaxIter', 500, 'Replicates', 10);
            
            probs = [sum(idx_c==1), sum(idx_c==2), sum(idx_c==3)] / n_scenarios;
            Ptot_cents = centers(:, 25:48);
            [~, sort_idx] = sort(sum(Ptot_cents, 2), 'descend');
            centers = centers(sort_idx, :); 
            probs = probs(sort_idx);
            
            Exp_Pg(s,:,:) = centers(:, 1:24); 
            Exp_Ptot(s,:,:) = centers(:, 25:48);
            Exp_Pd(s,:,:) = centers(:, 49:72); 
            Prob_Joint(s,:) = probs;
        end
        
        res.Expected_Pg = Exp_Pg; 
        res.Expected_Ptot = Exp_Ptot;
        res.Expected_Pd = Exp_Pd; 
        res.Prob_Joint = Prob_Joint;
        SensData_Monte{i,j} = res;
    end
end

if exist('mu_scale_list','var') && exist('sigma_scale_list','var')
    save('SensData_Monte.mat', 'SensData_Monte', ...
         'mu_list', 'sigma_list', 'mu_scale_list', 'sigma_scale_list', 'Nm', 'Ns');
else
    save('SensData_Monte.mat', 'SensData_Monte', ...
         'mu_list', 'sigma_list', 'Nm', 'Ns');
end

disp('Phase 2 done');
