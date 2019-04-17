% This file runs the EAD protocol for all the models in the population.
% For each simulation, this file analyzes the output to assess the presence
% of EADs or DADs, and then saves the result in the matrix 'all_outputs'.

% Uncomment line 96 for saving 'all_outputs' into the mat file
% 'SA_outputs_matrix_1000_s0p1'.

clear all
close all
clc

%% Initial conditions
% load ICs (obtained with the baseline model with [ACh] = 0.1 uM)
load yf_ham_ina_ran_ACh0p1_1Hz
y0 = yfinal;

%% Parameters
% load matrix all_parameters (columns: N parameters, rows: N trials)
load SA_par_matrix_1000_s0p1 % sigma 0.1
[N_trials N_par] = size(all_parameters);

%% Input parameters
prot_index = 2; % 2) 'pace_cc_ead';

prot_rate = 1; % (Hz)
prot_interval = 500; % (ms
prot_vm = yfinal(39); % (mV)

% Ranolazine parameters
drug_index = 0; drug_conc = 0; % Drug Free

% Other experimental conditions
exp_Temp = 310; % [K]
exp_Nao = 140; % [Na]o
exp_ISO = 1; % (boolean)
exp_Ach = 1; % (boolean, if 1 [ACh] = 0.1 uM)

% Parameter array for passing nondefault conditions
prot_par = [prot_index prot_rate prot_interval prot_vm];    % 1 2 3 4
drug_par = [drug_index drug_conc];                          % 5 6
exp_par = [exp_Temp exp_Nao exp_ISO exp_Ach];               % 7 8 9 10
p = [prot_par drug_par exp_par]; 

% Sensitivity analysis parameters
%p_SA = ones(1,19);

duration = 25.510e3;
tspan = [0 duration];
options = odeset('RelTol',1e-5,'MaxStep',1,'Stats','on'); 

%% Run cycle

% Output definition
tin = [21e3 22e3 23e3 24e3 25e3]; % assessment of EAD presence in 5 beats
N_outputs = length(tin)+3; % additional 3 outputs are for DAD properties
all_outputs = zeros(N_trials,N_outputs);

tic
parfor ii=1:N_trials,
    X = sprintf('Run %d on %d',ii,N_trials); disp(X)
    p_SA = all_parameters(ii,:); % 19 parameters
    [t,y] = ode15s(@morotti_et_al_ham_ina_ran_model_SA,tspan,y0,options,p,p_SA);
    
    time = t;
    Vm = y(:,39);
    
    EAD_index = zeros(1,length(tin));
    for i=1:length(tin),
        EAD_index(i) = EAD_occurrence(time,Vm,tin(i));
    end
    
    Ca = y(:,38);
    tin_roi=find(time>19.9e3); tin_idx=tin_roi(1)-1;
    tfin_roi=find(time>21.1e3); tfin_idx=tfin_roi(1);
    Ca_max1 = max(Ca(tin_idx:tfin_idx));
    tin_roi=find(time>20.1e3); tin_idx=tin_roi(1)-1;
    tfin_roi=find(time>21e3); tfin_idx=tfin_roi(1);
    Ca_max2 = max(Ca(tin_idx:tfin_idx));
    if Ca_max2 < Ca_max1,
        DAD_index = 0;
    else
        DAD_index = 1;
    end
    Vmax2 = max(Vm(tin_idx:tfin_idx));
    Vmin2 = min(Vm(tin_idx:tfin_idx));
    
    all_outputs(ii,:) = [EAD_index, DAD_index, Vmax2, Vmin2];
end
toc

all_outputs
% columns: N outputs
% rows: N trials

%% Saving
%save SA_outputs_matrix_1000_s0p1 all_outputs