% This file generates the random perturbations to model parameters.
% It creates the matrix 'all_parameters' with the perturbations values,
% and the array 'parameter_names' with their definitions (strings).

% Uncomment line 59 for saving 'all_parameters' and 'parameter_names' into
% the mat file 'SA_par_matrix_1000_s0p1'.

close all
clear all
clc

%% Definitions of Parameters
% 1) GNa
% 2) GNaB
% 3) IbarNaK
% 4) Gtof
% 5) GKr
% 6) GKs
% 7) GKur
% 8) GKp
% 9) GK1
% 10) GK,ACh
% 11) GClCa
% 12) GClB
% 13) GCa (PCa, PK, PNa)
% 14) GCaB
% 15) IbarPMCA
% 16) IbarNCX
% 17) VmaxSERCA
% 18) RyR
% 19) SR_leak

parameter_names = {'GNa','GNaB','vNKA','Gtof','GKr','GKs',...
    'GKur','GKp','GK1','GKACh','GClCa','GClB',...
    'GCa','GCaB','vPMCA','vNCX','vSERCA','vRyR',...
    'vSRleak'} ;

n_parameters = length(parameter_names);
baseline_parameters = ones(1,n_parameters);

%% Random variations
variations = 1000; % number of trials

sigmaG = 0.1*ones(1,n_parameters); % standard deviation for parameters
% in this example, sigma is the same for all parameters

all_parameters = zeros(variations,n_parameters);
for ii = 1:n_parameters,
    scaling = exp(sigmaG(ii)*randn(1,variations)) ;
    newparams = baseline_parameters(ii)*scaling ;
    all_parameters(:,ii) = newparams ;
end

all_parameters %size(all_parameters)
% columns: N parameters
% rows: N trials

%% Saving
%save SA_par_matrix_1000_s0p1 all_parameters parameter_names