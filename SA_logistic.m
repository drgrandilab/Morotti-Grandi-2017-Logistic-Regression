% This file performs the regression analysis and plots the results.
% It loads the matrix with the perturbations of model parameters
% (in 'SA_par_matrix_1000_s0p1') and the matrix with the simulation results
% (in 'SA_outputs_matrix_1000_s0p1').

clear all
close all
clc

color = [0 0 0];

%% Load parameters
load SA_par_matrix_1000_s0p1
[N_trials N_pars] = size(all_parameters);

%% Load outputs
load SA_outputs_matrix_1000_s0p1 
% all_outputs: 1000 x 8 (1-5 EAD index, 6 DAD index, 7-8 Vmax & Vmin during pause)
all_outputs_dad = all_outputs(:,6);
all_outputs_ead = all_outputs(:,1:5);

aos = size(all_outputs_ead); % 1000 x 5
N_beat = aos(2);
N_cell = aos(1);

%% Analysis DAD/no DAD
dad_presence = all_outputs_dad; % 1 for DAD occurrence, 0 for no DAD
fraction_dad = sum(dad_presence)/N_cell;

EAD_in_DAD_cells = all_outputs(find(dad_presence>0.5),1:5);
Vmax_in_DAD_cells = all_outputs(find(dad_presence>0.5),7);

%% Analysis EAD/no EAD
all_outputs_ead_sum = sum(all_outputs_ead'); % 1 x 1000
ead_presence = (all_outputs_ead_sum>1/2); % 1 for EAD occurrence, 0 for no EAD
disp('Fraction of simulations with at least 1 EAD:');
fraction_ead = sum(ead_presence)/N_cell

%% Logistic regression - EAD/no EAD
parameter_names_b0 = {'b0',...
    'GNa','GNaB','vNKA','Gtof','GKr','GKs',...
    'GKur','GKp','GK1','GKACh','GClCa','GClB',...
    'GCa','GCaB','vPMCA','vNCX','vSERCA','vRyR',...
    'vSRleak'} ;

allpars_LOGISTIC = all_parameters;
X_LOG = log(allpars_LOGISTIC) ;
for ii=1:N_pars, % z-score
    X_LOGISTIC(:,ii)=(X_LOG(:,ii)-mean(X_LOG(:,ii)))/std(X_LOG(:,ii));
end

Y_LOGISTIC = 1-(ead_presence-1); % positive integer!
% ead_presence: 0 with no EAD, 1 with EADs 
% Y_LOGISTIC: 1 with EADs, 2 with no EADs % positive integer!

[B_LOGISTIC,dev,stats] = mnrfit(X_LOGISTIC,Y_LOGISTIC);

figure; set(gcf,'color','w') % with b0
bar(B_LOGISTIC,'FaceColor',color)
set(gca,'box','off','tickdir','out','fontsize',10)
title('Probability EAD development')
set(gca,'XTick',1:N_pars+1)
set(gca,'XTickLabel',parameter_names_b0)
set(gca,'XLim',[0 N_pars+1+1])
rotateXLabels( gca(), 90)

%% EAD Probability
% PEAD in the baseline model
B0 = B_LOGISTIC(1);
pval_LOGISTIC = stats.p;
disp('PEAD in the baseline model:');
P_ead_mean_B0 = 1/(1+exp(-(B0)))

% PEAD in each model of the population
P_ead_array = zeros(1,N_trials);
for iii=1:N_trials,
    P_ead_array(iii) = 1/(1+exp(-(B0+sum(B_LOGISTIC(2:end).*X_LOGISTIC(iii,:)'))));
end

figure,set(gcf,'color','w')
plot((1:N_trials),P_ead_array,'*','Color',color)
set(gca,'box','off','tickdir','out','fontsize',10)
title('Probability EAD development')
xlabel('Trial')
ylabel('Probability (-)')
P_ead_mean = mean(P_ead_array);
P_ead_std = std(P_ead_array);

%% Plot PEAD as function of modulation in model parameters (Fig. 4B)
rangeM = (0.5:0.01:1.5); % modulation range

% GNa - index in B 2
indexGNa = 2;
muGNa = mean(X_LOG(:,indexGNa-1));
sigmaGNa = std(X_LOG(:,indexGNa-1));
rangeGNa = (log(rangeM)-muGNa)/sigmaGNa;
PGNa = 1./(1+exp(-(B0 + B_LOGISTIC(indexGNa)*rangeGNa)));
% GK,ACh - index in B 11
indexGKACh = 11;
muGKACh = mean(X_LOG(:,indexGKACh-1));
sigmaGKACh = std(X_LOG(:,indexGKACh-1));
rangeGKACh = (log(rangeM)-muGKACh)/sigmaGKACh;
PGKACh = 1./(1+exp(-(B0 + B_LOGISTIC(indexGKACh)*rangeGKACh))); 
% vNCX - index in B 17
indexVncx = 17;
muVncx = mean(X_LOG(:,indexVncx-1));
sigmaVncx = std(X_LOG(:,indexVncx-1));
rangeVncx = (log(rangeM)-muVncx)/sigmaVncx;
PVncx = 1./(1+exp(-(B0 + B_LOGISTIC(indexVncx)*rangeVncx))); 
% GKur - index in B 8
indexGKur = 8;
muGKur = mean(X_LOG(:,indexGKur-1));
sigmaGKur = std(X_LOG(:,indexGKur-1));
rangeGKur = (log(rangeM)-muGKur)/sigmaGKur;
PGKur = 1./(1+exp(-(B0 + B_LOGISTIC(indexGKur)*rangeGKur)));
% GKs - index in B 7
indexGKs = 7;
muGKs = mean(X_LOG(:,indexGKs-1));
sigmaGKs = std(X_LOG(:,indexGKs-1));
rangeGKs = (log(rangeM)-muGKs)/sigmaGKs;
PGKs = 1./(1+exp(-(B0 + B_LOGISTIC(indexGKs)*rangeGKs)));

% Figure
figure, set(gcf,'color','w'), hold on, grid on
plot(rangeM,PGNa,rangeM,PGKACh,rangeM,PVncx,rangeM,PGKur,rangeM,PGKs)
ylabel('PEAD (-)'), xlabel('Scale Factor (-)')
title('Probability EAD development')
legend('GNa','GK,ACh','vNCX','GKur','GKs')
set(gca,'box','off','tickdir','out','fontsize',12)

%% Boxplot (Fig. 4C)

% Identification of the 2 sub-populations EAD/no EAD
groupEAD = find(ead_presence>0.5);
groupNoEAD = find(ead_presence<0.5);

disp('Mean perturbation values in the 2 sub-populations EAD/no EAD:');

% Analysis GNa
GNa_perturbations = all_parameters(:,1);
GNa_groupEAD = GNa_perturbations(groupEAD);
GNa_groupNoEAD = GNa_perturbations(groupNoEAD);
GNa_groupEAD_mean = mean(GNa_groupEAD)
GNa_groupEAD_std = std(GNa_groupEAD);
GNa_groupNoEAD_mean = mean(GNa_groupNoEAD)
GNa_groupNoEAD_std = std(GNa_groupNoEAD);

% Analysis GK,ACh in the 2 sub-populations EAD/no EAD
GKACh_perturbations = all_parameters(:,10);
GKACh_groupEAD = GKACh_perturbations(groupEAD);
GKACh_groupNoEAD = GKACh_perturbations(groupNoEAD);
GKACh_groupEAD_mean = mean(GKACh_groupEAD)
GKACh_groupEAD_std = std(GKACh_groupEAD);
GKACh_groupNoEAD_mean = mean(GKACh_groupNoEAD)
GKACh_groupNoEAD_std = std(GKACh_groupNoEAD);

% Analysis vNCX in the 2 sub-populations EAD/no EAD
vNCX_perturbations = all_parameters(:,16);
vNCX_groupEAD = vNCX_perturbations(groupEAD);
vNCX_groupNoEAD = vNCX_perturbations(groupNoEAD);
vNCX_groupEAD_mean = mean(vNCX_groupEAD)
vNCX_groupEAD_std = std(vNCX_groupEAD);
vNCX_groupNoEAD_mean = mean(vNCX_groupNoEAD)
vNCX_groupNoEAD_std = std(vNCX_groupNoEAD);

% Analysis GKur in the 2 sub-populations EAD/no EAD
GKur_perturbations = all_parameters(:,7);
GKur_groupEAD = GKur_perturbations(groupEAD);
GKur_groupNoEAD = GKur_perturbations(groupNoEAD);
GKur_groupEAD_mean = mean(GKur_groupEAD)
GKur_groupEAD_std = std(GKur_groupEAD);
GKur_groupNoEAD_mean = mean(GKur_groupNoEAD)
GKur_groupNoEAD_std = std(GKur_groupNoEAD);

% Analysis GKs in the 2 sub-populations EAD/no EAD
GKs_perturbations = all_parameters(:,6);
GKs_groupEAD = GKs_perturbations(groupEAD);
GKs_groupNoEAD = GKs_perturbations(groupNoEAD);
GKs_groupEAD_mean = mean(GKs_groupEAD)
GKs_groupEAD_std = std(GKs_groupEAD);
GKs_groupNoEAD_mean = mean(GKs_groupNoEAD)
GKs_groupNoEAD_std = std(GKs_groupNoEAD);

% Figures
figure,set(gcf,'color','w')
set(gcf, 'Position', [15, 100, 1650, 400]);
subplot(1,6,1),boxplot(GNa_groupEAD,'Notch','on','Labels','EAD','Symbol','ro')
ylim([0.7 1.45]),ylabel('GNa perturbations (-)')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(1,6,2),boxplot(GNa_groupNoEAD,'Notch','on','Labels','No EAD','Symbol','ro')
ylim([0.7 1.45]),ylabel('GNa perturbations (-)')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(1,6,3),boxplot(GKACh_groupEAD,'Notch','on','Labels','EAD','Symbol','ro')
ylim([0.7 1.45]),ylabel('GK,ACh perturbations (-)')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(1,6,4),boxplot(GKACh_groupNoEAD,'Notch','on','Labels','No EAD','Symbol','ro')
ylim([0.7 1.45]),ylabel('GK,ACh perturbations (-)')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(1,6,5),boxplot(vNCX_groupEAD,'Notch','on','Labels','EAD','Symbol','ro')
ylim([0.7 1.45]),ylabel('vNCX perturbations (-)')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(1,6,6),boxplot(vNCX_groupNoEAD,'Notch','on','Labels','No EAD','Symbol','ro')
ylim([0.7 1.45]),ylabel('vNCX perturbations (-)')
set(gca,'box','off','tickdir','out','fontsize',12)

figure,set(gcf,'color','w')
set(gcf, 'Position', [15, 100, 1650, 400]);
subplot(1,6,1),boxplot(GKur_groupEAD,'Notch','on','Labels','EAD','Symbol','ro')
ylim([0.7 1.45]),ylabel('GKur perturbations (-)')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(1,6,2),boxplot(GKur_groupNoEAD,'Notch','on','Labels','No EAD','Symbol','ro')
ylim([0.7 1.45]),ylabel('GKur perturbations (-)')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(1,6,3),boxplot(GKs_groupEAD,'Notch','on','Labels','EAD','Symbol','ro')
ylim([0.7 1.45]),ylabel('GKs perturbations (-)')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(1,6,4),boxplot(GKs_groupNoEAD,'Notch','on','Labels','No EAD','Symbol','ro')
ylim([0.7 1.45]),ylabel('GKs perturbations (-)')
set(gca,'box','off','tickdir','out','fontsize',12)

%% Plot X, Y, B matrices (Fig. 1)
% X - z-score
figure; set(gcf,'color','w')
imagesc(X_LOGISTIC); colormap jet;
set(gca,'box','off','tickdir','out','fontsize',10)
set(gca,'YDir','normal')
title('Parameters (X) - z-score');
xlabel('Parameters');
ylabel('Trials');
%set(gca,'YTick',(1:N_trials))
set(gca,'XTick',(1:N_pars))
set(gca,'XTickLabel',parameter_names)
rotateXLabels( gca(), 90)
colorbar

% Y (EAD)
figure; set(gcf,'color','w')
imagesc(ead_presence'); colormap jet;
set(gca,'box','off','tickdir','out','fontsize',10)
set(gca,'YDir','normal')
title('Outputs (Y)');
xlabel('Output (yes/no)');
ylabel('Trials');
%set(gca,'YTick',(1:N_trials))
set(gca,'XTick',1)
set(gca,'XTickLabel','EAD')
%rotateXLabels( gca(), 90)
colorbar

% B (EAD)
figure; set(gcf,'color','w')
imagesc(B_LOGISTIC); colormap jet;
set(gca,'box','off','tickdir','out','fontsize',10)
set(gca,'YDir','normal')
title('Regression coefficients (B)');
xlabel('Output (yes/no)');
ylabel('Coefficients');
set(gca,'YTick',(1:N_pars+1))
set(gca,'YTickLabel',parameter_names_b0)
set(gca,'XTick',1)
set(gca,'XTickLabel','EAD')
%rotateXLabels( gca(), 90)
colorbar