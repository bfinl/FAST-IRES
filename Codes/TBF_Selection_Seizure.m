% This code cleans up the seizure data from noisy artefacts and selects the
% components that are relevant to seizure (those that demonstrate 
% correlation with seizure onset time).
% data_noise Sohrabpour
% 1/9/2018

clear all; close all; clc

% This part of the code was to load the data for our data base. We have
% saved data here so you can directly load your raw data, i.e. data_sz_tot.

% pause(0.01)
% Loop_Num                                = 1;
% data_sz_tot                            = [];
% for i_loop = 1 : Loop_Num
%     LoadCurryDataFile_AS
%     data_sz_tot = [data_sz_tot, data];
% end

cd('../Raw 273')
Elec_loc           = readlocs('273_ele_ERD.xyz');
load('data_sz_tot.mat')
clearvars -except    data_sz_tot Elec_loc

cd('../Codes')
Samp_Rate            = 500;
data_noise           = data_sz_tot(:,1:end);
[A1, A2, A3]         = svd(data_noise,'econ');
figure; plot(diag(A2)); hold on; plot(ones(size(A2,1))*mean(diag(A2)),'r')
sum(diag(A2)>=mean(diag(A2)))

f = msgbox('Press OK to continue');
uiwait(f);close all

Num_ICA              = 25;
[weights, sphere,compvars,bias,signs,lrates,activations] ...
                   = runica(data_noise,'extended',-1, 'pca', Num_ICA, 'sphering','on','maxsteps',2048);
Topo_ICA           = pinv(weights*sphere);

fig_pos = [331 680 560 420; 893 680 560 420; 49 174 560 420; 611 174 560 420; 1173 174 560 420];
for i_fig=0:5:Num_ICA-5
    figure; set(gcf,'position',fig_pos(i_fig/5+1,:));
    subplot (5,2,1); topoplot(Topo_ICA(:,i_fig+1),Elec_loc,'plotrad',0.55,'electrodes','on'); subplot (5,2,2);plot(activations(i_fig+1,:)); %line([1+wind_ind_num 1+wind_ind_num],[-5 5],'color','red','linewidth',2,'LineStyle','--'); line([1+wind_ind_num+67 1+wind_ind_num+67],[-5 5],'color','red','linewidth',2,'LineStyle','--')
    subplot (5,2,3); topoplot(Topo_ICA(:,i_fig+2),Elec_loc,'plotrad',0.55,'electrodes','on'); subplot (5,2,4);plot(activations(i_fig+2,:)); %line([1+wind_ind_num 1+wind_ind_num],[-5 5],'color','red','linewidth',2,'LineStyle','--'); line([1+wind_ind_num+67 1+wind_ind_num+67],[-5 5],'color','red','linewidth',2,'LineStyle','--')
    subplot (5,2,5); topoplot(Topo_ICA(:,i_fig+3),Elec_loc,'plotrad',0.55,'electrodes','on'); subplot (5,2,6);plot(activations(i_fig+3,:)); %line([1+wind_ind_num 1+wind_ind_num],[-5 5],'color','red','linewidth',2,'LineStyle','--'); line([1+wind_ind_num+67 1+wind_ind_num+67],[-5 5],'color','red','linewidth',2,'LineStyle','--')
    subplot (5,2,7); topoplot(Topo_ICA(:,i_fig+4),Elec_loc,'plotrad',0.55,'electrodes','on'); subplot (5,2,8);plot(activations(i_fig+4,:)); %line([1+wind_ind_num 1+wind_ind_num],[-5 5],'color','red','linewidth',2,'LineStyle','--'); line([1+wind_ind_num+67 1+wind_ind_num+67],[-5 5],'color','red','linewidth',2,'LineStyle','--')
    subplot (5,2,9); topoplot(Topo_ICA(:,i_fig+5),Elec_loc,'plotrad',0.55,'electrodes','on'); subplot (5,2,10);plot(activations(i_fig+5,:));%line([1+wind_ind_num 1+wind_ind_num],[-5 5],'color','red','linewidth',2,'LineStyle','--'); line([1+wind_ind_num+67 1+wind_ind_num+67],[-5 5],'color','red','linewidth',2,'LineStyle','--')
end


% How much data variance is not explained

data_recon                  = Topo_ICA*weights*data_noise;
100*norm(data_recon - data_noise)/norm(data_noise)

% Select Components that are noisy, leave components are time locked with
% seizure onset and are not clearly noisy artefacts.

Noise_comp                  = Topo_ICA;
ind_noise                   = [1:2,5:9,12:Num_ICA];
Noise_comp (:,ind_noise)    = 0;
data_no_blink               = Noise_comp*weights*data_noise;

% Select the Non-noise components, those that are not ind_noise, and plot
% for inspection
% runica can change from run-2-run, as in our tested case, we selected the
% following components according to their morphology and activations
ind_select = [3 4 10 11];
% ask for ICA components
prompt = {'Enter component numbers (separate numbers with comma, no space): We selected components [3 4 10 11], but ICA might change from run to run and you might have to select other components to match our results; try 12 instead of 11 as well.'};
title = 'TBF Selection - Seizure Components';
dims = [1 50];
definput = {''};
answer = inputdlg(prompt,title,dims,definput);
% convert input information into channel array
ind_select = strsplit(answer{1},{','});
ind_select = cellfun(@str2num,ind_select)';
close all
for i_sel=1:numel(ind_select)
    
    figure; subplot 211; topoplot(Topo_ICA(:,ind_select(i_sel)),...
        Elec_loc,'plotrad',0.55,'electrodes','on');
    subplot 212; plot(abs(activations(ind_select(i_sel),:)))
    
end

% Plot the clean data

figure; subplot 211; plot(data_noise.'); ...
        subplot 212; plot(data_no_blink.')
    
% q_factor - Only to select components that demonstrate good correlation
% with overal seizure onset time - correlation with clean/denoised EEG

signal_st           = Samp_Rate*10; % 10 second-interval before seizure onset
signal_end          = size(data_sz_tot,2);
for i_fac = 1:Num_ICA
    q_factor(i_fac) = std(abs(activations(i_fac,signal_st:signal_end)))/std(abs(activations(i_fac,[1:signal_st-1])));    
end
q_sel = q_factor(ind_select)
% Is it above one?

% Select relevant Components after inspection - Finalize selection and save
% results

Comp_Sel = [3 4 10 11];
Comp_Sel = ind_select;

EEG_Sz_3_First.data = data_noise;
EEG_Sz_3_First.data_clean = data_no_blink;
EEG_Sz_3_First.Elec = Elec_loc;
EEG_Sz_3_First.Topo = Topo_ICA;
EEG_Sz_3_First.activations = activations;
EEG_Sz_3_First.weights = weights;
EEG_Sz_3_First.Comp_Sel = Comp_Sel;
EEG_Sz_3_First.comp_num = numel(Comp_Sel);

cd('../Denoised Seizures')
save('EEG_Sz_3_First.mat','EEG_Sz_3_First')
cd('../Codes')

%% Dominant Frequency Selection - Useful to determine dominant Seizure Frequency
    
Phi         = data_no_blink;
Phi_ab_1    = Phi(:,1:1000);
Phi_ab_2    = Phi(:,1001:2000);
Phi_ab_3    = Phi(:,2001:3000);
Phi_ab_4    = Phi(:,3001:4000);
Phi_ab_5    = Phi(:,4001:5000);
Phi_ab_pre  = (Phi_ab_1+Phi_ab_2+Phi_ab_3+Phi_ab_4+Phi_ab_5)/5;
Phi_ab_1    = Phi(:,5001:6000);
Phi_ab_2    = Phi(:,6001:7000);
Phi_ab_3    = Phi(:,7001:8000);
Phi_ab_4    = Phi(:,8001:9000);
Phi_ab_5    = Phi(:,9001:10000);
Phi_ab_post = (1*Phi_ab_1+1*Phi_ab_2+1*Phi_ab_3+1*Phi_ab_4+1*Phi_ab_5)/5;
NFFT        = 2^10;
figure; 
plot(linspace(-250,250,NFFT),fftshift(sum(abs(fft(Phi_ab_pre',NFFT)),2)),'r')
hold on
plot(linspace(-250,250,NFFT),fftshift(sum(abs(fft(Phi_ab_post',NFFT)),2)),'g')
xlim([0 50])

close all