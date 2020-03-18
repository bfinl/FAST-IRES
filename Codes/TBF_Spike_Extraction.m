% This code cleans up the seizure data from noisy artefacts and selects the
% components that are relevant to seizure (those that demonstrate 
% correlation with seizure onset time).
% data_noise Sohrabpour
% 1/9/2018

clear all; close all; clc

% This part of the code was to load the data for our data base. We have
% saved data here so you can directly load your raw data, i.e. data_spk_tot.

pause(0.01)
SPK = 'Spike_I';
Spk_Folder = ['../Spikes/',SPK];
% Loop_Num                                = 2;
% data_spk_tot                            = [];
% for i_loop = 1 : Loop_Num
% %     Name_var = ['data_',num2str];
% %     eval([Name_var,'=data;']);
%     LoadCurryDataFile_AS
%     data_spk_tot = [data_spk_tot, data];
% end

cd('../Raw 273')
Elec_loc                                = readlocs('273_ele_ERD.xyz');
load('data_spk_tot.mat')
clearvars -except                       data_spk_tot Elec_loc Spk_Folder
cd('../Codes')
Samp_Rate                               = 500;
data_noise                              = data_spk_tot;

% Component Analysis
[A1, A2, A3]                            = svd(data_noise,'econ');
figure; plot(diag(A2)); hold on; plot(ones(size(A2,1))*mean(diag(A2)),'r')
sum(diag(A2)>=mean(diag(A2)))
pause(1)
close all
Num_ICA            = 15;
[weights, sphere,compvars,bias,signs,lrates,activations] ...
                   = runica(data_noise,'extended',-1, 'pca', Num_ICA, 'sphering','on','maxsteps',2048);
Topo_ICA           = pinv(weights*sphere);

for i_fig=0:5:Num_ICA-5
    figure;
    subplot (5,2,1); topoplot(Topo_ICA(:,i_fig+1),Elec_loc,'plotrad',0.55,'electrodes','on'); subplot (5,2,2);plot(activations(i_fig+1,:)); %line([1+wind_ind_num 1+wind_ind_num],[-5 5],'color','red','linewidth',2,'LineStyle','--'); line([1+wind_ind_num+67 1+wind_ind_num+67],[-5 5],'color','red','linewidth',2,'LineStyle','--')
    subplot (5,2,3); topoplot(Topo_ICA(:,i_fig+2),Elec_loc,'plotrad',0.55,'electrodes','on'); subplot (5,2,4);plot(activations(i_fig+2,:)); %line([1+wind_ind_num 1+wind_ind_num],[-5 5],'color','red','linewidth',2,'LineStyle','--'); line([1+wind_ind_num+67 1+wind_ind_num+67],[-5 5],'color','red','linewidth',2,'LineStyle','--')
    subplot (5,2,5); topoplot(Topo_ICA(:,i_fig+3),Elec_loc,'plotrad',0.55,'electrodes','on'); subplot (5,2,6);plot(activations(i_fig+3,:)); %line([1+wind_ind_num 1+wind_ind_num],[-5 5],'color','red','linewidth',2,'LineStyle','--'); line([1+wind_ind_num+67 1+wind_ind_num+67],[-5 5],'color','red','linewidth',2,'LineStyle','--')
    subplot (5,2,7); topoplot(Topo_ICA(:,i_fig+4),Elec_loc,'plotrad',0.55,'electrodes','on'); subplot (5,2,8);plot(activations(i_fig+4,:)); %line([1+wind_ind_num 1+wind_ind_num],[-5 5],'color','red','linewidth',2,'LineStyle','--'); line([1+wind_ind_num+67 1+wind_ind_num+67],[-5 5],'color','red','linewidth',2,'LineStyle','--')
    subplot (5,2,9); topoplot(Topo_ICA(:,i_fig+5),Elec_loc,'plotrad',0.55,'electrodes','on'); subplot (5,2,10);plot(activations(i_fig+5,:));%line([1+wind_ind_num 1+wind_ind_num],[-5 5],'color','red','linewidth',2,'LineStyle','--'); line([1+wind_ind_num+67 1+wind_ind_num+67],[-5 5],'color','red','linewidth',2,'LineStyle','--')
end

close all
% How much data variance is not explained
data_recon                              = Topo_ICA*weights*data_noise;
100*norm(data_recon - data_noise)/norm(data_noise)

% Select Components that are noisy, leave components are time locked with
% spike onset and are not clearly noisy artefacts.
Noise_comp                              = Topo_ICA;
Noise_comp (:,[9:end])                  = 0;
data_spk_clean                          = Noise_comp*weights*data_noise;

% Average the epochs
Win_Len                                 = 2001;
Num_Spk                                 = size(data_spk_tot,2)/Win_Len;
Act                                     = zeros(Num_ICA,Win_Len);
for i_epoch = 1:Num_Spk
    Act                                 = Act + ...
        activations(:,(i_epoch-1)*Win_Len+1:i_epoch*Win_Len);
end
Act                                     = Act/Num_Spk;
figure;
plot(Act')

% Plot IC with averaged time-courses (averaged time-locked to spike epochs)
close all
ind_select = [1:8];
q_factor   = zeros(numel(ind_select),1);
signal_st  = 800;
signal_end = 1250;
for i_sel=1:numel(ind_select)
    figure; subplot 211; topoplot(Topo_ICA(:,ind_select(i_sel)),...
        Elec_loc,'plotrad',0.55,'electrodes','on');
    subplot 212; plot(Act(ind_select(i_sel),:))
end

% q_factor - See which components are active around spike peak point
for i_fac = 1:Num_ICA
    q_factor(i_fac) = std(Act(i_fac,signal_st:signal_end))/std(Act(i_fac,[1:signal_st-1,signal_end:end]));    
end
q_sel = q_factor(ind_select)
q_gen = q_factor([9:end]);
mean(q_gen)
std(q_gen)
mean(q_gen)+std(q_gen)
figure; hist(q_gen,20)

% Is it strong enough?

% Select relevant Components after inspection - Finalize selection and save
% results

TBF_sel = [1,2,3,5,6];
TBF     = Act(TBF_sel,:);

% Plot Average/De-noised Spikes
Win_Len                                 = 2001;
Num_Spk                                 = size(data_spk_tot,2)/Win_Len;
Data_Avg                                = zeros(size(data_spk_tot,1),Win_Len);

for i_epoch = 1:Num_Spk
    Data_Avg                            = Data_Avg + ...
        data_spk_clean(:,(i_epoch-1)*Win_Len+1:i_epoch*Win_Len);
end
Data_Avg                                = Data_Avg/Num_Spk;
figure; 
plot(Data_Avg')

Data_Avg_proj          = Data_Avg*TBF.'*pinv(TBF*TBF.')*TBF;
N_wht                  = Data_Avg - Data_Avg_proj;
figure; plot(Data_Avg_proj')
figure; plot (N_wht')

% Save Spike Types and Parameters

Data_Avg_clean  = Data_Avg_proj;
ICA.Number      = Num_ICA;
ICA.act         = activations;
ICA.Act         = Act;
ICA.Topo        = Topo_ICA;
ICA.q_sel       = q_sel;
ICA.q_gen       = q_gen;
ICA.Noise_comp  = [9:Num_ICA];
ICA.TBF_comp    = TBF_sel;
ICA.Win_Len     = Win_Len;
cd('../Spikes')
 if ~exist(Spk_Folder, 'dir')
       mkdir(Spk_Folder)
 end
cd(Spk_Folder)

save('data_spk_tot.mat','data_spk_tot')
save('Data_Avg.mat','Data_Avg')
save('Data_Avg_clean.mat','Data_Avg_clean')
save('TBF.mat','TBF')
save('ICA.mat','ICA')
cd('../../Codes')