% Code to put the partioned data together
% Abbas Sohrabpour
% 18/3/2020

cd('../Denoised Seizures')
load('EEG_Sz_3_First_pt_1.mat')
load('EEG_Sz_3_First_pt_2.mat')

Data_pt_1 = EEG_Sz_3_First_pt_1.data;
Data_pt_2 = EEG_Sz_3_First_pt_2.data;

Data_clean_pt_1 = EEG_Sz_3_First_pt_1.data_clean;
Data_clean_pt_2 = EEG_Sz_3_First_pt_2.data_clean;

Act_pt_1 = EEG_Sz_3_First_pt_1.activations;
Act_pt_2 = EEG_Sz_3_First_pt_2.activations;

Data = [Data_pt_1, Data_pt_2];
Data_clean = [Data_clean_pt_1, Data_clean_pt_2];
Act = [Act_pt_1, Act_pt_2];

EEG_Sz_3_First = EEG_Sz_3_First_pt_1;
EEG_Sz_3_First.data = [];
EEG_Sz_3_First.data = Data;
EEG_Sz_3_First.data_clean = [];
EEG_Sz_3_First.data_clean = Data_clean;
EEG_Sz_3_First.activations = [];
EEG_Sz_3_First.activations = Act;

save('EEG_Sz_3_First.mat','EEG_Sz_3_First')

