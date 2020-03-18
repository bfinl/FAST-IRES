% This code loads Results and looks at the peak of spikes and determines
% IZ from spike imaging results
% Abbas Sohrabpour
% 13/8/2018


clear all;close all;clc
bad_chan                = [];
cd('../Grid Location and Defined Parameters')
load('UnconstCentLFD.mat')
[~, Num_Dip] = size(K);
K([bad_chan],:)         =[];                     % Removing bad channels if need be
Norm_K                  = svds(K,1)^2;
load('Gradient.mat')
load('Laplacian.mat')
load('NewMesh.mat')
Ori                     = New_Mesh(5:7,:);
load('LFD_Fxd_Vertices.mat')
load('Edge.mat')
V                       =  kron(V,eye(3));
load('Neighbor.mat')
Number_dipole           = numel(K(1,:));
Location                = New_Mesh(1:3,:);
TRI                     = currytri.';
Vertice_Location        = curryloc;
Number_edge             = numel(V(:,1));

perspect                = [-1 0.25 0.5];
cmap1                   = jet(64);
part1                   = cmap1(1:31,:);
part2                   = cmap1(34:end,:);
mid_tran                = gray(64);
mid                     = mid_tran(57:58,:);
cmap                    = [part1;mid;part2];

cd('../Raw 273')
Elec_loc                = readlocs('273_ele_ERD.xyz');

SPK = 'Spike_I';
Spk_Folder = ['../Spikes/',SPK];
cd(Spk_Folder)
load('Data_Avg.mat')
load('TBF.mat')
cd('../')
Phi = Data_Avg;
Phi_noisy = Phi(:,:);
[Number_sensor,Number_Source] = size(Phi_noisy);
Num_TBF                = size(TBF,1);
Res_Folder = ['../Results/',SPK];
cd(Res_Folder)
load('J_sol_1st.mat')
load('TBF_1st.mat')
J                       = squeeze(J_sol(:,:,end));
Num_TBF                 = size(J,2);
J_T                     = J*TBF;
load('Phi_noisy.mat')
cd('..')

%%
cd('../Codes')
J_abbas                 = reshape(J, [3, Number_dipole/3, Num_TBF]);
J_abbas_1               = squeeze(norms(J_abbas,1));
J_abbas_2               = J_abbas_1.^2*var(TBF(:,250:350),[],2);
J_init                  = (J_abbas_2);
figure
h1 = trisurf(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),(J_init)), colorbar
set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','phong');
light_position = [3 3 1];
light('Position',light_position);
light_position = [-3 -3 -1];
light('Position',light_position);
colorbar;
x_max = max(sum(abs(J_init),2));
caxis([-x_max x_max]);
view(perspect)
grid off
colormap(cmap)
%%
cd('../Codes')    
J_col                   = J_abbas_2;
Thr                     = 0.016;
J_col(abs(J_col)<Thr*max(abs(J_col)))=0;
IND                     = Find_Patch_V2(Edge, 1, J_col);
J_init                  = IND;
figure
h1 = trisurf(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),(J_init)), colorbar
set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','phong');
light_position = [3 3 1];
light('Position',light_position);
light_position = [-3 -3 -1];
light('Position',light_position);
colorbar;
x_max                   = max(sum(abs(J_init),2));
caxis([-x_max x_max]);
view(perspect)
grid off
colormap(cmap)

J_seg = IND;
cd('../Resection or SOZ data')
save('J_sol_Spk_Seg.mat','J_seg')
save('J_sol_Spk.mat','J_abbas_2')
%% Overlap
% Load clinical data from Folder and compare to J_seg to calculate NORs

