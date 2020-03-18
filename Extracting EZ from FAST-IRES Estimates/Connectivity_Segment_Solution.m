% This code loads Results from Connectivity Analysis and compares it to
% resection and SOZ if quantative clinical data is available
% Abbas Sohrabpour
% 19/3/2018

clear all; close all; clc
bad_chan                = [];
cd('../Grid Location and Defined Parameters')
load('UnconstCentLFD.mat')
[~, Num_Dip] = size(K);
K([bad_chan],:)         =[];                     % Removing Nz due to incomplete MRI
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

cd('../Connectivity')
load('IND_src.mat')

J_seg                   = IND_src;


figure
h1 = trisurf(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),(J_seg)), colorbar
set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','phong');
light_position = [3 3 1];
light('Position',light_position);
light_position = [-3 -3 -1];
light('Position',light_position);
colorbar;
x_max = max(sum(abs(J_seg),2));
caxis([-x_max x_max]);
view(perspect)
grid off
colormap(cmap)


cd('../Resection or SOZ data')
save('J_sol_Con_Seg.mat','J_seg')
%% Overlap
% Load clinical data from Folder and compare to j_seg to calculate NORs ... 


