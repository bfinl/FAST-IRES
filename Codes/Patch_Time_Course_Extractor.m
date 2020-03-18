% This codes extracts patches from FAST-IRES solution and then extracts
% time courses from it and then performs connectiviy and if need be a
% second round of source localization. 
% Abbas Sohrabpour
% 25/10/2018 - New Version automatically chooses/groups regions (I hereby
% acknowledge Chris Cline for informing me about linkage function in MATLAB
% and using this idea for grouping).
% 30/10/2018 - In this new version I have incorporated a newer version of
% finding patch which was modified by Zhengxianf Cai. Zhengxiang has also
% modified the format of this version to some extnt. Thanks to Zhengxiang.

%% cleaning and loading
clear all; close all; clc
% setup data folders
cd('../Codes')

cd('../Grid Location and Defined Parameters')
load('UnconstCentLFD.mat')
[~, Num_Dip] = size(K);
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

cd('../Results')
load('J_sol_1st.mat')
load('TBF_1st.mat')
J                       = squeeze(J_sol(:,:,end));
Num_TBF                 = size(J,2);
J_T                     = J*TBF;
load('Phi_noisy.mat')

cd('../Codes')

for i_tbf = 1:Num_TBF
    J_init              = J(:,i_tbf);
    J_init              = norms(reshape(J_init, [3 Num_Dip/3]))';
    figure; 
    h1 = trisurf(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),(J_init)); colorbar
    set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','phong');
    light_position = [3 3 1];
    light('Position',light_position);
    light_position = [-3 -3 -1];
    light('Position',light_position);
    colorbar;
    x_max               = max(sum(abs(J_init),2));
    caxis([-x_max x_max]);
    view(perspect)
    grid off
    colormap(cmap)
end

J_init                  = zeros(size(J_init));
for i_tbf = 1:size(J,2)
    J_tran              = J(:,i_tbf);
    J_tran              = norms(reshape(J_tran, [3 Num_Dip/3]))';
    J_init              = J_init + J_tran/max(J_tran);
end
figure; 
h1 = trisurf(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),(J_init)); colorbar
set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','gouraud');
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

%% Now How Many Patches for Each Solution
cd('../Codes')
J_col                   = squeeze(norms(reshape(J, [3 Num_Dip/3 Num_TBF])));
Num_Patch_max           = 50;           % Max Number of Patches(not really important as the find patch will break after it becomes empty, due to threshold and IRES piecewise continuous nature)
Thr                     = 0.016;        % Discard weaker patches below, i.e. background (1/64 - gray background)
IND_tot                 = cell(1,Num_TBF);
nPatch                  = zeros(1,Num_TBF);

% Eliminate background activity
for i_tbf = 1:Num_TBF
    J_col(J_col(:,i_tbf) < Thr*max(J_col(:,i_tbf)),i_tbf) = 0;
end
% Find num of patches for each column of J (solution), automatically
for i_tbf = 1:Num_TBF
    [IND_tot{i_tbf},nPatch(i_tbf)] = Find_Patch(Edge, Num_Patch_max, J_col(:,i_tbf));
end

for i_tbf = 1:Num_TBF
    
    for iIND = 1:size(IND_tot{i_tbf},2)
        J_init                  = IND_tot{i_tbf}(:,iIND);
        figure
        h1 = trisurf(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),...
            (J_init)); colorbar;
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
    end
    
end

%% Find Patches (for each level)
IND                     = cat(2,IND_tot{:});
IND_sum                 = sum(IND,2);
Max_Lev                 = max(IND_sum);
IND_sep                 = zeros(Num_Dip/3,Max_Lev);
close all
for i_ind = 1:size(IND_sep,2)
    
    IND_sep(IND_sum==i_ind,i_ind)=1;
    J_plot              = IND_sep(:,i_ind);
    figure; 
    h1 = trisurf(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),(J_plot)); colorbar
    set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','phong');
    light_position = [3 3 1];
    light('Position',light_position);
    light_position = [-3 -3 -1];
    light('Position',light_position);
    colorbar;
    x_max               = max(sum(abs(J_plot),2));
    caxis([-x_max x_max]);
    view(perspect)
    grid off
    colormap(cmap)
    
end
%% find num of patches for each level
J_col                           = IND_sep;
Num_Patch_max                   = 100;     % Max Number of Patches (not really important as the find patch will break after it becomes empty, due to threshold and IRES piecewise continuous nature)
IND_tot                         = cell(1,Max_Lev);
nPatch                          = zeros(1,Max_Lev);

% Find patches for each column of J (solution)
for i_lev = 1:Max_Lev
    [IND_tot{i_lev},nPatch(i_lev)] = Find_Patch(Edge, Num_Patch_max, J_col(:,i_lev));
end

for i_lev = 1:Max_Lev
    
    for iIND = 1:size(IND_tot{i_lev},2)
        J_init                  = IND_tot{i_lev}(:,iIND);
        figure
        h1 = trisurf(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),...
            (J_init)); colorbar;
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
    end
    
end

%% Save Results if you want - Go from one dimensional to full 3D/ just in case

IND_tot = cat(2,IND_tot{:});
IND_tot_rot = circshift(upsample(IND_tot,3),[0 0]) + ...
              circshift(upsample(IND_tot,3),[1 0]) + ...
              circshift(upsample(IND_tot,3),[2 0]);

% cd(',,/Connectivity')
% save('IND_sep.mat','IND_sep')
% save('IND_sum.mat','IND_sum')
% save('IND_tot.mat','IND_tot')
% cd('../Codes')
%% Extracting Time Course for these regions

for i_ind = 1:size(IND_tot,2)
    J_tran          = J_T(IND_tot_rot(:,i_ind)>0,:);
    T_mean(:,i_ind) = mean(J_tran)';
end

%% Time Series Formation
close all
Ab = corrcoef(T_mean);
figure; imagesc(abs(Ab)); colormap(jet); colorbar

% Automatically regroup and cluster regions with high correlation together
% so that we will not have too many regions with almost the same
% time-course of activity which is an issue when calculating DTF or any
% connectivity - Abbas 25/10/2018. I have been inspired by a code Chris
% Cline passed to me in March and some reading of the linkage and cluster
% function in MATLAB and partially used his code with some minor changes
% so many thanks to Chris Cline - Abbas 25/10/2018.

Thr_con       = 0.5;
Max_Clust_Num = Num_TBF;
L = linkage(T_mean','single',{@(xi, xj) 1-abs(corr(xi',xj'))});
% C = cluster(L,'cutoff',Thr_con,'criterion','distance')
C = cluster(L,'maxclust',Max_Clust_Num)
[groups,ind]=sort(C)
figure; imagesc(abs(Ab(ind,ind))); colormap(jet); colorbar

IND_grp                 = zeros(size(IND_tot,1),Max_Clust_Num); 
for i_grp=1:Max_Clust_Num
    IND_grp_tran        = sum(IND_tot(:,ind(groups==i_grp)),2);
    IND_grp(:,i_grp)    = IND_grp_tran;
    figure; 
    h1 = trisurf(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),(IND_grp_tran)); colorbar
    set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','phong');
    light_position = [3 3 1];
    light('Position',light_position);
    light_position = [-3 -3 -1];
    light('Position',light_position);
    colorbar;
    x_max               = max(IND_grp_tran);
    caxis([-x_max x_max]);
    view(perspect)
    grid off
    colormap(cmap)
end
 
% Old version for the code which did this by hand, i.e. the ind had to be
% manually put in by looking at absolute value of correlation.

% ind = [1, 2, 4,  3,6,8, 5,7];            
% Ab = Ab(ind,ind);
% figure; imagesc(abs(Ab)); colormap(jet)

% Potetnially I could remove areas that have very small area as they are
% not physiologically possible. Most, based on experience, will regroup
% with larger areas near by. This might be interesting in the initial
% portion of the code when we are finding patches, as sometimes extremely
% small regions may have to be extracted before getting to larger regions.
% Point that can be considered for further automatization of the code, in
% the future - Abbas.

T_norm = repmat(norms(T_mean),[size(T_mean,1) 1]);
T_n = T_mean./T_norm;

TBF_connect = zeros(size(T_mean,1), Num_TBF);
for i_grp=1:Max_Clust_Num
    TBF_connect(:,i_grp) = sum(T_n(:,ind(groups==i_grp)),2);
end
%% Connectivity Analysis
t_start             = 250;
delta_t             = 250;
t_end               = t_start + delta_t;
% Order Selection
ts                  = (TBF_connect(t_start:t_end,:));
pmin                = 1;
pmax                = floor(size(ts,1)/(size(ts,2)+1))-1
selector            = 'sbc';
no_const            = [];
[~,~,~,sbc,fpe,~]   = arfit(ts,pmin,pmax,selector);
order_sbc           = find(sbc==min(sbc))
order_fpe           = find(fpe==min(fpe))
pause(2)
% Actual Analysis
p                   = order_fpe;
fs                  = 500;
low_freq            = 1;
high_freq           = 15;
shufftimes          = 1000;
siglevel            = 0.05;

ts                  = (TBF_connect(t_start:t_end,:));
gamma2_set_tran     = DTF(ts,low_freq,high_freq,p,fs);
new_gamma2_tran     = DTFsigvalues(ts,low_freq,high_freq,p,fs,shufftimes,siglevel,[]);
gamma2_sig_tran     = DTFsigtest(gamma2_set_tran,new_gamma2_tran);

figure; imagesc(mean(gamma2_sig_tran,3)); colormap(jet); colorbar

[n_ROI, ~, Num_Freq] = size(gamma2_sig_tran);
figure;
for i_roi=1:n_ROI
    for j_roi=1:n_ROI
        subplot (n_ROI,n_ROI,j_roi+n_ROI*(i_roi-1)), plot(squeeze(gamma2_sig_tran(i_roi,j_roi,:)))
        ylim([0 1])
    end
end

%% Determine the driving node by the column that causally drives other nodes and/or has the strongest outflow
IND_src                 = IND_grp(:,3);

Parameter.t_start       = t_start;
Parameter.delta_t       = delta_t;
Parameter.p_max         = pmax;
Parameter.order_sbc     = order_sbc;
Parameter.order_fpe     = order_fpe;
Parameter.select        = p;
Parameter.low_freq      = low_freq;
Parameter.high_freq     = high_freq;
Parameter.shufftimes    = shufftimes;
Parameters.siglevel     = siglevel;

cd('../Connectivity')
save('IND_grp.mat'      ,'IND_grp')
save('IND_src.mat'      ,'IND_src')
save('Parameters.mat'   ,'Parameters')