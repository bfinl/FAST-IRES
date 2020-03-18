% Seizure Imaging Using FATS-IRES 
% 3D - Rotational Dipole
% Abbas Sohrabpour
% 5/8/2018

% Load the lead field, etc.
close all; clear all; clc
cd('../Grid Location and Defined Parameters')
load('UnconstCentLFD.mat')
Norm_K = svds(K,1)^2;
load('Gradient.mat')
load('Laplacian.mat')
load('NewMesh.mat')
load('LFD_Fxd_Vertices.mat')
load('Edge.mat')
V  =  kron(V,eye(3));
load('Neighbor.mat')
Number_dipole = numel(K(1,:));
Location  = New_Mesh(1:3,:);
TRI = currytri.';
Vertice_Location = curryloc;
Number_edge = numel(V(:,1));

Name_Sz = 'EEG_Sz_3_First';
Var_Sz  = [Name_Sz,'.mat'];
Var_Name = [Name_Sz,'.data_clean'];

perspect = [-1 0.25 0.5];
cmap1 = jet(64);
part1 = cmap1(1:31,:);
part2 = cmap1(34:end,:);
mid_tran = gray(64);
mid = mid_tran(57:58,:);
cmap = [part1;mid;part2];

cd('../Raw 273')
Elec_loc           = readlocs('273_ele_ERD.xyz');
%% Loading The Noisy Scalp Potential and The True Underlying Source and
% Other Needed Variables, Parameters, etc.

cd('../Denoised Seizures')
load(Var_Sz)
Phi = eval(Var_Name);

Phi_noisy = Phi(:,:);
[Number_sensor,Number_Source] = size(Phi_noisy);

% Defining Parameters
% Estimating the TBF
Activ    = eval([Name_Sz,'.activations']);
Weights  = eval([Name_Sz,'.weights']);
Topo_ICA = eval([Name_Sz,'.Topo']);
Sel_Comp = eval([Name_Sz,'.Comp_Sel']);

TBF                    = Activ(Sel_Comp,:);
% TBF = ones(1,15001);
Num_TBF                = size(TBF,1);

Noise_Topo             = Topo_ICA;
Noise_Topo(:,Sel_Comp) = [];
Noise_Act              = Activ;
Noise_Act(Sel_Comp,:)  = [];
Noise_Space            = Noise_Topo*Noise_Act;

for i_top=1:numel(Sel_Comp)
    figure;
    subplot 211; topoplot(Topo_ICA(:,Sel_Comp(i_top)), Elec_loc)
    subplot 212; plot(Activ(Sel_Comp(i_top),:))
end

%% Estimating Noise
close all
cd('../Codes')
Samp_Rate  = 500;

Noise_only = Phi_noisy(:,1:4400); % Pre-seizure period, i.e. baseline
Sigma_inv_half = diag(1./std(Noise_only(:,:),[],2));
Sigma_inv_half = 1*diag(repmat(mean(diag(Sigma_inv_half)),[Number_sensor 1]));

% Initializing the Solution
Phi_noisy = Phi_noisy(:,4500:6000);
TBF       = TBF (:,4500:6000);
psi       = Phi_noisy*(TBF.')/(TBF*(TBF.'));
SNR_avg  = mean(std(Sigma_inv_half*Phi_noisy(:,1:end),[],2).^2)/mean(std(Sigma_inv_half*Noise_only(:,1:end),[],2).^2);
D_W      = ((norms(K+10^-20).^-1));
K_n      = (Sigma_inv_half*K).*repmat(D_W,[Number_sensor 1]);
psi_n    = Sigma_inv_half*psi;
[~, gam] = eig((Sigma_inv_half*K.*repmat(D_W,[Number_sensor 1]))*((Sigma_inv_half*K.*repmat(D_W,[Number_sensor 1])).'));
Smf      = 10^1;
J_ini    = (K_n.'*pinv(K_n*K_n.' + mean(diag(gam))*(Smf/SNR_avg)*eye(Number_sensor))*psi_n).*repmat(D_W',[1 Num_TBF]);

%% Reweighting Parameter
[Number_sensor,Number_Source] = size(Phi_noisy);
epsilon = 0.05;


% Noise Power Parameter (Chi_2 Dist Theory) - 
prob = 0.95;
beta = icdf('chi2',prob,Number_sensor-1);  % Spike Avg
Number_iteration = 10;

Phi_norm = Sigma_inv_half*(Phi_noisy(:,:));
K_norm   = Sigma_inv_half*(K);

power = sum((Sigma_inv_half*Noise_only).^2,1);
SNR   = norm(Phi_noisy,'fro')^2/norm(Noise_only,'fro')^2;
[X,B] = hist(power,100);
% mean(power);
% mean(power) + std(power);
% figure; histogram(power,100)
Sum_X = cumsum(X)/sum(X);
% figure; plot(Sum_X)
ind_90 = find(Sum_X > 0.90); B(ind_90(1,1));
ind_95 = find(Sum_X > 0.95); B(ind_95(1,1));


CF_min = beta/(norm(Phi_norm,'fro')^2/Number_Source) % Correcting Factor 
% for noise power as noise might not be whitened enough

CF     = max(beta/B(ind_90(1,1)),CF_min) ;

% L-curve or some other method to select alpha
% alpha_vec = [ 1  0.5 0.35 0.25 0.2 0.1 0.05 0.025 0.01 0.075 0.005 0.0025 0.001]';
alpha_vec = 0.08;
Num_alpha = numel(alpha_vec);

% Variables to save results
J_sol     = zeros(Number_dipole,Num_TBF,Number_iteration,Num_alpha);
Y_sol     = zeros(Number_edge  ,Num_TBF,Number_iteration,Num_alpha);
J_sol_2   = J_sol;
Y_sol_2   = Y_sol;

% Weighting the depth - initializing weighting matrix
M   = size(K,2);
N   = size(V,1);
W   = ones(M,Num_TBF);
W_d = ones(N,Num_TBF);

% Initialize optimization parameters - set number of iterations, etc.
Lambda_min = Norm_K*sqrt(Number_sensor)/max(norms(K.'*Phi_noisy));

C_t         = (TBF*(TBF.'));
TBF_norm    = TBF.'*(C_t\TBF);
T_norm      = (TBF.')/C_t;
alpha       = alpha_vec;
lambda      = 10.01*max(1, Lambda_min);
W_v         = (V.')*V;
L_v         = 1.1*svds(W_v,1);
x_0         = zeros(Number_dipole,Num_TBF);
eps         = 10^-3;
betta_tilda = Number_Source*(beta)/CF;
Abbas       = K_norm*(K_norm.');
[U, D, ~]   = svd(Abbas);
K_U         = U.'*K_norm;
num_it_x    = 20;
num_it_y    = 20;
num_it_new  = 20;
stop_crt    = 10^-4;

% Initialize the optimization problem, Loop through till convergence,
% update the weights using the iterative re-weighting schema and solve
% agian until convergence or counter overflow.

stop_itr            = 0;
weight_it           = 0;
max_weight_itr_num  = Number_iteration;
t1 = tic;
for i_alpha = 1:Num_alpha
    alpha = alpha_vec(i_alpha)
    x_0 = J_ini;
    stop_itr = 0;
    weight_it = 0;
    W = ones(M,Num_TBF);
    W_d = ones(N,Num_TBF);
    while (~stop_itr) 
        weight_it = weight_it + 1;
        if weight_it > max_weight_itr_num
            stop_itr = 1;
        end

        [J,Y] = FISTA_ADMM_IRES (Phi_norm, TBF, TBF_norm, T_norm, alpha, lambda, W_v, V, L_v, x_0, ...
            W, W_d, eps, betta_tilda, U, D, K_norm, K_U, num_it_x,num_it_y,num_it_new);
        x_0                          = J;
        J_sol(:,:,weight_it,i_alpha) = J;
        Y_sol(:,:,weight_it,i_alpha) = Y;
        
        J_n     = reshape(J, [3, Number_dipole/3, Num_TBF]);
        Ab_J    = squeeze(norms(J_n))+10^-20;
        Y_n     = reshape(V*J, [3, Number_edge/3, Num_TBF]);
        Ab_Y    = squeeze(norms(Y_n))+10^-20;
        
        W_old   = W;
        W_tr    = 1./(0 + Ab_J./(repmat(max(Ab_J),[Number_dipole/3 1])) + (epsilon+10^-16) );
        W       = circshift(upsample(W_tr,3),[0 0]) + circshift(upsample(W_tr,3),[1 0]) + circshift(upsample(W_tr,3),[2 0]); 
        W_d_old = W_d;
        W_d_tr  = 1./(0 + Ab_Y./(repmat(max(Ab_Y),[Number_edge/3 1]))+(epsilon+10^-16));
        W_d     = circshift(upsample(W_d_tr,3),[0 0]) + circshift(upsample(W_d_tr,3),[1 0]) + circshift(upsample(W_d_tr,3),[2 0]); 


        if (norm(W-W_old)/norm(W_old) < stop_crt && norm(W_d-W_d_old)/norm(W_d_old) < stop_crt)
            stop_itr = 1;
        end
        clc
        figure
                h1 = trisurf(TRI,Vertice_Location(1,:),Vertice_Location(2,:),Vertice_Location(3,:),sum((Ab_J),2)); colorbar
                set(h1,'EdgeColor','None', 'FaceAlpha',1,'FaceLighting','phong');
                 light_position = [3 3 1];
                 light('Position',light_position);
                 light_position = [-3 -3 -1];
                 light('Position',light_position);
                 colorbar;
                 x_max = max(abs(sum(J,2)));
                 if isnan(x_max)
                     x_max = 1;
                 end
                 caxis([-x_max x_max]);
                 view(perspect)
                 grid off
             colormap(cmap)
             cd('../Figures')
         name_fig = ['TBF_1st_Iteration_',num2str(weight_it),'.fig'];
         saveas(gcf,name_fig)
         name_fig = ['TBF_1st_Iteration_',num2str(weight_it),'.jpeg'];
         saveas(gcf,name_fig)
         cd('../Codes')
    end
    
%     close all; 
end
t_elaps = toc(t1);
t_elaps/3600

cd('../Results')
save('J_sol_1st.mat','J_sol')
save('TBF_1st.mat','TBF')
save('Phi_norm_1st.mat','Phi_norm')
save('Phi_noisy.mat','Phi_noisy')
save('Sigma_inv_half.mat','Sigma_inv_half')
Parameters. alpha         = alpha;
Parameters. betta         = betta_tilda;
Parameters. num_it_weight = max_weight_itr_num;
Parameters. alpha         = alpha_vec;
Parameters. lambda        = lambda;
Parameters. num_it_x      = num_it_x;
Parameters. num_it_y      = num_it_y;
Parameters. eps           = eps;
Parameters. stop_crt_wt   = stop_crt;
Parameters. epsilon       = epsilon;
Parameters. time          = t_elaps;
Parameters. iteration     = weight_it;
save('Parameters.mat','Parameters')
cd('../Codes')
