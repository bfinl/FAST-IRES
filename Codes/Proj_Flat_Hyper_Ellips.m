% Code to perform flat hyperellipsoid projection. In this code we take a
% current density input which can be the output of a gradient descent
% procedure and project it inside (on the boundary) of an ellipsoid to
% satisfy the constraint of the problem. This is the generalization of
% projecting to an ellipsoid in a badly conditioned, underdetermined case.
% Abbas Sohrabpour
% 3/5/2017

function [E] = Proj_Flat_Hyper_Ellips (betta_tilda, Phi_i, phi_hat, U, D, K_U, Tp_norm,num_it)

% Output is basically the correction, or E, once added to the initial j we
% can find the projected value J_proj = J + E.
% The inputs are the shifted center or Phi_hat = Phi - K*J*T, the noise
% corrected power or betta_tilda (effect of unrelated time dynamics
% subtracted from calculated noise power, i.e. b~^2 = b^2 - n^2 where n is
% the noise defined as the independent components of the scalp potential
% not medeled by the TBF. U and D are such that K*K.' = U*D*U.' where D is
% diagnal and U is orthonormal. K_U is defined as U.'*K. if the time basis
% function (TBF) is denoted by T then TBF_norm = T.'*C_t^-1*T where C_t =
% T*T.' or the covariance. T_norm is defined as T.'*C_t^-1. Refer to my
% notes/paper for more details. The reason we pre-calculate U, D, K_U,
% TBF_norm, T_norm is to speed up the process of projection as these will
% be same for each individual problem. Phi_i is defined as the norm of a
% vector where the vector is the rows of Phi_tilda (basically summing
% effects over time). Phi_tilda is defined as Phi_tilda =
% U.'*Phi_hat*T.'*C_t^-1. num_it is the number of Newton method's
% iterations to find the optimal lambda. Refer to My notes/papers for full
% discussion.

% Elec_num = numel(diag(D));
% lambda = max(0, sqrt(Elec_num/betta_tilda)*sum(Phi_i)/sum(diag(D)) - Elec_num/sum(diag(D)));
% for it = 1:num_it
% %     lambda = (lambda - (sum((Phi_i.^1)./(1+lambda*diag(D)).^1)-sqrt(betta_tilda))/sum(-1*(Phi_i.^1).*diag(D)./((1+lambda*diag(D)).^2)));
%     lambda = (lambda + (sum((Phi_i.^1)./(1+lambda*diag(D)).^2)-(betta_tilda))/sum(-2*(Phi_i.^1).*diag(D)./((1+lambda*diag(D)).^3)));
% end
% lambda_ini = lambda
myfun = @(lambda,Phi,D_0,Betta) sum(Phi./(1+lambda*D_0).^2) - Betta;

Phi       = Phi_i;
D_0       = diag(D);
Betta     = betta_tilda;
Init_zero = 0; 

fun   = @(lambda) myfun(lambda,Phi,D_0,Betta);
lambda_ini = max(fzero(fun,Init_zero),0);

E = (K_U.'*diag(lambda_ini./(1+lambda_ini*diag(D)))*U.')*Tp_norm;

end
