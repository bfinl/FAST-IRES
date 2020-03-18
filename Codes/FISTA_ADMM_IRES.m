% Code to solve the spatio=temporal IRES. Basically a block coordinate
% projected gradient descent algorithm where the edges and the solution
% itself will be updated sequentially.
% Abbas Sohrabpour
% 5/5/2017

function [x, y] = FISTA_ADMM_IRES (Phi, TBF, TBF_norm, T_norm, alpha, lambda, W_v, V, L_v, x_0, W, W_d, epsilon, betta_tilda, U, D, K, K_U, num_it_x,num_it_y,num_it_new)

Alpha = diag(sum(abs(x_0))/max(sum(abs(x_0))));
y_cond = 1;

y_it = 0;
x_tild = x_0;
x_old = x_0;
y = (V*x_tild);
u = zeros(size(y));

e_rel  = 10^-4;
e_abs  = 10^-8;
tau    = 2;
mu     = 10;

while (y_cond)
    
    y_it = y_it + 1
    
    if (y_it>num_it_y)
        y_cond = 0;
    end
    
    t_old = 1;
    x_it = 0;
    x_cond = 1;
    
    while (x_cond)
        
        x_it = x_it + 1;
        
        if (x_it>num_it_x)
            x_cond = 0;
        end
        
        x_new_tran = wthresh ((x_tild - (W_v*x_tild-V.'*(y+u))/L_v),'s',(W*Alpha)*(alpha/lambda/L_v));
        phi_hat = Phi - K*x_new_tran*TBF;
        Phi_i = (norms ( ((U.')*phi_hat*TBF_norm).' ).').^2;
        Tp_norm = phi_hat*T_norm;
        N_Phi =  phi_hat - phi_hat*TBF_norm;
        N_F   = norm(N_Phi,'fro')^2;
        
        if (norm(phi_hat,'fro')^2<=betta_tilda)
            x_new_er = zeros(size(x_new_tran));
        else
            x_new_er = Proj_Flat_Hyper_Ellips (betta_tilda-N_F, Phi_i, phi_hat, U, D, K_U, Tp_norm,num_it_new);
        end
        x_new = x_new_tran + x_new_er;
        
        t_new = 0.5*(1+sqrt(1+4*t_old^2));
        
        x_tild = x_new + ((t_old-1)/t_new)*(x_new-x_old);
        if norm(x_new-x_old,'fro')/norm(x_new,'fro') < epsilon
            x_cond =0;
        end
       
        x_old = x_new;
        t_old = t_new;

    end
    
    y_new = wthresh((V*x_new)-u,'s',W_d/lambda);
    if norm(y-y_new,'fro')/norm(y_new,'fro') < epsilon
            y_cond =0;
    end
    
    prim_res = norm(y_new - V*x_new,'fro');
    dual_res = norm(lambda*(V.')*(y - y_new),'fro');
    
    e_prim = sqrt(numel(y))*e_abs + ...
        e_rel*max(sum(norms(V*x_new)),sum(norms(y_new)));
    e_dual = sqrt(numel(x_new))*e_abs + ...
        e_rel*sum(norms(V.'*y_new));
    
    if (prim_res <= e_prim && dual_res <= e_dual)
        y_cond = 0;    
    end
        
    y = y_new;
    u = u + y - V*x_new;
    
%     if     (prim_res > mu*dual_res)
%         lambda = lambda*tau;
%         u      = u/tau;
%     elseif (dual_res > mu*prim_res)
%         lambda = lambda/tau;
%         u      = u*tau;
%     end
    
end

x = wthresh(x_new,'s',(W*Alpha)*(alpha/lambda/L_v));
y = wthresh((V*x_new)-u,'s',W_d/lambda);
    
    