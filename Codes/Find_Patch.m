% Finding Patches from IRES Solution
% Abbas Sohrabpour
% 20/2/2018
% Zhengxiang Cai
% 30/10/2018
% This version has been vastly improved by removing extra for loops by
% using cell arrays and using the unique command and more importantly by
% defining a global pointer that avoids checking all global indices by 
% defining a pointer that remember how far in the list it has checked.
% Thanks to Zhengxiang for these changes.

function [IND,Num_TBF] = Find_Patch(Edge, Num_TBF, J_col)


Edge_mod        = Edge(:,[3,4]);
Number_src      = max(Edge_mod(:));
IND             = zeros(Number_src, Num_TBF);
J_amp           = J_col;

for i_clst = 1:Num_TBF
    
    [~, ind_max]        = max(J_amp); 
    IND(ind_max,i_clst) = 1;
    ind_neigh_local     = ind_max;
    ind_neigh_global    = zeros(Number_src,1);
    pntInd_neigh_global = 1;
    
    while (true)
        
        nNeighLocal     = length(ind_neigh_local);
        ind_neigh       = cell(nNeighLocal,1);
        
        for i_ind       = 1:nNeighLocal
            % find immediate neighbours of this node
            ind_neigh{i_ind} = ...
                [Edge_mod(Edge_mod(:,1)==ind_neigh_local(i_ind),2);...
                 Edge_mod(Edge_mod(:,2)==ind_neigh_local(i_ind),1)];
        end
        ind_neigh = unique(cat(1,ind_neigh{:}));
        
        % discard locations where solution amplitude is zero
        ind_neigh(J_amp(ind_neigh)==0)  = [];
        nNeigh = length(ind_neigh);
        %remove repeated indices, comparing with the global index
        for i_nei = nNeigh:-1:1
            if ~isempty(find(ind_neigh_global(1:pntInd_neigh_global-1) == ind_neigh(i_nei), 1))
                ind_neigh(i_nei) = [];
            end
        end
        
        nNeigh = length(ind_neigh);
        % store local neighbour results to global index/registry
        ind_neigh_global(pntInd_neigh_global-1+(1:nNeigh)) = ind_neigh;
        pntInd_neigh_global = pntInd_neigh_global + nNeigh;
        % set boundary nodes to current local neighbours so that in next
        % loop we look for the neighbors of these "frontier" nodes
        ind_neigh_local = ind_neigh;
        
        % if no new neighbors were found, break
        if isempty(ind_neigh)
            break
        else
            IND(ind_neigh,i_clst)     = 1;
        end
    end
    J_amp(IND(:,i_clst)>0)            = 0;
    
    % return if no more large patches of activity were detected
    if isempty(find(J_amp,1)) && i_clst<Num_TBF
        Num_TBF = i_clst;
        IND(:,i_clst+1:end) = [];
        break
    end
    
end
  