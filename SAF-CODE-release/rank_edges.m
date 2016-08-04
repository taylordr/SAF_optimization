%% function [rank,Q] = rank_edges(net,w,potential_edges,visualization)
%
% Rank potential edges so that the SAF decreases the most. This maximizes the
% order parameter R
%
% Input:
%           net = struct containing the network properties
%           w = the natural frequencies
%           edges = a matrix in which the nonzero entries indicate which
%                   edges (p,q) to compute Q_{pq} for
%
% OUTPUT:   
%           Q = matrix of linear approximation terms Q_{pq}
%           rank = matrix where entry (p,q) is the rank of potential edge
%                  (p,q)
%
% Dane Taylor - August 27, 2016

function [rank,Q] = rank_edges(net,w,potential_edges,visualization)

	%disp('ranking edges')
   potential_edges = triu(potential_edges); %matrix is symmetric so only consider upper triangle

   %only consider adding unweighted edges
   epsilon = 1;%edge weight
      
   % approximate change to SAF for adding edges (p,q)
   Q = compute_Q_matrix(net,w,potential_edges);      

   %rank the entries in matrix Q

   ids = find(potential_edges); %matrix locations of potential edges
   Q_temp = Q(ids);%vector of Q values
   [ranked_Q,rank_temp] = sort(Q_temp);%sort the vector entries
   ranking = zeros(length(Q_temp),1);
   temp = length(Q_temp);
   for i=1:temp
      ranking(rank_temp(i)) = i;
   end
   
   %organize rankings into matrix
   rank = spalloc(net.N,net.N,nnz(Q));
   rank(ids) = ranking;
     
   if visualization
      figure; 
      imagesc(Q+Q');
      h = colorbar;
      ylabel(h,'$Q_{pq}$','interpreter','latex')
      colormap hot
      xlabel('node index','interpreter','latex')
      ylabel('node index','interpreter','latex')
   end
         
end