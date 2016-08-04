%% function [Q_approx] = compute_Q_matrix(net,w,potential_edges)
%
% Create the matrix Q, wherein entry Q_{pq} approximates the change to the
% SAF after adding a new edge (p,q). Note that removing edge (p,q)
% corresponds to -Q_{pq}.
%
% INPUT:    
%           net = struct containing the network properties
%           w = the natural frequencies
%           potential_edges = a matrix in which the nonzero entries indicate for 
%                   which edges (p,q) to compute Q_{pq} 
%           visualization = 0/1 make figures?
%
% OUTPUT:   
%           Q_approx = matrix of linear approximation terms Q_{pq}
%
% Dane Taylor - July 27, 2016


function [Q_approx] = compute_Q_matrix(net,w,potential_edges)

   [row,col] = find(triu(potential_edges)); %(p,q) for potential edges   
   
   Q_approx = spalloc(net.N,net.N,length(row));%allocate memory for Q matrix
   
   %compute Q_{pq} for each potential edge (p,q)
   for i = 1:length(row)  
      disp(i/length(row));
      Q_approx(row(i),col(i)) = ...
         change_in_SAF_under_edge_addition(w,net,row(i),col(i));       
   end    


end



