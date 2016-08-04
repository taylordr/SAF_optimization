%% function [deltaSAF] = change_in_SAF_under_edge_addition(w,net,p,q)
%
% Create the matrix Q, wherein entry Q_{pq} approximates the change to the
% SAF after adding a new edge (p,q). Note that removing edge (p,q)
% corresponds to -Q_{pq}.
%
% INPUT:    
%           net = struct containing the network properties
%           w = the natural frequencies
%           p = node 1 for edge
%           q = node 2 for edge
%
% OUTPUT:   
%           deltaSAF = change in the synchrony alignment function 
%
% Dane Taylor - July 27, 2016

function [deltaSAF] = change_in_SAF_under_edge_addition(w,net,p,q)
      
      
%% First compute the variables that appear several times
   
   %inner product between frequency vector and eigenvector
   a = zeros(net.N,1);
   for i=1:net.N;
      a(i) = w'*net.v(:,i);
   end
   
   %terms describing v'\Delta L v when edge pq is added
   b = zeros(net.N);
   for i=1:(net.N-1);
      for j=(i+1):net.N;
         b(i,j) = (net.v(p,j)-net.v(q,j)) * (net.v(p,i)-net.v(q,i));
      end
   end
   b=b+b';
          
%% Now compute Eq. (4.3)

   temp2 = 0;
   for i=2:net.N
      temp1 = 0;
      for j =setdiff(1:net.N,i)
         temp1 = temp1 + a(j) * b(i,j) ...
            / (1-net.lambdas(j)/net.lambdas(i));
      end      
      temp1 = temp1 - a(i) * (net.v(p,i)-net.v(q,i))^2;
      temp2 = temp2 + temp1 * (a(i) / net.lambdas(i)^3);
   end
   
   deltaSAF = temp2 * 2/net.N;

        
   
 
end


