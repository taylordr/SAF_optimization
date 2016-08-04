%function [SAF_approx,SAF_actual] = algorithm_6_2(net,w,max_iter,SAF_0)
%
% Create an undirected chain containing N nodes
%
% Input:
%           net = struct containing the network properties
%           w = the natural frequencies
%           max_iter = number of edges to add
%           SAF_0 = original synchrony alignment function
%
% Output:
%        SAF_approx = estimated SAF using perturbation theory
%        SAF_approx = the actual new SAF after adding an edge
%
% Dane R. Taylor - July 27, 2016

function [SAF_approx,SAF_actual] = algorithm_6_2(net,w,max_iter,SAF_0)

disp('running Algorithm 6.2');

%% Define parameters and allocate memory

visualization = 0;

SAF_approx = zeros(1,max_iter+1);
SAF_approx(1) = SAF_0;

SAF_actual = zeros(1,max_iter+1);
SAF_actual(1) = SAF_0;

net_B = net;
clear net;


%% Iteratively add top-ranked edge

for i=1:max_iter
   %disp(i/max_iter)

   potential_edges = triu(ones(net_B.N) - net_B.A - eye(net_B.N));

   [rank,Q_approx] = rank_edges(net_B,w,potential_edges,visualization);% rank the potential edges
   [r2,c2] = find(rank==1);% find top ranked edge
   net_B = add_edge(net_B,[r2,c2]);% add top-ranked edge to network
   
   SAF_actual(1+i) = compute_SAF(w,net_B.L);
   SAF_approx(1+i) = SAF_actual(i) +  Q_approx(r2,c2);

   disp(['added edge (',num2str(r2),',',num2str(c2),')'])
end                  
      

end
      

