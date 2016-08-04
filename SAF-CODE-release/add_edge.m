%% net = add_edge(net,edges)
%
% Use this script to add an edge or edges to a network, and then update its
% properties
%
% INPUT: 
%        net - the struct that has the network properties
%        edges - a num_edges x 2 matrix in which a row denotes the edge to add
%
% OUTPUT: 
%        net - the struct that has the updated network properties
%
% Dane Taylor - July 27, 2016


function net = add_edge(net,edges)

   num_edges = size(edges,1);

   pp = edges(:,1);
   qq = edges(:,2);

   for i = 1:length(pp)
      p=pp(i);
      q=qq(i);
      
      net.A(p,q) = net.A(p,q) +1;
      net.A(q,p) = net.A(q,p) +1;

      net.L(p,p) = net.L(p,p) + 1;
      net.L(q,q) = net.L(q,q) + 1;
      
      net.L(q,p) = net.L(q,p) - 1;
      net.L(p,q) = net.L(p,q) - 1;

      net.D(p,p) = net.D(p,p) + 1;
      net.D(q,q) = net.D(q,q) + 1;
   end
   
   %update number of edges
   net.M = net.M + 2*num_edges;
   
   %update new spectral properties of Laplacian
   [net.v,net.lambdas] = eig(net.L);
   net.lambdas = diag(net.lambdas);


end