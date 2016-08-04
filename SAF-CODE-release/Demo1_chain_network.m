%% Pedagogical Example of a Chain Network
%
% Optimize phase synchronization of a chain network, by iteratively 
% adding edges based on gradient-descent algorithms that use the 
% perturbation of the synchrony alignment function (SAF).
%
% Dane R. Taylor - July 27, 2016

clear;clc;
 

%% Build the chain system

   visualization = 1;%visualize network? 1=yes
   network_size=9;%number of nodes

   %construct network
   net_chain = create_chain(network_size,visualization);

   
   
   %draw frequencies from normal distribution
   net_chain.w = randn(net_chain.N,1);

   %compute original SAF
   SAF_0 = compute_SAF(net_chain.w,net_chain.L);


%% Rank potential edges according to perturbation of the SAF

   epsilon = 1;%only consider unweighted edges

   %possible new edges (disallowing self-edges and repeat edges)  
   potential_edges = ones(net_chain.N) - net_chain.A - eye(net_chain.N);

   %create Q_{pq}, which give linear approximation for SAF perturbation
   [ranks,Q_chain] = rank_edges(net_chain,net_chain.w,potential_edges,visualization);

%% Find five top-ranked potential new edges

   max_iter = 5;% number of edge additions
   [top_edges(:,1),top_edges(:,2)] = find( (ranks <= 5).*(ranks > 0));
   
   %print the top-ranked edges
   %disp('top edges according to Algorithm 6.1')
   %top_edges

   [SAF_approx_6_1,SAF_actual_6_1] = algorithm_6_1(net_chain,net_chain.w,max_iter,SAF_0);
    
%% Optimize the system using Algorithm 6.2

   [SAF_approx_6_2,SAF_actual_6_2] = algorithm_6_2(net_chain,net_chain.w,max_iter,SAF_0);

   K = SAF_0*5;%coupling strength chosen so that R=0.9 before any edge additions
   R = @(SAF) 1 - SAF/(2*K); %variance order parameter given by Eq. (5.3)

   %NOTE THAT THE 1ST ORDER APPROXIMATIONS ARE ACCURATE FOR SMALL
   %PERTURBATIONS, WHICH IS NOT THE CASE FOR THIS SMALL CHAIN NETWORK.
   %NEVERTHELESS, THE OPTIMIZATION ALGORITHMS ARE VERY EFFECTIVE.

   f_algorithm = figure;
   plot(0:max_iter,R(SAF_actual_6_1),0:max_iter,R(SAF_actual_6_2));
   legend('Algorithm 6.1','Algorithm 6.2')
   xlabel('edges added','interpreter','latex')
   ylabel('variance order parameter $R$','interpreter','latex')
   title('Performance of Algorithm 6.1 and 6.2 for chain','interpreter','latex')

   %save('chain_experiment_data.mat')


