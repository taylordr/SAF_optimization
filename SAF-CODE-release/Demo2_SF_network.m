%% Example of a Scale-Free Network
%
% Optimize phase synchronization of a power-law network, by iteratively 
% adding edges based on gradient-descent algorithms that use the 
% perturbation of the synchrony alignment function (SAF).
%
% Dane R. Taylor - July 27, 2016

clear;clc;
 

%% Build the chain system

   visualization = 1;%visualize network? 1=yes
   network_size=200;%number of nodes
   gamma=2.6;
   dmin=20; 
   
   %construct network
   net_SF = create_SF(network_size,gamma,dmin,visualization);
     
   %draw frequencies from normal distribution
   net_SF.w = randn(net_SF.N,1);

   %compute original SAF
   SAF_0 = compute_SAF(net_SF.w,net_SF.L);
 

%% Rank potential edges according to perturbation of the SAF

   epsilon = 1;%only consider unweighted edges

   %possible new edges (disallowing self-edges and repeat edges)  
   potential_edges = ones(net_SF.N) - net_SF.A - eye(net_SF.N);

   %create Q_{pq}, which gives the linear approximation for SAF perturbation
   [ranks,Q_chain] = rank_edges(net_SF,net_SF.w,potential_edges,visualization);

   
%% Find five top-ranked potential new edges

   max_iter = 5;% number of edge additions
   [top_edges(:,1),top_edges(:,2)] = find( (ranks <= 5).*(ranks > 0));
   
   %print the top-ranked edges
   disp('top-ranked edges according to Q')
   top_edges

   [SAF_approx_6_1,SAF_actual_6_1] = algorithm_6_1(net_SF,net_SF.w,max_iter,SAF_0);
    
   
%% Optimize the system using Algorithm 6.2

   [SAF_approx_6_2,SAF_actual_6_2] = algorithm_6_2(net_SF,net_SF.w,max_iter,SAF_0);

   K = SAF_0*5;%coupling strength chosen so that R=0.9 before any edge additions
   R = @(SAF) 1 - SAF/(2*K); %variance order parameter given by Eq. (5.3)


   f_algorithm = figure;
   plot(0:max_iter,R(SAF_approx_6_1),0:max_iter,R(SAF_actual_6_1));
   hold on
   plot(0:max_iter,R(SAF_approx_6_2),0:max_iter,R(SAF_actual_6_2));
   legend('Algorithm 6.1 approx','Algorithm 6.1 actual','Algorithm 6.2 approx','Algorithm 6.2 actual')
   xlabel('edges added','interpreter','latex')
   ylabel('variance order parameter $R$','interpreter','latex')
   title('Performance of Algorithm 6.1 and 6.2 for power-law netowrk','interpreter','latex')

   %save('SF_experiment_data.mat')


