%function net = create_chain(size,visualization)
%
% Create an undirected chain containing N nodes
%
% Input:
%        size = number of nodes
%        visualization = 1/0 indicates whether or not to make figures
%
% Output:
%        net = a struct that contains all the info about the chain
%
% Dane R. Taylor - July 27, 2016

   
function net = create_chain(size,visualization)

net.N = size;%number of nodes

net.A = zeros(net.N,net.N);% create adjacency matrix and 
for n=1:(net.N-1)
   net.A(n,n+1) = 1;
   net.A(n+1,n) = 1;
end
net.M = sum(sum(net.A));%number of edges
net.D = diag(sum(net.A));%node degrees
net.L = net.D - net.A;%unnormalized Laplacian

% compute spectra of unnormalized Laplacian matrix
[net.v,net.lambdas] = eig(full(net.L));
net.lambdas = diag(net.lambdas);


%% visualize network and eigenvectors for largest and smallesst nonzero eigenvalues

if visualization
   figure;

   %define nodes coordinates for plotting
   net.coordinates = [(1:net.N)',zeros(net.N,1)];%node coordinates

   
   subplot(1,2,1)
   title('frequency vector maximizing synchrony','interpreter','latex')
   scatter(net.coordinates(:,1),net.coordinates(:,2),200*ones(net.N,1),net.v(:,net.N),'filled')
   colormap winter
   %q = max(max(abs(net.v)));
   caxis([-1,1])
   axis square
   title('$\omega=v^{(N)}$','interpreter','latex')
   xlim([0,net.N+1])   
   ylim([-1,1])
   h = colorbar;
   ylabel(h,'frequency, $\omega_m$','interpreter','latex')

   
   subplot(1,2,2)
   title('frequency vector minimizing synchrony','interpreter','latex')
   scatter(net.coordinates(:,1),net.coordinates(:,2),200*ones(net.N,1),net.v(:,2),'filled')
   colormap winter
   caxis([-1,1])
   h = colorbar;
   ylabel(h,'frequency, $\omega_m$','interpreter','latex')
   axis square
   title('$\omega=v^{(2)}$','interpreter','latex')
   xlim([0,net.N+1])
   ylim([-1,1])
   
   
   %plot the oscillator eigenvectors
   fun = @(n,m) sqrt(2/net.N)*cos(pi*(n-1)*(2*(m)*net.N-1)/2/net.N);
   MM=(net.N+1)*100;
   for n=1:net.N
      V(:,n)= fun(n,(1:MM)/(MM-100));                  
   end      
   figure;
   plot((1:MM)/MM*(net.N+1),V(:,2:end));%predicted eigenvector entries
   hold on;
   plot(net.v(:,2:end),'< ');%observed eigenvector entries
   title('eigenvectors of $L$ for chain','interpreter','latex')
   
end

end
