%function net = create_SF(type,size,visualization)
%
% Create a a scale free network using Chung Lu model
%
% Input:
%        size = number of nodes
%        visualization = 1/0 indicates whether or not to make figures
%
% Output:
%        net = a struct that contains all the info about the chain
%
% Dane R. Taylor - July 27, 2016

      
function net = create_SF(size,gamma,dmin,visualization)
      
   net.N = size;%number of nodes
   net.gamma = gamma;%power law exponent
   
   %generate an expected degree sequance
   d = sort(floor(dmin*(1-rand(net.N,1)).^(-1/(gamma-1))),'descend');            
   d(find(d>net.N-1)) = net.N-1;%upper bound degree at N-1
   
   M = 2*sum(d);%expected number of edges
   
   %adjacency matrix
   net.A = zeros(size);
   for n=1:(size-1) 
     m = (n+1):size;
         scalars = d(n)*d(m)'./ (M/2)*(size)/(size-1);
         net.A(n,m) = rand(1,length(m))<scalars;
   end
    
   net.A = net.A+net.A'; 
   net.M = sum(sum(net.A));
   net.D = diag(sum(net.A));
   net.L = net.D-net.A;
   [net.v,net.lambdas] = eig(net.A);
   net.lambdas=diag(net.lambdas);
   
   if visualization
      figure;
      spy(net.A);
   end
   
end

