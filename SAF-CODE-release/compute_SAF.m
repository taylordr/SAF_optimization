%function SAF = compute_SAF(w,L)
%
% Compute the synchrony alignment function (SAF) for given w and L
%
% Input:
%        w = frequency vector
%        L = unnormalized Laplacian matrix
%
% Output:
%        SAF = synchrony alignment function
%
% Dane R. Taylor - July 27, 2016

function SAF = compute_SAF(w,L)
   
   %calculate eigenvectors and eigenvalues of L
   [v,lambdas] = eig(full(L));
   lambdas = diag(lambdas);

   
   SAF = 0;
   for n=2:length(L)
      SAF = SAF + lambdas(n)^-2 * (w'*v(:,n))^2;
   end
   SAF = SAF/length(L);

end


