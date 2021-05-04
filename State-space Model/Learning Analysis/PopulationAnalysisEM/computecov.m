function [covterm3, covterm2] = computecov(Wnew, A, Bstarnew)
%computeCov is a helper function that computes the covariance terms R{K|K},
%which is squares, and R{k-1,k|K}, which is covariances between consecutive
%estimates
%
%variables:
%   covest       covariance estimate W{k,u|K} (equation A19)*
%   covterm2     covariance term R{K|K} (equation A21)*
%   covterm3     covariance term R{k-1,k|K} (equation A22)*

N    = length(Bstarnew);
covterm2 = []; 
covterm3 = []; 


for time = 3:N              

  covest = mtimes(A(:,:,time-1), Wnew(:,:,time));
  
  temp1   = diag(Wnew(:,:,time)) + Bstarnew(:,time).*Bstarnew(:,time);
  covterm2      = [covterm2 temp1];
  
  temp2 = diag(covest) + Bstarnew(:,time-1).*Bstarnew(:,time);
  covterm3    = [covterm3 temp2]; 

end