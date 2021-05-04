function [Bnew, Wnew, A] = backwardfilter(B, Bold, W, Wold);
%backwardfilter is a helper function that implements the fixed-interval backward 
%filter smoothing algorithm to estimate the augmented learning state
%process as a Gaussian random variable with mean BETA*{k|K} (Bnew) and
%W{k|K} (Wnew) at each trial for each subject.  
% 
%variables:
%   B            BETA*{k|k}, posterior mode (equation A11)*
%   Bold         BETA*{k|k-1}, one-step prediction (equation A9)*
%   W            W{k|k}, posterior variance (equation A12)*
%   Wold         W{k|k-1}, one-step prediction variance (equation A10)*
%   T            total number of posterior mode estimates (K + 1)
%   A(i)         A(k), (equation A17)*
%   Bnew         BETA*{k|K}, backward estimate of the augmented learning state given all the data (equation A16)*
%   Wnew         W{k|K}, backward estimate of the augmented learning state variance (equation A18)*

T = size(B,2);

%Initial conditions: use values of posterior mode and posterior variance
Bnew      = zeros(size(B));
Wnew      = zeros(size(W));

A             = zeros(size(W,1), size(W,2), size(W,3)-1);

Bnew(:,T)     = B(:,T);
Wnew(:,:,T) = W(:,:,T);

for i = T-1 :-1: 1
    %for each posterior mode prediction, compute new estimates given of all 
    %the data from the experiment (estimates from ideal observer)
    
    A(:,:,i)        = mtimes( W(:,:,i), inv(Wold(:,:,i+1)) );                
    
    temp1             = mtimes( A(:,:,i), Bnew(:,i+1) - Bold(:,i+1) );
	Bnew(:,  i)     = B(:, i) + temp1; 
    
    temp2            = mtimes( A(:,:,i), Wnew(:,:,i+1)-Wold(:,:,i+1));
	Wnew(:,:,i) = W(:,:,i) + mtimes(temp2, A(:,:,i)'); 

end