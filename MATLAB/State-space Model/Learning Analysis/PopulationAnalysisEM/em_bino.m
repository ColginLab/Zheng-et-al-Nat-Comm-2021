function [sigsqEnew, Bnew, sigsqBnew] = em_bino(covterm3, covterm2, Bstarnew, Wnew, J);
%em_bino is a helper function that computes the estimated learning process
%variance SIG_EPSILON^2 (sigsqEnew), BETA{0} (Bnew), and SIG_BETA^2
%(sigsqBnew).  
%
%variables:
%   Bstarnew     BETA*{k|K}, backward estimate of learning state
%   Wnew         SIG^2{k|K}, backward estimate of learning state variance  
%   term1        R{k|K}       (equation A20)*
%   term2        R{k-1,k|K}   (equation A22)*
%   term3        R{1|K}     (applies equation A20)* 
%   term4        R{K|K}     (applies equation A21)*
%   sigsqEnew    SIG_EPSILON^2, estimate of learning state variance from EM (equation A23)*
%   Bnew         BETA{0} (equation A24)*
%   sigsqBnew    SIG_BETA^2 (equation A25)*

K           = size(covterm2(1,:),2) + 1;

term1      = covterm2(1,:);
term2      = covterm3(1,:);
term3      = Wnew(1,1,2) + Bstarnew(1,2)^2;
term4      = Wnew(1,1,end) + Bstarnew(1,end)^2;

sigsqEnew   = ((2*(sum(term1) - sum(term2))) + 2.0*term3 - term4)/K;

Bnew         = sum(Bstarnew(2:end,end))/J;

sigsqBnew    = (    sum(covterm3(2:end,end))             ...
                    -   2*Bnew*sum(Bstarnew(2:end,end))  ...
                    +   J*Bnew^2  ) / J;