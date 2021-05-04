function [Q] = getQ(Bstarnew, I, sigE, sigB, Bnew, covterm3, covterm2, background_prob, Wnew)
%getQ is a helper function that computes Q, the expectation of the complete
%data log-likelihood of the SSRE model, given the observations and the
%estimated state-space random effects model parameters.  
%
%variables:
%   pkj          probability of a correct response on trial k from subject j (equation 2.2)*      
%   part1        the first part of equation A2*
%   part2        the second part of equation A2*
%   part3        the third part of equation A2*
%   Q            the expectation of the complete data log-likelihood of the SSRE model

muone = log(background_prob/(1-background_prob));
K     = size(Bstarnew,2) - 1;
J     = size(Bstarnew,1) - 1;

part1 = 0;

for k = 2:K+1
 for j = 1:J

  jval = j+1;

  pjk   = exp(muone + Bstarnew(1, k)*Bstarnew(jval, k))./  ...
          (1 + exp(muone + Bstarnew(1, k)*Bstarnew(jval, k)));
  part1 = part1 + I(j,k-1)*log(pjk) + (1-I(j, k-1))*log(1-pjk);

 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part2a = - (K+1)/2*log(2*pi*sigE^2)  ;


term1      = 2*sum(covterm2);
term2      = 2*sum(covterm3);
temp        = (term1 - term2 + 2.0*Bstarnew(1,2)^2  + 2*Wnew(1,1,2) ...
               - Bstarnew(1,end)^2 - Wnew(1,1,end));





part2b = - 0.5/(sigE^2)*temp; 

part2  = part2a + part2b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

part3a  = -J/2*log(2*pi*sigB^2);

part3b  = -0.5/(sigB^2) * (sum(covterm3(2:end,end))      ...
                   -2*Bnew*sum(Bstarnew(2:end,end))  ...
                   + J*Bnew^2);


part3 = part3a + part3b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q = part1 + part2 + part3;