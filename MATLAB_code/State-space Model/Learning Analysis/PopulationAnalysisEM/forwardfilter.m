function [Bhat, W, Bhatold, Wold, number_fail] = forwardfilter(I, sigE, muone, Bstart, Wstart)
%forwardfilter is a helper function that implements the forward recursive 
%filtering algorithm to estimate the augmented learning state process at 
%trial k for subject j as the Gaussian random variable with mean 
%BETA*{k|k-1} (Bhatold) and variance W{k|k-1} (Wold) and as the Gaussian random 
%variable with mean BETA*{k|k} (Bhat) and variance W{k|k} (W).  In
%order to calculate Bhat, the function calls another function, newtonsolve,
%which estimates its value using Newton's method.
% 
%variables:
%   Bhatold         BETA*{k|k-1}, one-step prediction (equation A9)*
%   Wold            W{k|k-1}, one-step prediction variance (equation A10)*
%   Bhat            BETA*{k|k}, posterior mode (equation A11)* 
%   W               W{k|k}, posterior variance (equation A12)*
%   pkj             p{k|k}, observation model probability (equation 2.2)*
%   Winv            W{k|k-1}^-1 (used in equation A12)*
%   G               G{k}  (equation A14)*
%   noise_var       W{BETA*}, covariance matrix (used in equation A10)*
%   J               total number of subjects
%   length_trials   total number of trials
%   number_fail     saves the time steps if Newton's Method fails

length_trials = size(I,2);
J             = size(I,1);

Bnew = zeros(J+1, length_trials);
Snew = zeros(J+1, J+1, length_trials);

%Initial conditions: use values from previous iteration
Bhat(1:J+1,1)     = Bstart;                  
W(1:J+1,1:J+1, 1) = Wstart; 

temp        = zeros(1, J);
noise_var   = diag([sigE^2 temp]);
number_fail = [];

for i = 2:length_trials+1
   %for each trial (i = k + 1) compute estimates of the one-step prediction, the
   %posterior mode (using Newton's Method), and the posterior variance
   %(estimates from subject's POV)
   
   jval = 1:J;

   %Compute the one-step prediction estimate of mean and variance         
   Bhatold(1:J+1, i)          = Bhat(1:J+1, i-1);
   Wold(1:J+1, 1:J+1,  i) = W(1:J+1,  1:J+1, i-1) + noise_var;

   %Assign newtonsolve input variables 
   startB   = Bhatold(1:J+1, i);
   startW = Wold(1:J+1, 1:J+1, i);
   
   %Use Newton's Method to compute the nonlinear posterior mode estimate
   [Breturn, flagfail] = ...
                newtonsolve(i, muone, startB, startW, I(1:J, i-1), jval);
     
   if flagfail>0 %if Newton's Method fails, number_fail saves the time step
        number_fail = [number_fail i];
		flagfail;
   end

   %calculate and assign variables
   Bhat(1:J+1, i) = Breturn;
   Winv = inv(Wold(1:J+1, 1:J+1,  i));
   G = zeros(J+1,J+1); %Initialize matrix of A14 to zeros
   pjk  = exp(muone + Bhat(jval+1,i).*Bhat(1,i))./(1+exp(muone + Bhat(jval+1,i).*Bhat(1,i)));
   G(1, 1) = sum(- Bhat(jval+1, i).^2 .* pjk.*(1-pjk)); %1st condition of A14
   G(1, jval+1)      = (I(jval, i-1) - pjk - Bhat(1,i).*Bhat(jval+1,i).*pjk.*(1-pjk))'; %1st part of 2nd condition of A14*
   G(jval+1, 1)      = I(jval, i-1) - pjk - Bhat(1,i).*Bhat(jval+1,i).*pjk.*(1-pjk); %2nd part of 2nd condition of A14* 
   for dd = 1:J 
    G(dd+1, dd+1) = -Bhat(1,i)^2*pjk(dd)*(1-pjk(dd)); %3rd condition of A14
   end
   W(1:J+1, 1:J+1, i) = inv( Winv - G ); %A12
end
