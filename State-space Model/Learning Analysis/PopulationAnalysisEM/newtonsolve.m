function [B, flagfail] = ...
                newtonsolve(trial, muone, Bhatold, W, obs, jval);
%newtonsolve is a helper function that implements Newton's Method in order 
%to recursively estimate the posterior mode of the augmented learning state
%process, BETA*{k|k} (B).  Once the subsequent estimates sufficiently
%converge, the function returns the last estimate.  If, having never met
%this convergence condition, the function goes through all of the
%recursions, then a special flag - indicating the convergence failure - is
%returned along with the last posterior mode estimate.
%
%variables: 
%   F            vector of equation A13*
%   pkj          p{k|k}, observation model probability (equation 2.2)*
%   obs          vector of observations
%   B            BETA*{k|k}, the posterior mode
%   g(i)         derivative of the learning state process
%   gprime(i)    derivative of g
%   it(i)        estimate of posterior mode

J  =jval(end);

flagfail = 1;
it = zeros(J+1, 40);

%do first estimate
F         = zeros(J,1);
pkj          = exp(muone + Bhatold(1)*Bhatold(jval+1))./ ...
               (1 + exp(muone+ Bhatold(1)*Bhatold(jval+1)));
F(1)      = sum( Bhatold(jval+1).*(obs(jval) - pkj) );
temp1          = obs(jval) - pkj;
F(jval+1) = temp1.*Bhatold(1);
it(:, 1)     = Bhatold + mtimes(W, F);

%now iterate
for i = 1:200

  pkj             = exp(muone + it(1,i)*it(jval+1, i))./ ...
                    (1 + exp(muone+ it(1,i)*it(jval+1,i)));
  F            = zeros(J+1,1);
  F(1)         = sum(it(2:end, i).*(obs-pkj));
  F(jval+1)    = it(1, i)*(obs-pkj);
  g               = zeros(J+1,1);
  g               = it(:,i) - Bhatold - mtimes( W, F);
  gprime          = diag(ones(1, J+1));
%first col

  for col1 = 1: J+1
   tempsum = sum( it(2:end,i).*it(2:end,i).* pkj.*(1 - pkj) );
   xx     = W(col1, 2:end)';
   yy     = it(1, i)*it(2:end,i).*pkj.*(1-pkj) - (obs -pkj);
   gprime(col1, 1)    = gprime(col1, 1)  + W(col1,1)*tempsum  + sum( xx.* yy); 
  end

%top row and col below
for row  = 1:J+1
 for col = 2: J
  temp2 = obs(col-1) - pkj(col-1)  - it(col-1,i)*it(1,i)*pkj(col-1)*(1-pkj(col-1));
  gprime(row, col) =  gprime(row, col) ...  
                     - W(row, 1)*temp2  ...
					 + W(row, col+1)*it(1,i)*pkj(col-1)*(1-pkj(col-1));                                       
 end
end

  gprimeinv       = inv(gprime);
  h               = zeros(J+1,1);
  h               = - mtimes(gprimeinv, g);
  it(:,i+1)       = it(:,i) + h;
  B               = it(:,i+1);


  if( mean(abs(B - it(:,i))) < 1e-14)
%    fprintf(2,'cvrged (1e-12) in %d iters %f \n', i, mean(abs(B - it(:,i))));
   flagfail = 0;
   return
  end

end

if(i>=200)
    fprintf(2,'failed to converge %d \n',trial)
	B=it(:,1)
    flagfail = 1;
end

