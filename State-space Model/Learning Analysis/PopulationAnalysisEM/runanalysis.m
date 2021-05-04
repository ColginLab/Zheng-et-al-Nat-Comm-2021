%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to run the analysis of the population learning curve using
% Hierarchical EM for cohort binary behavioral data.
% Anne Smith, May 2004
% 
% updated Leo Walton, July 2006 - added comments, changed variable names
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Given a set of behavioral experiment trial data from multiple subjects,
% this program estimates an unobservable learning state process,
% defined as a random walk, for the population and the individual subjects.
% It uses a state-space random effects (SSRE) model and
% Expectation-Maximization algorithm to estimate learning curves for the
% population and for each subject that characterize the
% dynamics of the learning process as a function of trial number.  For a 
% thorough description of this method of analysis, see:  
%       Smith et al. (2005)  Analysis and Design of Behavioral Experiments
%           to Characterize Population Learning. Journal of
%           Neurophysiology.  93:1776-1792, 2005.
% Throughout this code, references to equations from this journal article 
% will be indicated with a "*":  
% 
% Input: 
%        Responses      (J x N vector) The number ofcorrect responses at
%                       each trial for each subject. J is the number of
%                       subjects. N is number of trials. -- required 
%
%        BackgroundProb     (single value) probabilty of correct by chance
%                                       -- required
%
%        SigE       (single value) SIG_EPSILON, sqrt(variance) of learning
%                               state process -- optional. Default 0.35
% 
%        SigB       (single value) SIG_BETA, sqrt(variance) of SSRE model
%                               -- optional. Default 0.3
%
% other variables
%        J              number of subjects
%        muone          parameter determined by the probability of a correct response by chance
%        Wstart         W(0) (equation A15)*
%        Bstarhat, W    (vectors) unobservable learning state process and its variance (forward estimate)
%        Bstarnew, Wnew (vectors) unobservable learning state process and its variance (backward estimate)
%        number_fail    (vectors) trial number that Newton convergenc failed, if it does
%        covterm2, covterm3 (vectors)  covariance terms needed for the Expectation step of the EM algorithm
%        sigsqEnew      (vectors) SIG_EPSILON^2, estimate of learning state process variance from  the EM algorithm 
%        sigesqBnew     (vectors) SIG_BETA^2, estimate of SSRE model variance
%        Bnew           (vectors) estimate of SSRE model mean
%        A              (vectors) A(k), (equation A17)*
%        betas          (vectors) learning modulation parameter
%        Q              (vectors) the expectation of the complete data log-likelihood of the SSRE model               

% helper functions
%       forwardfilter   solves the forward recursive filtering algorithm to
%                       estimate Bstarhat, W, Bstarold, Wold
%       backwardfilter  solves the backward filter smoothing algorithm to 
%                       estimate Bstarnew, Wnew, and A
%       computecov      estimates the Expectation step's covariance terms
%       em_bino         solves the Maximization step of the EM algorithm to
%                       estimate sigsqEnew, Bnew, and sigsqBnew
%       getQ            computes the expectation of the complete data 
%                       log-likelihood of the SSRE model, Q     

function runanalysis(Responses, BackgroundProb, SigE, SigB) 
 
if nargin < 3
    warning('Start value of individual variances (SigE): default set to 0.3');
    SigE = 0.3; %guess value of starting individual variances 
end
if nargin < 4
    warning('Start value of population variances (SigB): default set to 0.3');
    SigB = 0.3; %guess value of starting population variances 
end

betas =[]; 

muone           = log(BackgroundProb/(1-BackgroundProb)) ;
J = size(Responses,1);
 
 %Set initial values
Bstarstart          = [0 ones(1,size(Responses,1))]';
Wstart   = diag([SigE^2 SigB^2*ones(1,J)]); 
 
number_steps = 1500;
 
%--------------------------------------------------------------------------

%loop through EM algorithm: forward filter, backward filter, compute
%covariance terms, M-step, then get complete data log-likelihood of the
%model
for k = 1:number_steps

     muonenew = muone;

     %Compute the forward (filter algorithm) estimates of the learning states
     %and their variances: BETA*{k|k} and W{k|k}, BETA*{k|k-1} and W{k|k-1} 
     [Bstarhat, W, Bstarold, Wold, number_fail] = ...
                           forwardfilter(Responses, SigE, muonenew, Bstarstart, Wstart);

     if isempty(number_fail)<1
           fprintf(2,'Newton convergence failed at iteration %d  \n', k)
           break
     end

     %Compute the backward (smoothing algorithm) estimates of the learning 
     %state and its variance: BETA*{k|K} and W{k|K}
     [Bstarnew, Wnew, A]  = backwardfilter(Bstarhat, Bstarold, W, Wold);

     %get estimate of x(0) and B(0)'s first
     Bstarnew(:,1)        = Bstarstart;
     Bstarnew(1,1)        = 0.0;


     %Compute the covariance terms R{K|K} and R{k-1,k|K}
     [covterm3, covterm2]     = computecov(Wnew, A, Bstarnew);

     %Compute the EM estimate of the augmented learning state process and its variance  
     [sigsqEnew(k), Bnew(k), sigsqBnew(k)]    = em_bino(covterm3, covterm2, Bstarnew, Wnew, J);

     SigE            = sqrt(sigsqEnew(k));
     SigB            = sqrt(sigsqBnew(k));

     if(SigB<1e-5)
          fprintf(2,'STOPPED BECAUSE SigB got too small')
          break
     end 

     Bstarstart          = Bstarnew(:,1);

     Wstart      = diag([SigE^2 SigB^2*ones(1,J)]);  
     Wnew(:,:,1) = Wstart; 
     Bstarstartsave(k)   = Bstarstart(1);
     betas           = [betas Bstarnew(2:end, 2)];

     %Compute the expectation of the complete data log-likelihood of the
     %augmented learning process model
     Q(k)  = getQ(Bstarnew, Responses, SigE, SigB, Bnew(k), covterm3(1,:), covterm2(1,:), BackgroundProb, Wnew);

     %fprintf(2,'    %d    %f    %f       %f   %f\n', k, SigE, SigB, Bnew(k), Q(k));

    %check for convergence
     CvgceCrit = 1e-8;
     if(k>1)
          diffa = abs(sigsqEnew(k) - sigsqEnew(k-1));
          diffb = abs(Bnew(k)      - Bnew(k-1)); 
          diffc = abs(sigsqBnew(k) - sigsqBnew(k-1));
          diffd = abs(Bstarstartsave(k)- Bstarstartsave(k-1));
          if( diffa < CvgceCrit & diffb < CvgceCrit & diffc < CvgceCrit & diffd < CvgceCrit )
               fprintf(2,'converged in %d iterations \n',k)
               break
          end
     end

end
save('resultspopulation');
%--------------------------------------------------------------------------