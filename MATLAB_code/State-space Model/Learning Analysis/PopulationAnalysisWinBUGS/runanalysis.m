function runanalysis(Responses, BackgroundProb)
%Bayesian implementation of Dynamic Analysis of Learning algorithm in 
%Smith et al., 2007
%run this script within Matlab
%Responses is the data (binary vector)

if nargin<2
    warning('Initial Probability (BackgroundProb) set to 0.25');
    BackgroundProb = 0.25;  %expected initial probability
end

if nargin <1
    warning('Data will be loaded from file data.mat');
    [fid, message]=fopen('data.mat');
    if(fid==-1)
        error('Input parameters needed - Pass them in the function call or save them in the file ''data.mat''.');
    else
        load('data.mat'); %  data is saved as variable Responses
    end
end

%put matlab data in format for matbugs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataStruct = struct('n', Responses, 'T', size(Responses,2), 'J', size(Responses,1), 'startp', BackgroundProb);

%initial guesses for the MC chains%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3 chains
init1 = struct( 'x', randn(1,size(Responses, 2) ));
init2 = struct( 'x', randn(1,size(Responses, 2) ));
init3 = struct( 'x', randn(1,size(Responses, 2) ));

initStructs(1) =  init1;
initStructs(2) =  init2;
initStructs(3) =  init3;


%call Winbugs from in matlab%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[samples, stats, structArray] = matbugs(dataStruct, ...
        fullfile(pwd, 'Model_bern.txt'), ...
        'init', initStructs, ...
                'nChains', 3, ...
        'view', 1, 'nburnin', 1000, 'nsamples', 2000, ...
        'thin', 10, ...
        'monitorParams', {'p','x','pPop','sigesq','tauprior','beta','beta0'}, ...
                'Bugdir', 'C:/Program Files/WinBUGS14');
              
TOOBIG1  = find(stats.Rhat.beta > 1.2);
TOOBIG2  = find(stats.Rhat.x > 1.2);
TOOBIG3  = find(stats.Rhat.beta0> 1.2);


fprintf('Checking for MC convergence using Rhat which should be less than 1.2. \n Rhats above 1.2 are shown below.  \n')

if(~isempty(TOOBIG1)) 
     fprintf('WARNING: Monte Carlo convergence for beta is not great.\n')
     fprintf('Largest value of x is %f \n', max(stats.Rhat.beta(TOOBIG1)) )
end

if(~isempty(TOOBIG2)) 
     fprintf('WARNING: Monte Carlo convergence for x is poor.\n')
     fprintf('Largest value of p is %f \n', max(stats.Rhat.x(TOOBIG2) ) )
end

if(~isempty(TOOBIG3) )
     fprintf('WARNING: Monte Carlo convergence for beta0 is not great.\n')
     fprintf('Largest value of xb is %f \n', max(stats.Rhat.beta0(TOOBIG3)) )
end

save('tmp/results.mat');
