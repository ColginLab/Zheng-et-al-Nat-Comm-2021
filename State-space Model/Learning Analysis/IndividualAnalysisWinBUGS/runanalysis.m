function runanalysis(Responses, MaxResponses, BackgroundProb)
%Bayesian implementation of Dynamic Analysis of Learning algorithm in 
%Smith et al., 2004
%run this script within Matlab

%Responses is the data (binary vector)
if nargin<3
    warning('Initial Probability (BackgroundProb) set to 0.25');
    BackgroundProb = 0.25;  %expected initial probability
end

if nargin <2
    warning('Data will be loaded from file data.mat');
    [fid, message]=fopen('data.mat');
    if(fid==-1)
        error('Input parameters needed - Pass them in the function call or save them in the file ''data.mat''.');
    else
		load('data.mat'); %  data is saved as variable Responses and MaxResponses
    end
end

if(length(Responses) ~= length(MaxResponses))
    error('Error - Responses and MaxResponses vectors must be of equal length');
end

if(length(MaxResponses) == sum(MaxResponses == 1))
    %put matlab data in format for matbugs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dataStruct = struct('n', Responses, 'T', length(Responses), 'startp', BackgroundProb);

    %initial guesses for the MC chains%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %3 chains
    init1 = struct( 'tau', 0.5);
    init2 = struct( 'tau', 1);
    init3 = struct( 'tau', 1.5);

    initStructs(1) =  init1;
    initStructs(2) =  init2;
    initStructs(3) =  init3;

    %call Winbugs from in Matlab using matbugs.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matbugs.m is available free from
    %http://www.cs.ubc.ca/~murphyk/Software/MATBUGS/matbugs.html

    [samples, stats, structArray] = matbugs(dataStruct, ...
            fullfile(pwd, 'Model_bern.txt'), ...
            'init', initStructs, ...
                    'nChains', 3, ...
            'view', 0, 'nburnin', 1000, 'nsamples', 5000, ...
            'thin', 10, ...
            'monitorParams', {'p','x','tau','tauprior','sigesq'}, ...
                    'Bugdir', '../../WinBUGS14');
                
else
    %put matlab data in format for matbugs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dataStruct = struct('n', Responses, 'T', length(Responses), 'ntot', MaxResponses, 'startp', BackgroundProb);

    %initial guesses for the MC chains%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %3 chains
    init1 = struct( 'tau', 0.5);
    init2 = struct( 'tau', 1);
    init3 = struct( 'tau', 15);

    initStructs(1) =  init1;
    initStructs(2) =  init2;
    initStructs(3) =  init3;

    [samples, stats, structArray] = matbugs(dataStruct, ...
            fullfile(pwd, 'Model_bin.txt'), ...
            'init', initStructs, ...
                    'nChains', 3, ...
            'view', 0, 'nburnin', 1000, 'nsamples', 5000, ...
            'thin', 10, ...
            'monitorParams', {'p','x','tau','tauprior','sigesq'}, ...
                    'Bugdir', '../../WinBUGS14');
end
save('tmp/results.mat');