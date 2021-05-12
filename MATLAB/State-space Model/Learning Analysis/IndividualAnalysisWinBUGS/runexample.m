%Bayesian implementation of Dynamic Analysis of Learning algorithm in 
%Smith et al., 2004
%run this script within Matlab

%Ipaper is the data (binary vector)
clear

I = [zeros(1,10) 1 1 0 0 0 0 1 0 0 ones(1,11)];

startp = 0.25;  %expected initial probability

%put matlab data in format for matbugs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataStruct = struct('n', I, 'T', length(I), 'startp', startp);

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
                'Bugdir', 'I:/MATLAB/State-space Model/WinBUGS14');
                            
%plot Pr(correct response) (Figure 1 in Smith, Frank, et al. 2004)%%%%%%%%%

pdata =[];
for t = 1:length(I) 
        allsamples   = [samples.p(1,:,t) samples.p(2,:,t) samples.p(3,:,t)];
        sort_samples = sort(allsamples);
        total        = length(sort_samples);
        ll           = sort_samples(fix(0.05*total));  %lower 95%interval
        ml           = sort_samples(fix(0.5*total));
        ul           = sort_samples(fix(0.95*total));
        pdata = [pdata; t ll ml ul];
end

figure(1);
subplot(211)
plotI(I, ones(1,length(I))); hold on;
xlabel('Trial Number');
ylabel('Pr(Correct Response)')
plot(pdata(:,2),'g-'); hold on;
plot(pdata(:,3),'g-'); hold on;
plot(pdata(:,4),'g-');
line([ 0 length(I)],[startp startp])

