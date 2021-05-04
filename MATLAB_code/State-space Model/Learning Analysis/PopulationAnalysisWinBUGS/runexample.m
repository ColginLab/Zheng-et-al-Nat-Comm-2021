%Bayesian implementation of Dynamic Analysis of Learning algorithm in 
%Smith et al., 2007
%run this script within Matlab


I = [0,0,1,0,0,0,1,1,0,1,1,0,1,1,0,0,1,1,1,0,1,1,0,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,0,1,0,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,0,1,0,1,1,1;1,0,0,1,0,1,1,0,1,1,0,0,1,0,1,0,0,0,1,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,1,0,0,1,0,1,0,0,0,1,1,1,1,1,0,1,0,1,0,0,1,0,1,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1,1,1,1,0;1,0,0,0,0,0,0,1,0,0,1,1,0,1,0,1,1,1,0,0,1,0,0,1,0,1,0,1,0,1,1,0,1,0,0,1,1,1,1,1,0,0,0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;];

%put matlab data in format for matbugs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataStruct = struct('n', I, 'T', size(I,2), 'J', size(I,1));

%initial guesses for the MC chains%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3 chains
init1 = struct( 'x', randn(1,size(I, 2) ));
init2 = struct( 'x', randn(1,size(I, 2) ));
init3 = struct( 'x', randn(1,size(I, 2) ));

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
            
                
               
%plot individuals in blue/cyan

figure(1)
pdata =[];
for j = 1:size(I,1)
    pdata = [];
  for t = 1:size(I,2)
        allsamps = [samples.p(1,:,j,t) samples.p(2,:,j,t) samples.p(3,:,j,t)];
        sort_samples = sort(allsamps);
        total        = length(sort_samples);
        ll           = sort_samples(fix(0.05*total));  %lower 95%interval
        ml           = sort_samples(fix(0.5*total));
        ul           = sort_samples(fix(0.95*total));
        pdata = [pdata; t ll ml ul];
        plot(pdata(:,2),'c-'); hold on;
        plot(pdata(:,3),'b-'); hold on;
        plot(pdata(:,4),'c-');
      
  end
end
line([ 0 size(I, 2)+1],[0.5 0.5])
hold on;
xlabel('Trial Number');
ylabel('Pr(Correct Response)')

%plot population in red
pdata = [];
for t = 1:size(I,2)
        allsamps = [samples.pPop(1,:,t) samples.pPop(2,:,t) samples.pPop(3,:,t)];
        sort_samples = sort(allsamps);
        total        = length(sort_samples);
        ll           = sort_samples(fix(0.05*total));  %lower 95%interval
        ml           = sort_samples(fix(0.5*total));
        ul           = sort_samples(fix(0.95*total));
        pdata = [pdata; t ll ml ul];
        h1 = plot(pdata(:,2),'r-'); hold on; set(h1,'LineWidth',2);
        h2 = plot(pdata(:,3),'r.-'); hold on; set(h2,'LineWidth',2);
        h3 = plot(pdata(:,4),'r-'); set(h3,'LineWidth',2);
end
  
title('Red = population; Blue = individuals');

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


