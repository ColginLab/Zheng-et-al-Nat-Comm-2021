%plot the results of the analysis

%load the results of the analysis
load('tmp/results.mat', 'Responses','MaxResponses','samples','BackgroundProb');

%plot Pr(correct response) 
pdata =[];
for t = 1:length(Responses) 
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
plotI(Responses, MaxResponses); hold on;
xlabel('Trial Number');
ylabel('Pr(Correct Response)')
plot(pdata(:,2),'g-'); hold on;
plot(pdata(:,3),'g-'); hold on;
plot(pdata(:,4),'g-');
line([ 0 length(Responses)],[BackgroundProb BackgroundProb]);
