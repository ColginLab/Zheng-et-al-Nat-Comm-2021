%plot the results of the analysis

load('tmp/results.mat', 'Responses', 'samples','BackgroundProb');
%plot individuals in blue/cyan
figure(1)
pdata =[];
for j = 1:size(Responses,1)
    pdata = [];
  for t = 1:size(Responses,2)
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
line([ 0 size(Responses, 2)+1],[0.5 0.5])
hold on;
xlabel('Trial Number');
ylabel('Pr(Correct Response)')

%plot population in red
pdata = [];
for t = 1:size(Responses,2)
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
