function BehaviorPerformance(trackdata,AnalysisFD)

outdir = strcat(AnalysisFD,'\Behaviors');
if ~isdir(outdir)
    mkdir(outdir)
end
%% Extract relevant data
StopProb=[];
TrialDur=[];
StopLocAxis = -3:3;
StopCounts_test = [];
StopCounts_posttest = [];
TrialN = 0;
for tt = 1:length(trackdata)
    temp = load(trackdata{tt});
    ind = strfind(trackdata{tt},'Rat');
    ratID = trackdata{tt}(ind:ind+5);
    if ~isfield(StopProb,ratID)
        StopProb.(ratID).Test = [];
        StopProb.(ratID).posttest = [];
        TrialDur.(ratID).Duration = [];
        TrialDur.(ratID).Label = [];
        sign_correct.(ratID)=[];
        RewardLoc.(ratID)=[];
    end
    reward_ind = temp.Ind_rewardloc;
    RewardLoc.(ratID) = cat(1,RewardLoc.(ratID),reward_ind);
    ind_rewardloc_sample = temp.ind_rewardloc_sample;
    SampleInd = ind_rewardloc_sample-reward_ind;
    sampleDuration = (temp.Ts_sample{1, 1}(:,2)-temp.Ts_sample{1, 1}(:,1))/1e6;
    
    ind_rewardloc_test = temp.ind_rewardloc_test;
    TestInd = ind_rewardloc_test-reward_ind;
    StopProb_test = zeros(length(StopLocAxis),1);
    for ss = 1:length(StopLocAxis)
        if ss == 1
           StopProb_test(ss) = sum(TestInd<=StopLocAxis(ss));
        elseif ss == length(StopLocAxis)
           StopProb_test(ss) = sum(TestInd>=StopLocAxis(ss));
        else
           StopProb_test(ss) = sum(TestInd==StopLocAxis(ss));
        end
    end
    StopCounts_test = cat(1,StopCounts_test,StopProb_test');
    StopProb_test = StopProb_test/sum(StopProb_test);
    testDuration = (temp.Ts_test{1, 1}(:,2)-temp.Ts_test{1, 1}(:,1))/1e6;
    
    ind_rewardloc_posttest = temp.ind_rewardloc_posttest;
    PostTestInd = ind_rewardloc_posttest-reward_ind;
    StopProb_posttest = zeros(length(StopLocAxis),1);
    for ss = 1:length(StopLocAxis)
        if ss == 1
           StopProb_posttest(ss) = sum(PostTestInd<=StopLocAxis(ss));
        elseif ss == length(StopLocAxis)
           StopProb_posttest(ss) = sum(PostTestInd>=StopLocAxis(ss));
        else
           StopProb_posttest(ss) = sum(PostTestInd==StopLocAxis(ss));
        end
    end
    StopCounts_posttest = cat(1,StopCounts_posttest,StopProb_posttest');
    StopProb_posttest = StopProb_posttest/sum(StopProb_posttest);
    posttestDuration = (temp.Ts_posttest{1, 1}(:,2)-temp.Ts_posttest{1, 1}(:,1))/1e6;
    
    Duration = zeros(size(sampleDuration,1)+size(testDuration,1)+size(posttestDuration,1),1);
    Label = cell(size(sampleDuration,1)+size(testDuration,1)+size(posttestDuration,1),1);
    Duration(1:2:2*size(sampleDuration,1)-1,1) = sampleDuration;
    Label(1:2:2*size(sampleDuration,1)-1,1) = {'Sample'};
    Duration(2:2:2*size(sampleDuration,1),1) = testDuration;
    Label(2:2:2*size(sampleDuration,1),1) = {'Test'};
    Duration(2*size(sampleDuration,1)+1:end,1) = posttestDuration;
    Label(2*size(sampleDuration,1)+1:end,1) = {'Posttest'};
    StopProb.(ratID).Test = cat(2,StopProb.(ratID).Test,StopProb_test);
    StopProb.(ratID).posttest = cat(2,StopProb.(ratID).posttest,StopProb_posttest);
    TrialDur.(ratID).Duration = cat(2,TrialDur.(ratID).Duration,Duration);
    TrialDur.(ratID).Label = cat(2,TrialDur.(ratID).Label,Label);
    
    Responses = [temp.sign_correct_test; temp.sign_correct_posttest]';
    %Responses = temp.sign_correct_test(1:3)';
    %Responses(isnan(Responses)) = 0;
    sign_correct.(ratID) = cat(1,sign_correct.(ratID),Responses);
    TrialN = max(TrialN,size(Responses,2));
end

RatID = fieldnames(sign_correct);
sign_correct.all=[];
RewardLoc.all=[];
for rr = 1:length(RatID)
    behaviorOutcome = sign_correct.(RatID{rr});
    temp = NaN(size(behaviorOutcome,1),TrialN);
    temp(:,1:size(behaviorOutcome,2)) = behaviorOutcome;
    sign_correct.all = cat(1,sign_correct.all,temp);
    RewardLoc.all = cat(1,RewardLoc.all,RewardLoc.(RatID{rr}));
end

N = 5000;
RatID = fieldnames(sign_correct);
for rr = 1:length(RatID)
    %PerformanceSample.(RatID{rr}) = runanalysis(sign_correct.(RatID{rr}),
    %1/7);            IS THIS NOT NECESSARY????
    
    %% Bootstrap proportion of correct response
    behaviorOutcome = sign_correct.(RatID{rr});
    proportion = nansum(behaviorOutcome)./sum(~isnan(behaviorOutcome));
    proportion_sample = zeros(N,size(behaviorOutcome,2));
    for ii = 1:N
        sample = randi(size(behaviorOutcome,1),size(behaviorOutcome,1),1);
        behaviorOutcome_sample = behaviorOutcome(sample,:);
        proportion_sample(ii,:) = nansum(behaviorOutcome_sample)./sum(~isnan(behaviorOutcome_sample));
    end
    CI.(RatID{rr}) = [prctile(proportion_sample,2.5); proportion; prctile(proportion_sample,97.5)];
    %% Shuffle trial number to test if behavioral performance improve over time
    proportion_shu = zeros(N,size(behaviorOutcome,2));
    for ii = 1:N
        behaviorOutcome_shuffle = zeros(size(behaviorOutcome));
        for ss = 1:size(behaviorOutcome,1)
            shuffle = randperm(size(behaviorOutcome,2),size(behaviorOutcome,2));
            behaviorOutcome_shuffle(ss,:) = behaviorOutcome(ss,shuffle);
        end    
        proportion_shu(ii,:) = nansum(behaviorOutcome_shuffle)./sum(~isnan(behaviorOutcome_shuffle));
    end
    CI_shu.(RatID{rr}) = [prctile(proportion_shu,2.5); prctile(proportion_shu,50); prctile(proportion_shu,97.5)];
end

%% Plot results
h = figure;
set(h,'OuterPosition',[53,566,1836,507])
RatID = fieldnames(sign_correct);
ha = zeros(length(RatID),1);
for rr = 1:length(RatID)
    if strcmp(RatID{rr},'all')
        ha(length(RatID)) = subplot(1,length(RatID),length(RatID));
        hold on;
        ha1 = plot(CI.(RatID{rr})(2,:)','r','LineWidth',1);
        plot(CI.(RatID{rr})([1 3],:)','r--','LineWidth',1)
        
        ha2 = plot(CI_shu.(RatID{rr})(2,:)','k','LineWidth',1);
        plot(CI_shu.(RatID{rr})([1 3],:)','k--','LineWidth',1)
        titlestr = [RatID{rr},' (n=',num2str(size(sign_correct.(RatID{rr}),1)),')'];
        title(titlestr)
        axis tight square
        ha3 = line(xlim,[1/7 1/7]);
        legend([ha1 ha2 ha3],'Observed','Shuffled','Chance','Location','SouthEast')
    else
        ha(length(RatID)-rr) = subplot(1,length(RatID),length(RatID)-rr);
        hold on;
        plot(CI.(RatID{rr})(2,:)','r','LineWidth',1)
        plot(CI.(RatID{rr})([1 3],:)','r--','LineWidth',1)
        titlestr = [RatID{rr},' (n=',num2str(size(sign_correct.(RatID{rr}),1)),')'];
        title(titlestr)
        axis tight square
        line(xlim,[1/7 1/7]);
    end
end
set(ha,'YLim',[0 1.1])
ylabel(ha(1),'Proportion of correct response')
saveas(h,strcat(outdir,'\Performance'),'epsc')
saveas(h,strcat(outdir,'\Performance'),'png')
close(h)

colorLabel = lines; close(gcf)
h = figure;
set(h,'OuterPosition',[14,465,1895,608])
RatID = fieldnames(TrialDur);
ha = zeros(length(RatID),1);
for rr = 1:length(RatID)
        
    nTrials = size(TrialDur.(RatID{rr}).Label,1);
    Duration = TrialDur.(RatID{rr}).Duration;
    CI = bootci(5000,@nanmean,Duration');
    CI(1,:) = nanmean(Duration,2)'-CI(1,:);
    CI(2,:) = CI(2,:)-nanmean(Duration,2)';
    subplot(1,2,1); hold on
    errorbar(1:nTrials,nanmean(Duration,2),...
        CI(1,:),CI(2,:),'Color',colorLabel(rr,:));
    if rr == length(RatID)
        axis tight
        ylabel('Duration (sec)')
    end
    
    Duration_norm = Duration./repmat(nanmean(Duration),size(Duration,1),1);
    CI = bootci(5000,@nanmean,Duration_norm');
    CI(1,:) = nanmean(Duration_norm,2)'-CI(1,:);
    CI(2,:) = CI(2,:)-nanmean(Duration_norm,2)';
    subplot(1,2,2); hold on
    ha(rr) = errorbar(1:nTrials,nanmean(Duration_norm,2),...
        CI(1,:),CI(2,:),'Color',colorLabel(rr,:));
    if rr == length(RatID)
        axis tight
        ylabel('Normalized duration')
        xlabel('Lap#')
        legend(ha,RatID)
    end
end
saveas(h,strcat(outdir,'\TrialDuration'),'epsc')
saveas(h,strcat(outdir,'\TrialDuration'),'png')
close(h)

%% Plot correct rate as a function of angular position
h = figure;
set(h,'OuterPosition',[53,334.6,1400,381.6])
RatID = fieldnames(sign_correct);
ha = zeros(length(RatID),1);
for rr = 1:length(RatID)
    correctRate = nansum(sign_correct.(RatID{rr}),2)./sum(~isnan(sign_correct.(RatID{rr})),2)*100;
    
    if strcmp(RatID{rr},'all')
        ha(length(RatID)) = subplot(1,length(RatID),length(RatID));
        scatter(RewardLoc.(RatID{rr}),correctRate)
        [r,p] = corr(RewardLoc.(RatID{rr}),correctRate,'Type','Spearman');
        str = ['r=',num2str(r),' p=',num2str(p)];
        titlestr = [RatID{rr},' (n=',num2str(size(sign_correct.(RatID{rr}),1)),') ',str];
        title(titlestr)
        ylabel('% correct')
        axis tight square
    else
        ha(length(RatID)-rr) = subplot(1,length(RatID),length(RatID)-rr);
        scatter(RewardLoc.(RatID{rr}),correctRate)
        titlestr = [RatID{rr},' (n=',num2str(size(sign_correct.(RatID{rr}),1)),')'];
        title(titlestr)
        axis tight square
    end
    xlabel('Reward location index')
end
set(ha,'Box','off')
saveas(h,strcat(outdir,'\CorrectRate_position'),'epsc')
saveas(h,strcat(outdir,'\CorrectRate_position'),'png')
close(h)

%% Plot number of stops relative to reward
h = figure;
set(h,'OuterPosition',[481,73.8,401.6,768.8])
subplot(2,1,1)
StopCounts = sum(StopCounts_test);
CI = bootci(5000,@sum,StopCounts_test);
CI(1,:) = StopCounts-CI(1,:);
CI(2,:) = CI(2,:)-StopCounts;
bar(StopLocAxis,StopCounts); hold on
errorbar(StopLocAxis,StopCounts,CI(1,:),CI(2,:),'k','LineStyle','none')
title('Test')
ylabel('Stop counts')
axis square
set(gca,'Box','off')
subplot(2,1,2)
StopCounts = sum(StopCounts_posttest);
CI = bootci(5000,@sum,StopCounts_posttest);
CI(1,:) = StopCounts-CI(1,:);
CI(2,:) = CI(2,:)-StopCounts;
bar(StopLocAxis,StopCounts); hold on
errorbar(StopLocAxis,StopCounts,CI(1,:),CI(2,:),'k','LineStyle','none')
title('Post-test')
ylabel('Stop counts')
xlabel('Reward location #')
axis square
set(gca,'Box','off')
saveas(h,strcat(outdir,'\Stop_position'),'epsc')
saveas(h,strcat(outdir,'\Stop_position'),'png')
close(h)

%% Export analysis results and all the necessary scripts
FunctionPath = mfilename('fullpath');
FunctionName = mfilename;
% % %
% OutScript_dir = strcat(outdir,'\Scripts');
% if isdir(OutScript_dir)
%     rmdir(OutScript_dir,'s')
% end
% mkdir(OutScript_dir)
% [fList] = matlab.codetools.requiredFilesAndProducts(FunctionPath);          % LOOKS UNNECESSARY 
% for ff = 1:length(fList)
%     [~,fname,ext] = fileparts(fList{ff});
%     copyfile(fList{ff},strcat(OutScript_dir,'\',fname,ext));
% end
% % %
save(strcat(outdir,'\Results.mat'),'StopLocAxis','StopCounts','FunctionName')

function samples = runanalysis(Responses, BackgroundProb)
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
% Directory with the model text file for Winbugs
model_dir = 'C:\CA1Replayproject\Matlab\LearningAnalysis\PopulationAnalysisWinBUGS';
%call Winbugs from in matlab%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[samples, stats, structArray] = matbugs(dataStruct, ...
        fullfile(model_dir, 'Model_bern.txt'), ...
        'init', initStructs, ...
                'nChains', 3, ...
        'view', 0, 'nburnin', 1000, 'nsamples', 2000, ...
        'thin', 10, 'workingDir',[model_dir,'\tmp'],...
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
