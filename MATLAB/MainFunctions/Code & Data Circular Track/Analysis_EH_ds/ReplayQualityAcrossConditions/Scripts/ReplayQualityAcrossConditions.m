function ReplayQualityAcrossConditions(AnalysisFD)

outdir = strcat(AnalysisFD,'\ReplayQualityAcrossConditions');
if ~isdir(outdir)
    mkdir(outdir)
end

%% Load results
if exist(strcat(AnalysisFD,'\BayesianDecoding_ripple\BayesianDecodingResult.mat'),'file') ~= 2
    error('Please run RunBayesianDecoder_ripple.m first')
else
    post = load(strcat(AnalysisFD,'\BayesianDecoding_ripple\BayesianDecodingResult.mat'));
end

if exist(strcat(AnalysisFD,'\BayesianDecoding_ripple_PreRunTuning\BayesianDecodingResult.mat'),'file') ~= 2
    error('Please run RunBayesianDecoder_ripple_PreRunTuning.m first')
else
    pre = load(strcat(AnalysisFD,'\BayesianDecoding_ripple_PreRunTuning\BayesianDecodingResult.mat'));
end

%% Extract relevant data for plotting
[~,r_squared,~,sessionType] = Extract(post.BayesianDecodingResult);
[~,r_squared_pre,~,sessionType_pre] = Extract(pre.BayesianDecodingResult);

%% Plot results 
% compare different sessions using post-reward position tunning as template
h = figure;
sType = unique(sessionType);
sType = sType([4; 1; 2; 3]); %re-arrange the order
for ss = 1:length(sType)
    in = strcmp(sessionType,sType{ss});
    CI = bootci(5000,@nanmean,r_squared(in));
    CI(1) = nanmean(r_squared(in))-CI(1);
    CI(2) = CI(2)-nanmean(r_squared(in));
    
    plot(ss+randn(size(r_squared(in)))*.05,r_squared(in),'k.')
    hold on
    errorbar(ss,nanmean(r_squared(in)),CI(1),CI(2),'k-','LineWidth',1);
end
xlim([.5 ss+.5])
axis square
ylabel('r2')
set(gca,'XTick',1:ss,'XTickLabel',sType,'Box','off')
saveas(h,strcat(outdir,'\r2dist'),'epsc')
saveas(h,strcat(outdir,'\r2dist'),'png')
close(h)

% compare pre- and post-reward position template
r_squared_all = [r_squared;r_squared_pre];
sessionType_all = [sessionType;sessionType_pre];
templateType = [repmat({'Post-reward'},size(r_squared));...
                repmat({'Pre-reward'},size(r_squared_pre))];
h = figure;
set(h,'OuterPosition',[187.4,361,1200,452.8])
sType = unique(sessionType_all);
sType = sType([4; 1; 2; 3]); %re-arrange the order
ha = zeros(length(sType),1);
for ss = 1:length(sType)
    ha(ss) = subplot(1,length(sType),ss);
    tType = unique(templateType);
    for tt = 1:length(tType)
        in = strcmp(sessionType_all,sType{ss}) & strcmp(templateType,tType{tt});
        CI = bootci(5000,@nanmean,r_squared_all(in));
        CI(1) = nanmean(r_squared_all(in))-CI(1);
        CI(2) = CI(2)-nanmean(r_squared_all(in));
        
        plot(tt+randn(size(r_squared_all(in)))*.05,r_squared_all(in),'k.')
        hold on
        errorbar(tt,nanmean(r_squared_all(in)),CI(1),CI(2),'k-','LineWidth',1);
    end
    xlim([.5 tt+.5])
    axis square
    set(gca,'XTick',1:tt,'XTickLabel',tType,'Box','off')
    title(sType{ss})
    if ss == 1
        ylabel('r2')
    end
end
linkaxes(ha,'xy')
saveas(h,strcat(outdir,'\r2dist_PreVSPost'),'epsc')
saveas(h,strcat(outdir,'\r2dist_PreVSPost'),'png')
close(h)

%% Perform stats
[~,tbl,stats] = kruskalwallis(r_squared,sessionType,'off');
[posthoc] = multcompare(stats,'Display','off');
Stats.kruskalwallis = tbl;
Stats.posthoc = posthoc;

sType = unique(sessionType_all);
sType = sType([4; 1; 2; 3]); %re-arrange the order
pvalues = zeros(length(sType),1);
for ss = 1:length(sType)
    in = strcmp(sessionType_all,sType{ss});
    pre = strcmp(templateType,'Pre-reward');
    pvalues(ss) = ranksum(r_squared_all(in & pre),r_squared_all(in & ~pre));
end
sessionType = sType;
Stats_preVSpost = table(sessionType,pvalues);

%% Export analysis results and all the necessary scripts
FunctionPath = mfilename('fullpath');
FunctionName = mfilename;
OutScript_dir = strcat(outdir,'\Scripts');
if isdir(OutScript_dir)
    rmdir(OutScript_dir,'s')
end
mkdir(OutScript_dir)
[fList] = matlab.codetools.requiredFilesAndProducts(FunctionPath);
for ff = 1:length(fList)
    [~,fname,ext] = fileparts(fList{ff});
    copyfile(fList{ff},strcat(OutScript_dir,'\',fname,ext));
end
save(strcat(outdir,'\Results.mat'),'r_squared_all','sessionType_all','templateType','FunctionName')
save(strcat(outdir,'\Stats.mat'),'Stats','Stats_preVSpost')

% MemoryPerformance = [];
% r2_mean = [];
% r2_median = [];
% ratID = fieldnames(BayesianDecodingResult);
% for rr = 1:length(ratID)
%     dateID = fieldnames(BayesianDecodingResult.(ratID{rr}));
%     for dd = 1:length(dateID)
%         correct_test = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sign_correct_test;
%         correct_posttest = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sign_correct_posttest;
%         r_squared = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared;
%         rmv = strcmp(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sessionType,'PreRun');
%         in = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).NUniqueCell>8 & r_squared>.5;
%         
%         perform = (sum(correct_test==1)+sum(correct_posttest==1))/(length(correct_test)+length(correct_posttest));
%         MemoryPerformance = cat(1,MemoryPerformance,perform);
%         r2_mean = cat(1,r2_mean,nanmean(r_squared(~rmv & in)));
%         r2_median = cat(1,r2_median,nanmedian(r_squared(~rmv & in)));
%     end
% end

function [p_x_n,r_squared,slope,sessionType,Correct_test,stopLoc,rewardLoc] =...
         Extract(BayesianDecodingResult)
p_x_n = {};
r_squared = [];
slope = [];
sessionType = {};
Correct_test = [];
stopLoc = [];
rewardLoc = [];
NUniqueCell_thr = 8;
ratID = fieldnames(BayesianDecodingResult);
for rr = 1:length(ratID)
    dateID = fieldnames(BayesianDecodingResult.(ratID{rr}));
    for dd = 1:length(dateID)
        correct_sample = ~isnan(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sign_correct_sample);
        correct_sample = find(correct_sample);
        NUniqueCell = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).NUniqueCell;
        NTrial = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).NTrial;
        NTrial_post = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).NTrial_post;
        NTrial_preRun = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).NTrial_preRun;
        rewardLoc_session = mode(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).ang_sample_reward_ontrack);
        if ~isempty(NTrial)
            include1 = NTrial~=0;
            include2 = ~isnan(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared) & NUniqueCell > NUniqueCell_thr;
            correct_test = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sign_correct_test;
            in = correct_test==1;
            correct_test = find(in);

            p_x_n = cat(1,p_x_n,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).p_x_n(include1 & include2));
            r_squared = cat(1,r_squared,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared(include1 & include2));
            slope = cat(1,slope,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).slope(include1 & include2));
            sessionType = cat(1,sessionType,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sessionType(include1 & include2));
            Correct_test = cat(1,Correct_test,ismember(NTrial(include1 & include2),correct_test));
            stopLoc = cat(1,stopLoc,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).ang_test_reward_ontrack(NTrial(include1 & include2)));
            rewardLoc = cat(1,rewardLoc,repmat(rewardLoc_session,sum(include1 & include2),1));
        end

        if ~isempty(NTrial_post)
            include3 = NTrial_post~=0 & ~isnan(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared) & NUniqueCell > NUniqueCell_thr;
            correct_posttest = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sign_correct_posttest;
            in = correct_posttest==1;
            correct_posttest = find(in);
            
            p_x_n = cat(1,p_x_n,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).p_x_n(include3));
            r_squared = cat(1,r_squared,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared(include3));
            slope = cat(1,slope,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).slope(include3));
            sessionType = cat(1,sessionType,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sessionType(include3));
            Correct_test = cat(1,Correct_test,ismember(NTrial_post(include3),correct_posttest));
            stopLoc = cat(1,stopLoc,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).ang_posttest_reward_ontrack(NTrial_post(include3)));
            rewardLoc = cat(1,rewardLoc,repmat(rewardLoc_session,sum(include3),1));            
        end

        if ~isempty(NTrial_preRun)
            include4 = NTrial_preRun~=0 & ~isnan(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared) & NUniqueCell > NUniqueCell_thr;
            
            p_x_n = cat(1,p_x_n,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).p_x_n(include4));
            r_squared = cat(1,r_squared,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared(include4));
            slope = cat(1,slope,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).slope(include4));
            sessionType = cat(1,sessionType,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sessionType(include4));
            Correct_test = cat(1,Correct_test,false(sum(include4),1));
            stopLoc = cat(1,stopLoc,NaN(sum(include4),1));
            rewardLoc = cat(1,rewardLoc,repmat(rewardLoc_session,sum(include4),1));
        end        
    end
end