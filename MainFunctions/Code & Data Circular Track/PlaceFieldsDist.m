function PlaceFieldsDist(trackdata,Indir,Outdir)

%outdir = strcat(Outdir,'\PlaceFieldsDist_alignedToReward');
outdir = Outdir;
if ~isdir(outdir)
    mkdir(outdir)
end

if exist(strcat(Indir,'\LapFRcorr.mat'),'file') ~= 2
    error('Please run RunPlaceCellPropertiesAcrossLap_v6.m first')
else
    load(strcat(Indir,'\LapFRcorr.mat'))
end

PosTuning=[];
PosTuning_preRun=[];
RatID = fieldnames(LapFRcorr);
SessionID = [];
sessionCount = 0;
for rr = 1:length(RatID)
    DateID = fieldnames(LapFRcorr.(RatID{rr}));
    for ii = 1:length(DateID)
        sessionCount = sessionCount+1;
        if strcmp(DateID{ii}(end),'T')
            trackfileName = strcat(DateID{ii}(2:end-2),'_CT_tracking.mat');
        else
            trackfileName = strcat(DateID{ii}(2:end-3),'_CT_tracking_',DateID{ii}(end),'.mat');
        end
        ind = find(~isempty(strfind(trackdata,trackfileName)));
        trackInfo = load(trackdata{ind});
        posbinaxis = LapFRcorr.(RatID{rr}).(DateID{ii}).Test.posbinaxis;
        rewardLoc = LapFRcorr.(RatID{rr}).(DateID{ii}).Test.rewardLoc/180*pi;
        restBoxLoc = [0, trackInfo.Ang_StartZone_depart_ontrack,...
                      trackInfo.Ang_StartZone_arrive_ontrack, 2*pi];
        FR_SampleTestRun = nanmean(cat(3,LapFRcorr.(RatID{rr}).(DateID{ii}).Sample.lapFR,...
                                  LapFRcorr.(RatID{rr}).(DateID{ii}).Test.lapFR,...
                                  LapFRcorr.(RatID{rr}).(DateID{ii}).Postrunning.lapFR),3);
        FR_PreRun = nanmean(LapFRcorr.(RatID{rr}).(DateID{ii}).Prerunning.lapFR,3);
        % Remove firing rate within rest box 
%         restBoxLoc_ind = match(restBoxLoc,posbinaxis/180*pi);
%         FR_SampleTestRun(:,restBoxLoc_ind(1):restBoxLoc_ind(2))=NaN;
%         FR_SampleTestRun(:,restBoxLoc_ind(3):restBoxLoc_ind(4))=NaN;
%         FR_PreRun(:,restBoxLoc_ind(1):restBoxLoc_ind(2))=NaN;
%         FR_PreRun(:,restBoxLoc_ind(3):restBoxLoc_ind(4))=NaN;
        
        % Align reward position to middle
        ind = match(rewardLoc,posbinaxis/180*pi);
        ind_half = match(180,posbinaxis);
        shift = ind-ind_half;
        FR_SampleTestRun = circshift(FR_SampleTestRun,-shift,2);
        FR_PreRun = circshift(FR_PreRun,-shift,2);

        posaxis = posbinaxis-180;
        PosTuning = cat(1,PosTuning,FR_SampleTestRun);
        PosTuning_preRun = cat(1,PosTuning_preRun,FR_PreRun);
        SessionID = cat(1,SessionID,repmat(sessionCount,size(FR_SampleTestRun,1),1));
    end
end

% Sort cell by positions of maximal firing rate
[peakFR,maxInd] = nanmax(PosTuning,[],2);
in = peakFR>=1;
maxInd = maxInd(in);
PosTuning = PosTuning(in,:);
PosTuning_preRun = PosTuning_preRun(in,:);
SessionID = SessionID(in);
[~,sortInd] = sort(maxInd);
PosTuning_after = PosTuning(sortInd,:);
PosTuning_preRun_after  = PosTuning_preRun(sortInd,:);

% Get place field width before and after reward presentation
[FieldWidth,PeakFRLoc] = FindWidth(PosTuning,PosTuning_preRun,posaxis);

h = Makeplots(PosTuning_after,PosTuning_preRun_after,posaxis/180*pi);
saveas(h,strcat(outdir,'\PlaceFieldDis_after'),'epsc')
saveas(h,strcat(outdir,'\PlaceFieldDis_after'),'fig')
close(h)

%h = MakeRewardPerdictionPlots(PosTuning_preRun_after,maxInd_sort,posaxis,[-20 40]);
h = MakeRewardPerdictionPlots_v2(PosTuning,PosTuning_preRun,[-4 4],posaxis,1);
% saveas(h,strcat(outdir,'\PlaceFieldDisFromFutureReward_after'),'epsc')
% saveas(h,strcat(outdir,'\PlaceFieldDisFromFutureReward_after'),'png')
close(h)

% Sort cell by positions of maximal firing rate
[~,maxInd] = nanmax(PosTuning_preRun,[],2);
[~,sortInd] = sort(maxInd);
PosTuning_before = PosTuning(sortInd,:);
PosTuning_preRun_beofre  = PosTuning_preRun(sortInd,:);

h2 = Makeplots(PosTuning_before,PosTuning_preRun_beofre,posaxis/180*pi);
% saveas(h2,strcat(outdir,'\PlaceFieldDis_before'),'epsc')
% saveas(h2,strcat(outdir,'\PlaceFieldDis_before'),'png')
close(h2)

% Make correlation plot
[h3, rho] = MakeCorrelationPlot(PosTuning,PosTuning_preRun,posaxis,SessionID);
% saveas(h3,strcat(outdir,'\PlaceFieldDisCorr'),'epsc')
% saveas(h3,strcat(outdir,'\PlaceFieldDisCorr'),'png')
close(h3)

% Make spatial correlation against center of mass of a place cell
[h4,COM] = MakeCorrelationPlot_v2(PosTuning,PosTuning_preRun,posaxis);
% saveas(h4,strcat(outdir,'\SpatialCorrBetweenRewardPresentation'),'epsc')
% saveas(h4,strcat(outdir,'\SpatialCorrBetweenRewardPresentation'),'png')
close(h4)


%% Export analysis results and all the necessary scripts
% FunctionPath = mfilename('fullpath');
% FunctionName = mfilename;
% OutScript_dir = strcat(outdir,'\Scripts');
% if isdir(OutScript_dir)
%     rmdir(OutScript_dir,'s')
% end
% mkdir(OutScript_dir)
% [fList] = matlab.codetools.requiredFilesAndProducts(FunctionPath);
% for ff = 1:length(fList)
%     [~,fname,ext] = fileparts(fList{ff});
%     copyfile(fList{ff},strcat(OutScript_dir,'\',fname,ext));
% end
save(strcat(outdir,'\PlotOutput.mat'),'PosTuning_preRun','PosTuning','posaxis','rho')

%% Helper function
function [PosTuning,CellID,RippleFR,posaxis] = ExtractPosTuning(input)
PosTuning = [];
CellID = [];
RippleFR = [];
ratID = fieldnames(input.BayesianDecodingResult);
posbinsize = input.BayesianDecodingResult.Rat139.D20170220CT2.param.posbinsize/2/pi*360;
posaxis = 0:posbinsize:360-posbinsize;
for rr = 1:length(ratID)
    dateID = fieldnames(input.BayesianDecodingResult.(ratID{rr}));
    for dd = 1:length(dateID)
        posTuning = input.BayesianDecodingResult.(ratID{rr}).(dateID{dd}).PosTuning;
        cellID = input.BayesianDecodingResult.(ratID{rr}).(dateID{dd}).cellID;
        rewardPos = mode(input.BayesianDecodingResult.(ratID{rr}).(dateID{dd}).ang_sample_reward_ontrack);
        shift = match(rewardPos/2/pi*360,posaxis)-length(posaxis)/2;
        PosTuning = cat(1,PosTuning,circshift(posTuning,[0,-shift]));
        CellID = cat(1,CellID,strcat(ratID{rr},dateID{dd},cellID));
        RippleFR = cat(1,RippleFR,input.BayesianDecodingResult.(ratID{rr}).(dateID{dd}).RippleFR);
    end
end
posaxis = posaxis-180;

function h = Makeplots(PosTuning,PosTuning_preRun,posaxis)
FRthr = 0.5; % in Hz
% Normalize position tuning
PosTuning_norm = PosTuning./repmat(max(PosTuning,[],2),1,size(PosTuning,2));
PosTuning_preRun_norm = PosTuning_preRun./repmat(max(PosTuning_preRun,[],2),1,size(PosTuning_preRun,2));

% Make plots
h = figure('Units','normalized','Position',[0 0 .5 1]);
h1 = subplot(4,2,[1 3 5]);
imagesc(posaxis,1:length(PosTuning_preRun_norm),PosTuning_preRun_norm)
hold on
%plot([0,0],[1 length(PosTuning_preRun_norm)],'LineStyle','--','LineWidth',1,'Color','w')
colormap jet
title('Before reward presentation')
ylabel('CA1 cell indices')
axis tight xy
colorbar
set(gca,'XTick',[])

h2 = subplot(4,2,[2 4 6]);
imagesc(posaxis,1:length(PosTuning_norm),PosTuning_norm);
hold on
%plot([0,0],[1 length(PosTuning_preRun_norm)],'LineStyle','--','LineWidth',1,'Color','w')
title('After reward presentation')
axis tight xy
colormap jet
hb = colorbar;
ylabel(hb,'Normalized firing rate')
set(gca,'YTick',[],'XTick',[])
set(hb,'YTick',[0 1],'YTickLabel',{'min','max'})
% plot proportion 
proportion_preRun = sum(PosTuning_preRun>FRthr)./sum(~isnan(PosTuning_preRun));
h3 = subplot(4,2,7);
ha = plot(posaxis,proportion_preRun,'LineWidth',1);
%set(gca,'YScale','log')
axis tight
set(h3,'box','off');
set(ha,'LineWidth',1);
xlabel('Aligned angular position (rad)')
ylabel('Cell proportion')
colorbar
linkaxes([h1 h3],'x')

proportion = sum(PosTuning>FRthr)./sum(~isnan(PosTuning));
h4 = subplot(4,2,8);
ha = plot(posaxis,proportion,'LineWidth',1);
set(gca,'YTick',[])
axis tight
set(h4,'box','off');
set(ha,'LineWidth',1);
xlabel('Aligned angular position (rad)')
colorbar
linkaxes([h2 h4],'x')
linkaxes([h3 h4],'xy')
set(h3,'YLim',[0 .35])

function [FieldWidth,PeakFRLoc] = FindWidth(PosTuning,PosTuning_preRun,posaxis)
%% Compare the place field width before and after reward presentation within a cell
FRthr = 1; % in Hz
% align position tunning to the middle
[~,maxInd] = max(PosTuning,[],2);
shift = maxInd-ceil(size(PosTuning,2)/2);
PosTuning_s = PosTuning;
width = [];
peakind = [];
for cellid = 1:size(PosTuning,1)
    PosTuning_s(cellid,:) = circshift(PosTuning_s(cellid,:),-shift(cellid),2);    
    %xmax = findpeaks(PosTuning_s(cellid,:),FRthr);
    xmax.loc = maxInd(cellid)-shift(cellid);
    Bound = NaN(length(xmax.loc),2);
    for ii = 1:length(xmax.loc)
        bound_ind = [xmax.loc(ii) xmax.loc(ii)];
        while PosTuning_s(cellid,bound_ind(1)) >= FRthr && bound_ind(1)>1
            bound_ind(1) = bound_ind(1)-1;
        end
        while PosTuning_s(cellid,bound_ind(2)) >= FRthr && bound_ind(2)<size(PosTuning_s,2)-1
            bound_ind(2) = bound_ind(2)+1;
        end
        Bound(ii,:) = bound_ind;
    end
    if length(xmax.loc) > 1
        rmv = sum(diff(Bound)==0,2)==2;
        Bound(rmv,:)=[];
    end
    width = cat(1,width,Bound(:,2)-Bound(:,1));
    peakind = cat(1,peakind,repmat(maxInd(cellid),size(Bound,1),1));
end

% align position tunning to the middle
% [~,maxInd] = max(PosTuning_preRun,[],2);
% shift = maxInd-ceil(size(PosTuning_preRun,2)/2);
PosTuning_preRun_s = PosTuning_preRun;
width_pre = [];
peakind_pre = [];
for cellid = 1:size(PosTuning,1)
    PosTuning_preRun_s(cellid,:) = circshift(PosTuning_preRun_s(cellid,:),-shift(cellid),2);    
    %xmax = findpeaks(PosTuning_preRun_s(cellid,:),FRthr);
    xmax.loc = maxInd(cellid)-shift(cellid);
    Bound = NaN(length(xmax.loc),2);
    for ii = 1:length(xmax.loc)
        bound_ind = [xmax.loc(ii) xmax.loc(ii)];
        while PosTuning_preRun_s(cellid,bound_ind(1)) >= FRthr && bound_ind(1)>1
            bound_ind(1) = bound_ind(1)-1;
        end
        while PosTuning_preRun_s(cellid,bound_ind(2)) >= FRthr && bound_ind(2)<size(PosTuning_preRun_s,2)-1
            bound_ind(2) = bound_ind(2)+1;
        end
        Bound(ii,:) = bound_ind;
    end
    if length(xmax.loc) > 1
        rmv = sum(diff(Bound)==0,2)==2;
        Bound(rmv,:)=[];
    end
    width_pre = cat(1,width_pre,Bound(:,2)-Bound(:,1));
    peakind_pre = cat(1,peakind_pre,repmat(maxInd(cellid),size(Bound,1),1));
end
FieldWidth = [width_pre width];
PeakFRLoc = posaxis(peakind);
    
function h = MakeRewardPerdictionPlots(PosTuning_preRun,maxPosInd,posaxis,posRange)
% Normalize position tuning
PosTuning_preRun_norm = PosTuning_preRun./repmat(max(PosTuning_preRun,[],2),1,size(PosTuning_preRun,2));

posRangeInd = match(posRange,posaxis);
in = maxPosInd>=min(posRangeInd) & maxPosInd<=max(posRangeInd);
PosTuningMat = PosTuning_preRun_norm(in,:);

% Sort position tuning matrix
% [~,maxInd] = max(PosTuningMat,[],2);
% [~,sortInd] = sort(maxInd);
% PosTuningMat = PosTuningMat(sortInd,:);

% Obtain null distribution of place fields
nShuffle = 5000;
NullDist = zeros(nShuffle,size(PosTuning_preRun_norm,2));
for nn = 1:nShuffle
    shu = randperm(size(PosTuning_preRun_norm,1),sum(in));
    NullDist(nn,:) = median(PosTuning_preRun_norm(shu,:));
end
error = [prctile(NullDist,97.5)-median(NullDist); median(NullDist)-prctile(NullDist,2.5)];

h = figure;
set(h,'OuterPosition',[615,49,858,1032])
h1 = subplot(4,1,[1 2 3]); hold on
imagesc(posaxis,1:size(PosTuningMat,1),PosTuningMat)
title('Place cells near future reward locations')
ylabel('CA1 cell indices')
axis tight xy
pos = get(h1,'Position');
hb = colorbar;
ylabel(hb,'Normalized firing rate')
set(h1,'Position',pos)

h2 = subplot(4,1,4); hold on
ha1 = shadedErrorBar(posaxis,median(NullDist),error,'k');
ha = plot(posaxis,median(PosTuningMat));
set(gca,'YScale','log')
set([ha1.mainLine ha],'LineWidth',1);
legend([ha1.mainLine ha],'Null','Observe')
xlabel('Aligned angular position (deg)')
ylabel('Median normalized FR')
linkaxes([h1 h2],'x')

function h = MakeRewardPerdictionPlots_v2(PosTuning,PosTuning_preRun,rewardPosRange,posaxis,FRthreshold)
% Normalize position tuning
PosTuning_norm = PosTuning./repmat(max(PosTuning,[],2),1,size(PosTuning,2));
PosTuning_preRun_norm = PosTuning_preRun./repmat(max(PosTuning_preRun,[],2),1,size(PosTuning_preRun,2));

rewardPosRangeInd = match(rewardPosRange,posaxis);
range = min(rewardPosRangeInd):max(rewardPosRangeInd);
in = sum(PosTuning(:,range)>FRthreshold & PosTuning_preRun(:,range)<FRthreshold,2) > 0; %Select place cells with new fields appeared at reward location after reward presentation
PosTuning_subset = PosTuning_norm(in,:);
PosTuning_subset_preRun = PosTuning_preRun_norm(in,:);

% Sort position tuning matrix
[~,maxInd] = max(PosTuning_subset,[],2);
[~,sortInd] = sort(maxInd);
PosTuning_subset = PosTuning_subset(sortInd,:);
PosTuning_subset_preRun = PosTuning_subset_preRun(sortInd,:);

% Obtain null distribution of place fields
nShuffle = 5000;
NullDist = zeros(nShuffle,size(PosTuning,2));
NullDist_preRun = zeros(nShuffle,size(PosTuning_preRun,2));
for nn = 1:nShuffle
    shu = randperm(size(PosTuning_preRun,1),sum(in));
    NullDist(nn,:) = median(PosTuning_norm(shu,:));
    NullDist_preRun(nn,:) = median(PosTuning_preRun_norm(shu,:));
end
error = [prctile(NullDist,97.5)-median(NullDist); median(NullDist)-prctile(NullDist,2.5)];
error_preRun = [prctile(NullDist_preRun,97.5)-median(NullDist_preRun); median(NullDist_preRun)-prctile(NullDist_preRun,2.5)];

h = figure;
set(h,'OuterPosition',[345,128,1041,891])
h1 = subplot(2,2,1); hold on
imagesc(posaxis,1:size(PosTuning_subset_preRun,1),PosTuning_subset_preRun)
title('Before reward presentation')
ylabel('CA1 cell indices')
axis tight xy square
pos = get(h1,'Position');
hb = colorbar;
ylabel(hb,'Normalized firing rate')
set(h1,'Position',pos)

h2 = subplot(2,2,2); hold on
imagesc(posaxis,1:size(PosTuning_subset,1),PosTuning_subset)
title('After reward presentation')
ylabel('CA1 cell indices')
axis tight xy square
pos = get(h2,'Position');
hb = colorbar;
ylabel(hb,'Normalized firing rate')
set(h2,'Position',pos)

h3 = subplot(2,2,3); hold on
ha1 = shadedErrorBar(posaxis,median(NullDist_preRun),error_preRun,'k');
ha = plot(posaxis,median(PosTuning_subset_preRun));
axis square
set(gca,'YScale','log')
set([ha1.mainLine ha],'LineWidth',1);
xlabel('Aligned angular position (deg)')
ylabel('Median normalized FR')
linkaxes([h1 h3],'x')

h4 = subplot(2,2,4); hold on
ha1 = shadedErrorBar(posaxis,median(NullDist),error,'k');
ha = plot(posaxis,median(PosTuning_subset));
axis square
set(gca,'YScale','log')
set([ha1.mainLine ha],'LineWidth',1);
legend([ha1.mainLine ha],'Null','Observe')
xlabel('Aligned angular position (deg)')
ylabel('Median normalized FR')
linkaxes([h2 h4],'x')

function [h, rho] = MakeCorrelationPlot(PosTuning,PosTuning_preRun,posaxis,SessionID)
% session = unique(SessionID);
% corrData = zeros(length(session),length(posaxis));
% for ss = 1:length(session)
%     in = SessionID==session(ss);
%     rho = corr(PosTuning_preRun(in,:),PosTuning(in,:));
%     corrData(ss,:) = diag(rho);
% end

rho = corr(PosTuning_preRun,PosTuning);    
h = figure; hold on
imagesc(posaxis,posaxis,rho)
axis square tight xy
plot([0 0],ylim,'w--')
plot(xlim,[0 0],'w--')
ylabel('Before reward position')
xlabel('After reward position')
set(gca,'CLim',[0 1])
hb = colorbar;
ylabel(hb,'Pearson r')

function [h,com,rho] = MakeCorrelationPlot_v2(PosTuning,PosTuning_preRun,posaxis)
com = posaxis*PosTuning'./sum(PosTuning,2)';
rho = zeros(size(PosTuning,1),1);
for ii = 1:size(PosTuning,1)
    rho(ii) = corr(PosTuning(ii,:)',PosTuning_preRun(ii,:)');
end

h = figure; hold on
scatter(abs(com),rho,'.')
axis square tight
plot([0 0],ylim,'w--')
ylabel('Spatial correlation')
xlabel('COM location from reward position')

function [h] = MakeRippleFRPlot(RippleFR,COM)
h = figure;
set(h,'OuterPosition',[163,451,1632,622])
ha1 = subplot(1,3,1);
scatter(abs(COM),RippleFR(:,4),'.')
ylabel('Ripple FR before reward (Hz)')
xlabel('COM location from reward position')
axis square tight

ha2 = subplot(1,3,2);
scatter(abs(COM),mean(RippleFR(:,1:3),2),'.')
ylabel('Ripple FR after reward (Hz)')
xlabel('COM location from reward position')
axis square tight

ha3 = subplot(1,3,3);
scatter(abs(COM),mean(RippleFR(:,1:3),2)-RippleFR(:,4),'.')
ylabel('after reward - before reward (Hz)')
xlabel('COM location from reward position')
axis square tight

set([ha1 ha2],'YScale','log')