%% Plot correlation matrix between correct and error trials (include all cells)
% This version make pairs of trials into variables
function CorrErrPopVecCorr_v7(trackdata,AnalysisFD)

if exist(strcat(AnalysisFD,'\LapFRcorr\LapFRcorr.mat'),'file') ~= 2
    error('Please run RunPlaceCellPropertiesAcrossLap_v4.m first')
else
    load(strcat(AnalysisFD,'\LapFRcorr\LapFRcorr.mat'))
end

outdir = strcat(AnalysisFD,'\CorrErrPopVecCorr_allcells');
if ~isdir(outdir)
    mkdir(outdir)
end

nLapThr = 3;
RatID = fieldnames(LapFRcorr);
PairID = {'Sample',     'TestCorrect','SampleCorrect','SampleError';...
          'TestCorrect','TestCorrect','TestCorrect'  ,'TestCorrect';...
          'TestError'  ,'TestError'    ,'TestError'  ,'TestError'}; 
for pp = 1:size(PairID,2)
    outdir_sub = strcat(outdir,'\',PairID{1,pp},'_',PairID{2,pp},'_',PairID{3,pp});
    if ~isdir(outdir_sub)
        mkdir(outdir_sub)
    end      
    FRmat = [];
    FRmat_norm = [];
    FRmat_pre = [];
    CellType = {};
    SessionID = {};
    LapSpkind_preRun = [];
    LapSpkind_sample = [];
    LapSpkind_test = [];
    LapSpkind_correctInd = [];
    for rr = 1:length(RatID)
        DateID = fieldnames(LapFRcorr.(RatID{rr}));
        for ii = 1:length(DateID)

            if strcmp(DateID{ii}(end),'T')
                trackfileName = strcat(DateID{ii}(2:end-2),'_CT_tracking.mat');
            else
                trackfileName = strcat(DateID{ii}(2:end-3),'_CT_tracking_',DateID{ii}(end),'.mat');
            end
            ind = contains(trackdata,trackfileName);
            trackInfo = load(trackdata{ind});
            % Find reward range      
            rewardLoc = LapFRcorr.(RatID{rr}).(DateID{ii}).Test.rewardLoc/180*pi;
            posbinaxis = LapFRcorr.(RatID{rr}).(DateID{ii}).Test.posbinaxis;
            ind = match(rewardLoc,trackInfo.Ang_RewardLoc_ontrack);
            if ind == 1
                rewardPos_range = [rewardLoc-(trackInfo.Ang_RewardLoc_ontrack(ind+1)-rewardLoc)/2 (trackInfo.Ang_RewardLoc_ontrack(ind+1)+rewardLoc)/2];
            elseif ind == length(trackInfo.Ang_RewardLoc_ontrack)
                rewardPos_range = [(trackInfo.Ang_RewardLoc_ontrack(ind-1)+rewardLoc)/2 rewardLoc+(rewardLoc-trackInfo.Ang_RewardLoc_ontrack(ind-1))/2];
            else
                rewardPos_range = [(trackInfo.Ang_RewardLoc_ontrack(ind-1)+rewardLoc)/2 (trackInfo.Ang_RewardLoc_ontrack(ind+1)+rewardLoc)/2];
            end
            FR_SampleTestRun = nanmean(cat(3,LapFRcorr.(RatID{rr}).(DateID{ii}).Sample.lapFR,...
                                             LapFRcorr.(RatID{rr}).(DateID{ii}).Test.lapFR,...
                                             LapFRcorr.(RatID{rr}).(DateID{ii}).Postrunning.lapFR),3);
            FR_PreRun = nanmean(LapFRcorr.(RatID{rr}).(DateID{ii}).Prerunning.lapFR,3);
            % Select cells with place fields near reward position
            ind_reward = match(rewardPos_range,posbinaxis/180*pi);
    %         [peakFR,maxInd] = max(FR_PreRun,[],2);
    %         ActiveNearReward_preRun = peakFR >= 1 & maxInd>=ind_reward(1) & maxInd<=ind_reward(2);
            ActiveNearReward_preRun = mean(FR_PreRun(:,ind_reward(1):ind_reward(2)),2) >= 1;
    %         [peakFR,maxInd] = max(FR_SampleTestRun,[],2);
    %         ActiveNearReward_SampleTestRun = peakFR >= 1 & maxInd>=ind_reward(1) & maxInd<=ind_reward(2);
            ActiveNearReward_SampleTestRun = mean(FR_SampleTestRun(:,ind_reward(1):ind_reward(2)),2) >= 1;
            celltype = cell(size(ActiveNearReward_preRun));
            celltype(ActiveNearReward_preRun) = {'RewardPre'};
            celltype(~ActiveNearReward_preRun & ActiveNearReward_SampleTestRun) = {'RewardEmerged'};

            lapFR_sample = LapFRcorr.(RatID{rr}).(DateID{ii}).Sample.lapFR;
            lapFR_test = LapFRcorr.(RatID{rr}).(DateID{ii}).Test.lapFR;
            lapSpkind_preRun = LapFRcorr.(RatID{rr}).(DateID{ii}).Prerunning.lapSpkind;
            lapSpkind_sample = LapFRcorr.(RatID{rr}).(DateID{ii}).Sample.lapSpkind;
            lapSpkind_test = LapFRcorr.(RatID{rr}).(DateID{ii}).Test.lapSpkind;
            lapSpkind_preRun = cellfun(@circDist,lapSpkind_preRun,repmat({rewardLoc},size(lapSpkind_preRun)),'UniformOutput',false);
            lapSpkind_sample = cellfun(@circDist,lapSpkind_sample,repmat({rewardLoc},size(lapSpkind_sample)),'UniformOutput',false);
            lapSpkind_test = cellfun(@circDist,lapSpkind_test,repmat({rewardLoc},size(lapSpkind_test)),'UniformOutput',false);
            corrInd = trackInfo.sign_correct_test;
            stopLoc = trackInfo.ang_test_reward_ontrack;
            % Align reward position to middle
            ind = match(rewardLoc,posbinaxis/180*pi);
            ind_half = match(180,posbinaxis);
            shift = ind-ind_half;
            lapFR_sample = circshift(lapFR_sample,-shift,2);
            lapFR_test = circshift(lapFR_test,-shift,2);
            FR_PreRun = circshift(FR_PreRun,-shift,2);

            nCorrect = sum(corrInd==1);
            nError = sum(corrInd==0);
            if nCorrect >= nLapThr && nError >= nLapThr %some cells during the session with insufficient # of correct and error trials are excluded
                LapFR = [];
                for ff = 1:size(PairID,1) % loop through each condition
                    if contains(PairID{ff,pp},'Sample') 
                        if contains(PairID{ff,pp},'Correct')
                            lapFR = squeeze(mean(lapFR_sample(:,:,corrInd==1),3));
                        elseif contains(PairID{ff,pp},'Error')
                            lapFR = squeeze(mean(lapFR_sample(:,:,corrInd==0),3));
                        else
                            lapFR = squeeze(mean(lapFR_sample,3));
                        end
                    elseif contains(PairID{ff,pp},'Test')
                        if contains(PairID{ff,pp},'Correct')
                            lapFR = squeeze(mean(lapFR_test(:,:,corrInd==1),3));
                        elseif contains(PairID{ff,pp},'Error')
                            lapFR = squeeze(mean(lapFR_test(:,:,corrInd==0),3));
                        else
                            lapFR = squeeze(mean(lapFR_test,3));
                        end
                    end
                    LapFR = cat(3,LapFR,lapFR);
                end
                frmat_norm = LapFR./repmat(max(LapFR,[],2),1,size(LapFR,2),1);
                
                FRmat = cat(1,FRmat,LapFR);
                FRmat_norm = cat(1,FRmat_norm,frmat_norm);
                FRmat_pre = cat(1,FRmat_pre,FR_PreRun);
                CellType = cat(1,CellType,celltype);
                SessionID = cat(1,SessionID,repmat({[RatID{rr},DateID{ii}]},size(celltype)));
                if strcmp(RatID{rr},'Rat139') % Create empty cells to match number of trials to other rats
                    lapSpkind_preRun = cat(2,lapSpkind_preRun,cell(size(lapSpkind_preRun,1),2));
                end
                LapSpkind_preRun = cat(1,LapSpkind_preRun,lapSpkind_preRun);
                LapSpkind_sample = cat(1,LapSpkind_sample,lapSpkind_sample);
                LapSpkind_test = cat(1,LapSpkind_test,lapSpkind_test);
                LapSpkind_correctInd = cat(1,LapSpkind_correctInd,repmat(corrInd',size(lapSpkind_sample,1),1));
            end
        end
    end
    FRmat_pre_norm = FRmat_pre./repmat(max(FRmat_pre,[],2),1,size(FRmat_pre,2));
    
    [~,maxInd] = max(squeeze(FRmat(:,:,1)),[],2);
    %[~,maxInd] = max(FRmat_pre,[],2);
    [PeakFRInd,sortInd] = sort(maxInd);
    FRmat = FRmat(sortInd,:,:);
    FRmat_norm = FRmat_norm(sortInd,:,:);
    FRmat_pre = FRmat_pre(sortInd,:);
    FRmat_pre_norm = FRmat_pre_norm(sortInd,:);
    CellType = CellType(sortInd);
    LapSpkind_preRun = LapSpkind_preRun(sortInd,:);
    LapSpkind_sample = LapSpkind_sample(sortInd,:);
    LapSpkind_test = LapSpkind_test(sortInd,:);
    LapSpkind_correctInd = LapSpkind_correctInd(sortInd,:);

    in = ~isnan(sum(sum(FRmat_norm,2),3));
    %in = max(max(FRmat,[],2),[],3)>1 & max(FRmat_pre,[],2);
    FRmat = FRmat(in,:,:);
    FRmat_norm = FRmat_norm(in,:,:);
    PeakFRInd = PeakFRInd(in);
    FRmat_pre = FRmat_pre(in,:);
    FRmat_pre_norm = FRmat_pre_norm(in,:);
    CellType = CellType(in);
    SessionID = SessionID(in);
    LapSpkind_preRun = LapSpkind_preRun(in,:);
    LapSpkind_sample = LapSpkind_sample(in,:);
    LapSpkind_test = LapSpkind_test(in,:);
    LapSpkind_correctInd = LapSpkind_correctInd(in,:);
    
    in = true(size(CellType));
    frmat_pre_norm = FRmat_pre_norm(in,:);
    relativePosAxis = posbinaxis-posbinaxis(ind_half);
    PeakFRLoc = relativePosAxis(PeakFRInd(in));    
    %% Make raster plot and mean firing rate plot for each cell
    if strcmp(PairID{1,pp},'Sample') && strcmp(PairID{2,pp},'TestCorrect') && strcmp(PairID{3,pp},'TestError')
        nSubplot = 5;
        nFigures = 0;
        n = 0;
        for ii = 1:size(LapSpkind_sample,1)
            n = n+1;
            ha = zeros(4,1);
            if n == 1
                h = figure;
                set(h,'OuterPosition',[343.4,107.4,968,728])
                nFigures = nFigures+1;
            end
            
            SpkposInd=[];
            TrialID=[];
            for trial = 1:size(LapSpkind_preRun,2)
                nSpikes = length(LapSpkind_preRun{ii,trial});
                SpkposInd = cat(1,SpkposInd,LapSpkind_preRun{ii,trial});
                TrialID = cat(1,TrialID,repmat(trial,nSpikes,1));
            end
            ha(1) = subplot(4,nSubplot,n);
            yyaxis left
            plot(SpkposInd/pi*180,TrialID,'.')
            ylim([.5 6.5])
            yyaxis right
            plot(posbinaxis-180,FRmat_pre(ii,:))
            title(['Cell ',num2str(ii)])
            if n == 1
               yyaxis left
               ylabel('PreRun trials')
            end
            yyaxis right            
            
            SpkposInd=[];
            TrialID=[];
            for trial = 1:size(LapSpkind_sample,2)
                nSpikes = length(LapSpkind_sample{ii,trial});
                SpkposInd = cat(1,SpkposInd,LapSpkind_sample{ii,trial});
                TrialID = cat(1,TrialID,repmat(trial,nSpikes,1));
            end
            ha(2) = subplot(4,nSubplot,n+nSubplot);
            yyaxis left
            plot(SpkposInd/pi*180,TrialID,'.')
            ylim([.5 8.5])
            yyaxis right
            plot(posbinaxis-180,FRmat(ii,:,1))
            if n == 1
               yyaxis left
               ylabel('Sample trials')
            end
            yyaxis right
            
            SpkposInd_crt=[];
            TrialID_crt=[];
            SpkposInd_err=[];
            TrialID_err=[];
            nCorr = 0;
            nError = 0;
            for trial = 1:size(LapSpkind_test,2)
                if LapSpkind_correctInd(ii,trial)==1
                    nCorr = nCorr+1;
                    nSpikes = length(LapSpkind_test{ii,trial});
                    SpkposInd_crt = cat(1,SpkposInd_crt,LapSpkind_test{ii,trial});
                    TrialID_crt = cat(1,TrialID_crt,repmat(nCorr,nSpikes,1));
                elseif LapSpkind_correctInd(ii,trial)==0
                    nError = nError+1;
                    nSpikes = length(LapSpkind_test{ii,trial});
                    SpkposInd_err = cat(1,SpkposInd_err,LapSpkind_test{ii,trial});
                    TrialID_err = cat(1,TrialID_err,repmat(nError,nSpikes,1));
                end
            end            
            ha(3) = subplot(4,nSubplot,n+2*nSubplot);
            yyaxis left
            plot(SpkposInd_crt/pi*180,TrialID_crt,'.')
            ylim([.5 nCorr+.5])
            yyaxis right
            plot(posbinaxis-180,FRmat(ii,:,2))
            if n == 1
               yyaxis left
               ylabel('Test correct trials')
            end
            yyaxis right
            
            ha(4) = subplot(4,nSubplot,n+3*nSubplot);
            yyaxis left
            plot(SpkposInd_err/pi*180,TrialID_err,'.')
            ylim([.5 nError+.5])
            yyaxis right
            plot(posbinaxis-180,FRmat(ii,:,3))
            if n == 1
               yyaxis left
               ylabel('Test error trials')
            end
            yyaxis right
            yvalues = cell2mat(get(ha,'YLim'));
            set(ha,'YLim',[max(min(yvalues(:)),0) max(yvalues(:))])
            
            if n == nSubplot || ii == size(LapSpkind_sample,1)
               saveas(h,strcat(outdir_sub,'\PositionTuningEachCell',num2str(nFigures)),'epsc')
               saveas(h,strcat(outdir_sub,'\PositionTuningEachCell',num2str(nFigures)),'png')
               close(h)
               n = 0;
            end
        end
    end
    
    %% Plot raw firing rate during correct and error trials
    h = figure;
    set(h,'Outerposition',[272.2,327.4,923.2,508])
    subplot(1,size(PairID,1)+1,1)
        imagesc(posbinaxis-posbinaxis(ind_half),1:size(FRmat_pre,1),FRmat_pre)
        axis xy square
        colorbar
        ylabel('CA1 cell index')
        title('Pre-running trials')
    for ff = 1:size(PairID,1)
        subplot(1,size(PairID,1)+1,ff+1)
            imagesc(posbinaxis-posbinaxis(ind_half),1:size(FRmat_pre,1),squeeze(FRmat(:,:,ff)))
            axis xy square
            colorbar
            title(PairID{ff,pp})
            xlabel('Position relative to reward (deg)')
    end
    colormap jet
    saveas(h,strcat(outdir_sub,'\PositionTuning'),'epsc')
    saveas(h,strcat(outdir_sub,'\PositionTuning'),'png')
    close (h)
    %% Plot normalized firing rate during correct and error trials
    h = figure;
    set(h,'Outerposition',[272.2,327.4,923.2,508])
    subplot(1,size(PairID,1)+1,1)
        imagesc(posbinaxis-posbinaxis(ind_half),1:sum(in),frmat_pre_norm)
        axis xy square
        colorbar
        set(gca,'CLim',[0 1])
        ylabel('CA1 cell index')
        title('Pre-running trials')
    for ff = 1:size(PairID,1)
        subplot(1,size(PairID,1)+1,ff+1)
            imagesc(posbinaxis-posbinaxis(ind_half),1:sum(in),squeeze(FRmat_norm(in,:,ff)))
            axis xy square
            colorbar
            set(gca,'CLim',[0 1])
            title(PairID{ff,pp})
            xlabel('Position relative to reward (deg)')
    end
    colormap jet
    saveas(h,strcat(outdir_sub,'\PositionTuning_norm'),'epsc')
    saveas(h,strcat(outdir_sub,'\PositionTuning_norm'),'png')
    close (h)
    Output.posbinaxis = posbinaxis;

    %% Plot time lag Cross correlation between correct and error position tunnings 
    % align position tunning to the middle to avoid edege effect from xcorr
    [~,maxInd] = max(FRmat(:,:,1),[],2);
    shift = maxInd-ceil(size(FRmat,2)/2);
    FRmat_s = FRmat;
    for condition = 1:size(FRmat,3)
        for cellid = 1:size(FRmat,1)
            FRmat_s(cellid,:,condition) = circshift(FRmat_s(cellid,:,condition),-shift(cellid),2);
        end
    end
    % perform xcorr
    lagAxis = -100:4:100;
    XCorrMat = zeros(size(FRmat_s,1),length(lagAxis),size(FRmat,3)-1);
    for ff = 1:size(XCorrMat,3)
        for ii = 1:size(FRmat_s,1)
            XCorrMat(ii,:,ff) = xcorr(FRmat_s(ii,:,1),FRmat_s(ii,:,ff+1),floor(length(lagAxis)/2),'coeff');
        end
    end
    Output.XCorrMat = XCorrMat;
    h_xcorr = figure;
    set(gcf,'OuterPosition',[481,95.4,812.8,747.2])
    KernelMap = [];
    for ff = 1:size(XCorrMat,3)
        subplot(size(XCorrMat,3),2,ff)
        imagesc(lagAxis,1:size(XCorrMat,1),XCorrMat(:,:,ff));
        xlim([-50 50])
        ylabel('Cell index')
        title(PairID{ff+1,pp})
        colorbar
        axis xy square

        % Estimate 2D map of Cross correlation coefficient as a function of
        % position lag and peak firing position using kernel density estimation
        [kernelMap] = PFR(squeeze(XCorrMat(:,:,ff))',PeakFRLoc/180*pi,lagAxis,relativePosAxis/180*pi);
        ha = subplot(size(XCorrMat,3),2,ff+2);
        imagesc(lagAxis,relativePosAxis,kernelMap'); hold on
        line(ha,[0 0],ylim,'Color','white','LineStyle','--','LineWidth',1)
        xlim([-50 50])
        str = sprintf('Position shift relative to %s (deg)',PairID{1,pp});
        xlabel(str)
        ylabel('Position aligned to reward (deg)')
        cb = colorbar;
        ylabel(cb,'Cross correlation coefficient')
        axis xy square
        
        KernelMap = cat(3,KernelMap,kernelMap);
    end
    Output.KernelMap = KernelMap;
    Output.lagAxis = lagAxis;
    Output.relativePosAxis = relativePosAxis;

    saveas(h_xcorr,strcat(outdir_sub,'\XCorr_allcell'),'epsc')
    saveas(h_xcorr,strcat(outdir_sub,'\XCorr_allcell'),'png')
    close (h_xcorr)
    save(strcat(outdir_sub,'\Results.mat'),'Output')
    
    % Plot correlation between pairs of conditions
    c_value = zeros(size(FRmat_pre,1),6);
    for ii = 1:size(FRmat_pre,1)
        c_value(ii,1) = corr(FRmat_pre(ii,:)',squeeze(FRmat(ii,:,1))');
        c_value(ii,2) = corr(FRmat_pre(ii,:)',squeeze(FRmat(ii,:,2))');
        c_value(ii,3) = corr(FRmat_pre(ii,:)',squeeze(FRmat(ii,:,3))');
        c_value(ii,4) = corr(squeeze(FRmat(ii,:,1))',FRmat(ii,:,2)');
        c_value(ii,5) = corr(squeeze(FRmat(ii,:,1))',FRmat(ii,:,3)');
        c_value(ii,6) = corr(squeeze(FRmat(ii,:,2))',FRmat(ii,:,3)');
    end
    Pairs_Label = {['PreRun-',PairID{1}];...
                   ['PreRun-',PairID{2}];...
                   ['PreRun-',PairID{3}];...
                   [PairID{1},'-',PairID{2}];...
                   [PairID{1},'-',PairID{3}];...
                   [PairID{2},'-',PairID{3}]};
    CI = bootci(5000,@nanmean,c_value);
    CI(1,:) = nanmean(c_value)-CI(1,:);
    CI(2,:) = CI(2,:)-nanmean(c_value);
    h = figure;
    errorbar(1:length(Pairs_Label),nanmean(c_value),CI(1,:),CI(2,:),'k-','LineWidth',1)
    axis square
    xlim([.5 length(Pairs_Label)+.5])
    ylimt = get(gca,'YLim');
    ylimt(1) = 0;
    set(gca,'YLim',ylimt,'XTick',1:length(Pairs_Label),'XTickLabel',Pairs_Label,'Box','off')
    ylabel('Pearson r')
    saveas(h,strcat(outdir_sub,'\MeanCorrAcrossCondition'),'epsc')
    saveas(h,strcat(outdir_sub,'\MeanCorrAcrossCondition'),'png')
    close (h)
    
    [~,tbl,stats] = kruskalwallis(c_value,[],'off');
    posthoc = multcompare(stats,'Display','off');
    Stats.KKW = tbl;
    Stats.posthoc = posthoc;
    save(strcat(outdir_sub,'\Stats.mat'),'Stats')
end
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

function [pfr] = PFR(TFR,phase,freqVec,pbins)

kappa = 90;
delta = angle(exp(1i*repmat(phase',1,length(pbins))).*conj(exp(1i*repmat(pbins,length(phase),1))));
W = exp(kappa*cos(delta));
D = ones(length(freqVec),length(phase))*exp(kappa*cos(delta));
pfr = TFR*W./D;

function [delta] = circDist(x,y)
delta = angle(exp(1i*x).*conj(exp(1i*y)));

function [CI_u, CI_l] = CorrectionForMultipleComparsion(sample)
%Simultaneous bounds are wider than separate bounds, because it is more stringent to require that the entire 
%curve be within the bounds than to require that the curve at a single predictor value be within the bounds.
[sample_norm] = tiedrank(sample);

sample_min = prctile(min(min(sample_norm,[],2),[],3),2.5);
sample_max = prctile(max(max(sample_norm,[],2),[],3),97.5);

CI_u=zeros(size(sample,2),size(sample,3));
for ii = 1:size(sample,2)
    for i2 = 1:size(sample,3)
        if sum(sample_norm(:,ii,i2)>=sample_max) > 0
            CI_u(ii,i2) = min(sample(sample_norm(:,ii,i2)>=sample_max,ii,i2));
        else
            CI_u(ii,i2) = max(sample(:,ii,i2));
        end
    end
end

CI_l=zeros(size(sample,2),size(sample,3));
for ii = 1:size(sample,2)
    for i2 = 1:size(sample,3)
        if sum(sample_norm(:,ii,i2)<=sample_min) > 0
            CI_l(ii,i2) = max(sample(sample_norm(:,ii,i2)<=sample_min,ii,i2));
        else
            CI_l(ii,i2) = min(sample(:,ii,i2));
        end
    end
end