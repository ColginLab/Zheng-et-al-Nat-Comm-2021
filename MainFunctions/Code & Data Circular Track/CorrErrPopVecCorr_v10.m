%% Plot correlation matrix between correct and error trials (Separate reward preexiting and reward emerged cells)
% This version makes pairs of trials into variables
% This version aligns error trials to stop location
function CorrErrPopVecCorr_v10(trackdata,AnalysisFD)

if exist(strcat(AnalysisFD,'\LapFRcorr\LapFRcorr.mat'),'file') ~= 2
    error('Please run RunPlaceCellPropertiesAcrossLap_v4.m first')
else
    load(strcat(AnalysisFD,'\LapFRcorr\LapFRcorr.mat'))
end

outdir = strcat(AnalysisFD,'\CorrErrPopVecCorr_cellType_alignedToStopLoc');
if ~isdir(outdir)
    mkdir(outdir)
end

nLapThr = 3;
RatID = fieldnames(LapFRcorr);
PairID = {'SampleCorrect','SampleError','Sample';...
          'TestCorrect'  ,'TestCorrect','TestCorrect';...
          'TestError'    ,'TestError'  ,'TestError'}; 
for pp = 1:size(PairID,2)
    outdir_sub = strcat(outdir,'\',PairID{1,pp},'_',PairID{2,pp},'_',PairID{3,pp});
    if ~isdir(outdir_sub)
        mkdir(outdir_sub)
    end      
    FRmat = [];
    FRmat_norm = [];
    FRmat_pre = [];
    CellType = {};
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
            corrInd = trackInfo.sign_correct_test;
            stopLoc = trackInfo.ang_test_reward_ontrack;
            
            % Align stop position to middle        
            ind = match(stopLoc,posbinaxis/180*pi);
            ind_half = match(180,posbinaxis);
            shift = ind-ind_half;
            shift(isnan(stopLoc)) = 0; % don't shift the trials where rats didn't stop
            for ss = 1:length(shift)                           
                lapFR_test(:,:,ss) = circshift(lapFR_test(:,:,ss),-shift(ss),2);
            end
            % Align reward position to middle
            ind = match(rewardLoc,posbinaxis/180*pi);
            ind_half = match(180,posbinaxis);
            shift = ind-ind_half;
            FR_PreRun = circshift(FR_PreRun,-shift,2);
            lapFR_sample = circshift(lapFR_sample,-shift,2);

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
            end
        end
    end
    FRmat_pre_norm = FRmat_pre./repmat(max(FRmat_pre,[],2),1,size(FRmat_pre,2));
    
    [~,maxInd] = max(squeeze(FRmat(:,:,1)),[],2);
    [PeakFRInd,sortInd] = sort(maxInd);
    FRmat = FRmat(sortInd,:,:);
    FRmat_norm = FRmat_norm(sortInd,:,:);
    FRmat_pre = FRmat_pre(sortInd,:);
    FRmat_pre_norm = FRmat_pre_norm(sortInd,:);
    CellType = CellType(sortInd);
    

    in = ~isnan(sum(sum(FRmat_norm,2),3)) & (strcmp(CellType,'RewardPre')|strcmp(CellType,'RewardEmerged'));
    FRmat = FRmat(in,:,:);
    FRmat_norm = FRmat_norm(in,:,:);
    PeakFRInd = PeakFRInd(in);
    FRmat_pre = FRmat_pre(in,:);
    FRmat_pre_norm = FRmat_pre_norm(in,:);
    CellType = CellType(in);

    celltype = unique(CellType);
    for tt = 1:length(celltype)
        in = strcmp(CellType,celltype{tt});
        frmat_pre_norm = FRmat_pre_norm(in,:);
        frmat = FRmat(in,:,:);
        frmat_norm = FRmat_norm(in,:,:);
        relativePosAxis = posbinaxis-posbinaxis(ind_half);
        PeakFRLoc = relativePosAxis(PeakFRInd(in));
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
                imagesc(posbinaxis-posbinaxis(ind_half),1:sum(in),squeeze(frmat_norm(:,:,ff)))
                axis xy square
                colorbar
                set(gca,'CLim',[0 1])
                title(PairID{ff,pp})
                xlabel('Position relative to stop (deg)')
        end
        colormap hot
        saveas(h,strcat(outdir_sub,'\PositionTuning_',celltype{tt}),'epsc')
        saveas(h,strcat(outdir_sub,'\PositionTuning_',celltype{tt}),'png')
        close (h)
        Output.(celltype{tt}).posbinaxis = posbinaxis;

        %% Plot time lag Cross correlation between correct and error position tunnings 
        % align position tunning to the middle to avoid edege effect from xcorr
        [~,maxInd] = max(frmat(:,:,1),[],2);
        shift = maxInd-ceil(size(frmat,2)/2);
        FRmat_s = frmat;
        for condition = 1:size(frmat,3)
            for cellid = 1:size(frmat,1)
                FRmat_s(cellid,:,condition) = circshift(FRmat_s(cellid,:,condition),-shift(cellid),2);
            end
        end
        % perform xcorr
        lagAxis = -100:4:100;
        XCorrMat = zeros(size(FRmat_s,1),length(lagAxis),size(frmat,3)-1);
        for ff = 1:size(XCorrMat,3)
            for ii = 1:size(FRmat_s,1)
                XCorrMat(ii,:,ff) = xcorr(FRmat_s(ii,:,1),FRmat_s(ii,:,ff+1),floor(length(lagAxis)/2),'coeff');
            end
        end
        Output.(celltype{tt}).XCorrMat = XCorrMat;
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
            ylabel('Position aligned to stop (deg)')
            cb = colorbar;
            ylabel(cb,'Cross correlation coefficient')
            axis xy square

            KernelMap = cat(3,KernelMap,kernelMap);
        end
        Output.(celltype{tt}).KernelMap = KernelMap;
        Output.(celltype{tt}).lagAxis = lagAxis;
        Output.(celltype{tt}).relativePosAxis = relativePosAxis;

        saveas(h_xcorr,strcat(outdir_sub,'\XCorr_',celltype{tt}),'epsc')
        saveas(h_xcorr,strcat(outdir_sub,'\XCorr_',celltype{tt}),'png')
        close (h_xcorr)
        save(strcat(outdir_sub,'\Results.mat'),'Output')
    end
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

function [CI_u, CI_l] = CorrectionForMultipleComparsion(sample)
%Simultaneous bounds are wider than separate bounds, because it is more stringent to require that the entire 
%curve be within the bounds than to require that the curve at a single predictor value be within the bounds.
[sample_norm] = zscore(sample);

sample_min = prctile(min(min(sample_norm,[],2),[],3),5);
sample_max = prctile(max(max(sample_norm,[],2),[],3),95);

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