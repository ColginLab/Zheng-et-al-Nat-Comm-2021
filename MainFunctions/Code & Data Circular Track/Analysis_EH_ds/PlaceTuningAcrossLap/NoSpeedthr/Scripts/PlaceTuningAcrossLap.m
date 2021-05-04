function PlaceTuningAcrossLap(DataFD,AnalysisFD)
Required_files = {strcat(AnalysisFD,'\LapFRcorr\LapFRcorr.mat');...
                  strcat(AnalysisFD,'\LapFRcorr_noSpeedThr\LapFRcorr.mat')};
Required_scripts = {'RunPlaceCellPropertiesAcrossLap_v4.m';...
                    'RunPlaceCellPropertiesAcrossLap_v5.m'};
outdir_name = {strcat(AnalysisFD,'\PlaceTuningAcrossLap\Speedthr');...
              strcat(AnalysisFD,'\PlaceTuningAcrossLap\NoSpeedthr')};
          
for file = 1:length(Required_files)
    if exist(Required_files{file},'file') ~= 2
        error(['Please run ',Required_scripts{file},' first'])
    else
        load(Required_files{file})
    end

    outdir = outdir_name{file};
    if ~isdir(outdir)
        mkdir(outdir)
    end

    ratID = fieldnames(LapFRcorr);
    for rr = 1:length(ratID)
        dateID = fieldnames(LapFRcorr.(ratID{rr}));
        DateFd = translateDateID(dateID);
        lapFR_rat = [];
        for dd = 1:length(DateFd)
            % Load reward position
            if DateFd{dd}(end) == 'T'
                trackfile = strcat(DataFD,'\',ratID{rr},'\',DateFd{dd},'\',dateID{dd}(2:9),'_CT_tracking.mat');
            else
                trackfile = strcat(DataFD,'\',ratID{rr},'\',DateFd{dd},'\',dateID{dd}(2:9),'_CT_tracking_',dateID{dd}(end),'.mat');
            end
            trackdata = load(trackfile);
            trackdata.sign_correct_sample(isnan(trackdata.sign_correct_sample)) = false;
            rewardPos = mode(trackdata.ang_sample_reward_ontrack(logical(trackdata.sign_correct_sample)));
            rewardPos = rewardPos/2/pi*360;

            trialID = fieldnames(LapFRcorr.(ratID{rr}).(dateID{dd}));
            lapFR = [];
            active = false(size(LapFRcorr.(ratID{rr}).(dateID{dd}).(trialID{1}).lapFR,1),1);
            for tt = 1:length(trialID)
                active(LapFRcorr.(ratID{rr}).(dateID{dd}).(trialID{tt}).active) = true;
                if tt == 2
                    holder = NaN(size(LapFRcorr.(ratID{rr}).(dateID{dd}).(trialID{tt}).lapFR));
                    holder = repmat(holder,1,1,2);
                    holder(:,:,1:2:size(holder,3)-1) = LapFRcorr.(ratID{rr}).(dateID{dd}).(trialID{tt}).lapFR;
                    holder(:,:,2:2:size(holder,3)) = LapFRcorr.(ratID{rr}).(dateID{dd}).(trialID{tt+1}).lapFR;
                    lapFR = cat(3,lapFR,holder);
                elseif tt == 3
                    continue
                else
                    lapFR = cat(3,lapFR,LapFRcorr.(ratID{rr}).(dateID{dd}).(trialID{tt}).lapFR);
                end            
            end
            lapFR = lapFR(active,:,:);
            n = ceil(sqrt(size(lapFR,1)));
            nPreRunLap = size(LapFRcorr.(ratID{rr}).(dateID{dd}).(trialID{1}).lapFR,3);
            nSampleLap = size(LapFRcorr.(ratID{rr}).(dateID{dd}).(trialID{2}).lapFR,3);
            nPostTestLap = size(LapFRcorr.(ratID{rr}).(dateID{dd}).(trialID{4}).lapFR,3);
            temp = zeros(nPreRunLap+nSampleLap*2+nPostTestLap,1);
            temp(nPreRunLap+1:2:nPreRunLap+nSampleLap*2-1) = trackdata.sign_correct_sample;
            temp(nPreRunLap+2:2:nPreRunLap+nSampleLap*2) = trackdata.sign_correct_test;
            temp(nPreRunLap+nSampleLap*2+1:end) = trackdata.Sign_correct_posttest{1};
            ReceiveReward = false(nPreRunLap+nSampleLap*2+nPostTestLap,1);
            ReceiveReward(temp == 1) = true;

            posbinaxis = LapFRcorr.(ratID{rr}).(dateID{dd}).(trialID{tt}).posbinaxis;
            h = figure;
            set(h,'OuterPosition',[505,57,1240,1016])
            for ii = 1:size(lapFR,1)
                lapFR_cell = squeeze(lapFR(ii,:,:));
                subplot(n,n,ii)
                imagesc(1:size(lapFR_cell,2),posbinaxis,lapFR_cell)
                axis xy square tight
                hold on
                plot([nPreRunLap+1 size(lapFR_cell,2)+.5],[rewardPos rewardPos],'w-')
                if ii>1
                    set(gca,'TickDir','out','XTick',find(ReceiveReward),'XTickLabel',{},'YTickLabel',{})
                end
                title([num2str(round(max(max(lapFR_cell))*10)/10),' Hz'])
            end
            colormap jet
            saveas(h,strcat(outdir,'\PlaceTuning_',ratID{rr},'_',dateID{dd}),'epsc')
            saveas(h,strcat(outdir,'\PlaceTuning_',ratID{rr},'_',dateID{dd}),'png')
            close(h)
            % Align to reward position
            ind = match(rewardPos,posbinaxis);
            shift = round(ind-length(posbinaxis)/2);
            lapFR = circshift(lapFR,-shift,2);
            lapFR_rat = cat(1,lapFR_rat,lapFR);
        end

        [h,h2] = CellProportionPlot(lapFR_rat,posbinaxis,nPreRunLap,nSampleLap,ratID{rr});      
        saveas(h,strcat(outdir,'\CellPropAcrossLap_',ratID{rr}),'epsc')
        saveas(h,strcat(outdir,'\CellPropAcrossLap_',ratID{rr}),'fig')
        saveas(h2,strcat(outdir,'\CellPropAcrossLap_',ratID{rr},'_rawFR'),'epsc')
        saveas(h,strcat(outdir,'\CellPropAcrossLap_',ratID{rr}),'png')
        saveas(h2,strcat(outdir,'\CellPropAcrossLap_',ratID{rr},'_rawFR'),'png')
        close(h,h2)

        [h] = ExampleCells(lapFR_rat,posbinaxis,nPreRunLap,nSampleLap,ratID{rr});
        saveas(h,strcat(outdir,'\ExampleCells_',ratID{rr}),'epsc')
        saveas(h,strcat(outdir,'\ExampleCells_',ratID{rr}),'png')
        close(h)

    %     [~,peakFR_loc] = max(mean(lapFR_rat(:,:,1:nPreRunLap),3),[],2);
    %     [~,sortInd] = sort(peakFR_loc);
    %     lapFR_rat_sorted = lapFR_rat(sortInd,:,:);
    %     baseline = mean(lapFR_rat_sorted(:,:,1:nPreRunLap),3);
    %     n = ceil(sqrt(size(lapFR_rat_sorted,3)));
    %     h = figure;
    %     set(h,'OuterPosition',[505,57,1240,1016])
    %     for ll = 1:size(lapFR_rat_sorted,3)
    %         lapFR_pop = squeeze(lapFR_rat_sorted(:,:,ll));
    %         lapFR_pop_norm = squeeze((lapFR_pop-baseline)./(lapFR_pop+baseline));
    %         lapFR_pop_norm(isnan(lapFR_pop_norm)) = 0;
    %         subplot(n,n,ll)
    %         imagesc(posbinaxis-180,1:size(lapFR_rat_sorted,2),lapFR_pop_norm)
    %         axis xy square tight
    %         colormap(redblue)
    %         set(gca,'CLim',[-max(max(abs(lapFR_pop_norm))) max(max(abs(lapFR_pop_norm)))])
    %         if ll>1
    %             set(gca,'XTickLabel',{},'YTickLabel',{})
    %             title(['Lap',num2str(ll)])
    %         else
    %             title(ratID{rr})
    %             ylabel('Cell ID')
    %             xlabel('Position aligned to reward')
    %         end
    %     end
    %     saveas(h,strcat(outdir,'\PopFRAcrossLap_',ratID{rr}),'fig')
    %     saveas(h,strcat(outdir,'\PopFRAcrossLap_',ratID{rr}),'png')
    %     close(h)

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
    % save(strcat(outdir,'\Stats.mat'),'Stats')
end
function [h,h2] = CellProportionPlot(lapFR_rat,posbinaxis,nPreRunLap,nSampleLap,ratID)
    FRthr = 0.5; % in Hz
    nCell = sum(lapFR_rat>FRthr);
    CellProportion = squeeze(nCell)'/size(lapFR_rat,1);
    %nulldist = CellProportion_shu(lapFR_rat,FRthr);
    %CI = prctile(nulldist,97.5,3);
    
    h = figure;
    set(h,'OuterPosition',[1030,-35,494,1102])   
    h1 = subplot(3,1,1);
    imagesc(1:size(lapFR_rat,3),posbinaxis-180,CellProportion');
    %colormap(gca,'plasma');
    axis xy square
    set(gca,'CLim',[0 0.3],'YLim',[-100 100])
    hold on;
    plot([nPreRunLap+.5 nPreRunLap+.5],ylim,'w--','LineWidth',1)
    plot([nPreRunLap+nSampleLap*2+.5 nPreRunLap+nSampleLap*2+.5],ylim,'w--','LineWidth',1)
    ylabel('Position aligned to reward')
    c = colorbar;
    ylabel(c,'Cell proportion')
    title([ratID,' FR threshold:',num2str(FRthr),'Hz'])
    
    FRthr = 10; % in Hz
    nCell = sum(lapFR_rat>FRthr);
    CellProportion = squeeze(nCell)'/size(lapFR_rat,1);
    %nulldist = CellProportion_shu(lapFR_rat,FRthr);
    %CI = prctile(nulldist,2.5,3);
    
    h2 = subplot(3,1,2); 
    imagesc(1:size(lapFR_rat,3),posbinaxis-180,CellProportion');
    %colormap(gca,'plasma');
    axis xy square
    set(gca,'CLim',[0 0.1])
    hold on
    plot([nPreRunLap+.5 nPreRunLap+.5],ylim,'w--','LineWidth',1)
    plot([nPreRunLap+nSampleLap*2+.5 nPreRunLap+nSampleLap*2+.5],ylim,'w--','LineWidth',1)
    c = colorbar;
    ylabel(c,'Cell proportion')
    title(['FR threshold:',num2str(FRthr),'Hz'])
    
    preRun_meanFR = squeeze(mean(lapFR_rat(:,:,1:nPreRunLap),3));
    [~,sortInd] = sort(preRun_meanFR,'descend');
    nCell = round(size(lapFR_rat,1)*.1); % choose top 10% of cells
    NormalizeFRAcrossLap = zeros(size(lapFR_rat,2),size(lapFR_rat,3));
    RawFRAcrossLap = zeros(size(lapFR_rat,2),size(lapFR_rat,3));
    for ii = 1:size(NormalizeFRAcrossLap,1)
        lapFR = squeeze(lapFR_rat(sortInd(1:nCell,ii),ii,:));
        baseline = preRun_meanFR(sortInd(1:nCell,ii),ii);
        NormalizeFRAcrossLap(ii,:) = mean(lapFR./repmat(baseline,1,size(lapFR,2)));
        RawFRAcrossLap(ii,:) = mean(lapFR);
    end
    h3 = subplot(3,1,3);
    imagesc(1:size(lapFR_rat,3),posbinaxis-180,NormalizeFRAcrossLap*100);
    colormap(gca,'redblue');                   % redblue?
    axis xy square
    set(gca,'CLim',[0 200])
    hold on
    plot([nPreRunLap+.5 nPreRunLap+.5],ylim,'k--','LineWidth',1)
    plot([nPreRunLap+nSampleLap*2+.5 nPreRunLap+nSampleLap*2+.5],ylim,'k--','LineWidth',1)
    c = colorbar;
    ylabel(c,'Mean % firing rate from baseline')
    xlabel('Lap #')
    title(['Top ',num2str(nCell),' cells'])
    linkaxes([h1 h2 h3],'x')
    
    h2 = figure;
    imagesc(1:size(lapFR_rat,3),posbinaxis-180,RawFRAcrossLap);
    colormap(gca,'jet');
    axis xy square
    ylabel('Position aligned to reward')
    c = colorbar;
    ylabel(c,'Mean firing rate (Hz)')
    xlabel('Lap #')
    title(['Top ',num2str(nCell),' cells'])

function [h] = ExampleCells(lapFR_rat,posbinaxis,nPreRunLap,nSampleLap,ratID)    
    nCell = 5;
    FRthr = 1; % in Hz
    ind = match(0,posbinaxis-180);
    [~,sortInd_pre] = sort(squeeze(mean(mean(lapFR_rat(:,ind-1:ind+1,1:nPreRunLap),2),3)),'descend');
    Reward_pre = squeeze(mean(mean(lapFR_rat(sortInd_pre,ind-1:ind+1,1:nPreRunLap),2),3))>=FRthr;
    Reward_pre_notSorted = squeeze(mean(mean(lapFR_rat(:,ind-1:ind+1,1:nPreRunLap),2),3))>=FRthr;
    [~,sortInd_emerge] = sort(squeeze(mean(mean(lapFR_rat(:,ind-1:ind+1,nPreRunLap+1:end),2),3)),'descend');
    Reward_emerge = squeeze(mean(mean(lapFR_rat(sortInd_emerge,ind-1:ind+1,nPreRunLap+1:end),2),3))>=FRthr & ~Reward_pre_notSorted(sortInd_emerge);
    h = figure;
    set(h,'OuterPosition',[481,39,678,1004])
    for ii = 1:nCell
        subplot(nCell,2,2*ii-1)
        ind = find(Reward_pre,ii);
        temp = lapFR_rat(sortInd_pre,:,:);
        lapFR_cell = squeeze(temp(ind(end),:,:));
        imagesc(1:size(lapFR_cell,2),posbinaxis-180,lapFR_cell)
        axis xy square tight
        hold on
        plot([nPreRunLap+.5 nPreRunLap+.5],ylim,'w--')
        plot([nPreRunLap+nSampleLap*2+.5 nPreRunLap+nSampleLap*2+.5],ylim,'w--')
        if ii == 1
            titlestr = sprintf([ratID,'\nReward_p_r_e max FR:',num2str(max(max(round(lapFR_cell*10)/10))),'Hz']);
            title(titlestr);
            ylabel('Position aligned to reward')
            xlabel('Lap #')
        else
            title([num2str(max(max(round(lapFR_cell*10)/10))),'Hz'])
        end
        
        subplot(nCell,2,2*ii)
        ind = find(Reward_emerge,ii);
        temp = lapFR_rat(sortInd_emerge,:,:);
        lapFR_cell = squeeze(temp(ind(end),:,:));
        imagesc(1:size(lapFR_cell,2),posbinaxis-180,lapFR_cell)
        axis xy square tight
        hold on
        plot([nPreRunLap+.5 nPreRunLap+.5],ylim,'w--')
        plot([nPreRunLap+nSampleLap*2+.5 nPreRunLap+nSampleLap*2+.5],ylim,'w--')        
        if ii == 1
            title(['Reward_e_m_e_r_g_e max FR:',num2str(max(max(round(lapFR_cell*10)/10))),'Hz'])
        else
            title([num2str(max(max(round(lapFR_cell*10)/10))),'Hz'])
        end
    end
    colormap jet

function DateFd = translateDateID(dateID)
DateFd = cell(size(dateID));
for dd = 1:length(dateID)
    str = dateID{dd};
    str2 = strcat(str(2:5),'-',str(6:7),'-',str(8:9),'-',str(10:11));
    if length(str) > 11
        str2 = strcat(str2,'-',str(end));
    end
    DateFd{dd} = str2;
end

function nulldist = CellProportion_shu(lapFR,FRthr)
nShuffle = 5000;
nulldist = zeros(size(lapFR,3),size(lapFR,2),nShuffle);
for ss = 1:nShuffle
    shu = randperm(size(lapFR,2));   
    nCell = sum(lapFR(:,shu,:)>FRthr);
    nulldist(:,:,ss) = squeeze(nCell)'/size(lapFR,1);
end