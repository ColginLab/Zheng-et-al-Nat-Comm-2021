function RunPlaceCellPropertiesAcrossLap_v6(Runfd,trackdata,Outdir,cellList,spktime_fd)
% This version use histogram instead of kernel density estimate
Condition = {'Prerunning';'Sample';'Test';'Postrunning'};
outdir = strcat(Outdir,'\LapFRcorr');
if ~isdir(outdir)
    mkdir(outdir)
end

for rr = 1:length(Runfd)
    ind = strfind(Runfd{rr},'Rat');
    RatID = Runfd{rr}(ind:ind+5);    
    ind = strfind(Runfd{rr},'\');
    DateID = Runfd{rr}(ind(end-1)+1:end-1);
    ind = strfind(DateID,'-');
    DateID(ind)=[];
    DateID = strcat('D',DateID);
    spikeIndfile = strcat(Runfd{rr},'\',spktime_fd,'\spike_index.mat');
    eeg2posIndfile = strcat(Runfd{rr},'\',spktime_fd,'\eeg2pos.mat');
    posIndfile = strcat(Runfd{rr},'\Data_angle_ontrack.mat');
    spkfile = strcat(Runfd{rr},'\',cellList);
    trackInfo = load(trackdata{rr});
    lapFR_all = [];
    disp(Runfd{rr})
    for tt = 1:length(Condition)

        [lapFR,lapSpkind,cellID,posbinaxis,active,CorrectInd,corrMat] = ...
        PlaceCellPropertiesAcrossLap_v4(spkfile,spikeIndfile,eeg2posIndfile,posIndfile,Condition,trackInfo);
        lapFR_all = cat(3,lapFR_all,lapFR.(Condition{tt}));
        %% Construct output structure
        LapFRcorr.(RatID).(DateID).(Condition{tt}).lapFR = lapFR.(Condition{tt});
        LapFRcorr.(RatID).(DateID).(Condition{tt}).lapSpkind = lapSpkind.(Condition{tt});
        LapFRcorr.(RatID).(DateID).(Condition{tt}).cellID = cellID;
        LapFRcorr.(RatID).(DateID).(Condition{tt}).posbinaxis = posbinaxis;
        LapFRcorr.(RatID).(DateID).(Condition{tt}).active = active.(Condition{tt});
        LapFRcorr.(RatID).(DateID).(Condition{tt}).CorrectInd = CorrectInd.(Condition{tt});
        LapFRcorr.(RatID).(DateID).(Condition{tt}).rewardLoc = trackInfo.Ang_RewardLoc_ontrack(trackInfo.Ind_rewardloc)/pi*180;
        LapFRcorr.(RatID).(DateID).(Condition{tt}).corrMat = corrMat.(Condition{tt});          
    end
    CorrMat.(RatID).(DateID) = PopVecCorr(lapFR_all);
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
save(strcat(outdir,'\LapFRcorr.mat'),'LapFRcorr','CorrMat','FunctionName')

%% Plot correlation matrix
RatID = fieldnames(LapFRcorr);
for rr = 1:length(RatID)
    DateID = fieldnames(LapFRcorr.(RatID{rr}));
    nFigures = 0;
    nSubplot = 0; 
    for ii = 1:length(DateID)        
        nSubplot = nSubplot+1;
        condition = fieldnames(LapFRcorr.(RatID{rr}).(DateID{ii}));
        if nSubplot == 1
            nFigures = nFigures+1;
            h = figure;
            set(h,'OuterPosition',[110,269,1615,705])
            h2 = figure;
            set(h2,'OuterPosition',[93,400,1738,673])
            ha = zeros(length(condition),5);
            ha2 = zeros(1,5);
        end
        figure(h)
        for cc = 1:length(condition)
            corrMat= LapFRcorr.(RatID{rr}).(DateID{ii}).(condition{cc}).corrMat;

            ha(cc,nSubplot) = subplot(length(condition),5,(cc-1)*5+nSubplot);
            imagesc(corrMat); axis square xy
            %colormap hot;
            set(gca,'CLim',[0,1])
            if cc == 1
                title(DateID{ii},'Interpreter','none')
            elseif cc == 2
                xlabel('Lap number')
            end
            if nSubplot == 1
                label = sprintf([condition{cc},'\nLap number']);
                ylabel(label)
            end
            
            if  cc == length(condition)
                p = get(ha(cc,nSubplot),'position');
                c = colorbar;
                p2 = get(c,'position');
                set(ha(cc,nSubplot),'position',p)
                %p2(1) = .8; set(c,'position',p2);
                ylabel(c,'Population vector correlation')
            end
        end
        
        corrMat_all= CorrMat.(RatID{rr}).(DateID{ii});
        corrMat_all = corrMat_all-diag(NaN(length(corrMat_all),1)); % remove diagonal elements for visualization
        imAlpha=ones(size(corrMat_all));
        imAlpha(isnan(corrMat_all))=0;
        figure(h2)        
        ha2(1,nSubplot) = subplot(1,5,nSubplot);
        img = imagesc(corrMat_all);
        set(img,'AlphaData',imAlpha)
        axis square xy
        colorbar;
        %set(gca,'CLim',[0,1])
        title(DateID{ii},'Interpreter','none')
        xlabel('Lap number')
        if nSubplot == 1
            ylabel('Lap number')
        end
        
        if nSubplot == 5 || ii == length(DateID)
            saveas(h,strcat(outdir,'\LapFRcorr_',RatID{rr},'_',num2str(nFigures)),'fig')
            saveas(h,strcat(outdir,'\LapFRcorr_',RatID{rr},'_',num2str(nFigures)),'png')
            saveas(h2,strcat(outdir,'\LapFRcorr_all_',RatID{rr},'_',num2str(nFigures)),'fig')
            saveas(h2,strcat(outdir,'\LapFRcorr_all_',RatID{rr},'_',num2str(nFigures)),'png')
            close([h h2])
            nSubplot = 0;
        end 
    end
end

%% Plot correlation matrix averaged by rats
RatID = fieldnames(CorrMat);
h = figure;
set(h,'OuterPosition',[93,400,1738,673])
for rr = 1:length(RatID)
    DateID = fieldnames(CorrMat.(RatID{rr}));
    corrMat_all = CorrMat.(RatID{rr}).(DateID{1});    
    for ii = 2:length(DateID)
        corrMat_all= corrMat_all+CorrMat.(RatID{rr}).(DateID{ii});
    end
    corrMat_all = corrMat_all/length(DateID);
    corrMat_all = corrMat_all-diag(NaN(length(corrMat_all),1)); % remove diagonal elements for visualization
    imAlpha=ones(size(corrMat_all));
    imAlpha(isnan(corrMat_all))=0;
    
    subplot(1,length(RatID),rr)
    img = imagesc(corrMat_all);
    set(img,'AlphaData',imAlpha)
    axis square xy
    colorbar;
    title(RatID{rr},'Interpreter','none')
    xlabel('Lap number')
    if rr == 1
        ylabel('Lap number')
    end
end
saveas(h,strcat(outdir,'\LapFRcorr_all_averaged'),'epsc')
saveas(h,strcat(outdir,'\LapFRcorr_all_averaged'),'png')
close(h)

%% Plot place cell emergence
RatID = fieldnames(LapFRcorr);
for rr = 1:length(RatID)
    DateID = fieldnames(LapFRcorr.(RatID{rr}));
    for ii = 1:length(DateID)   
        condition = fieldnames(LapFRcorr.(RatID{rr}).(DateID{ii}));
        % sort the cells using position tunning in sample, test, and
        % Postrunning
        spkposFR=[];
        condition_selected = {'Sample','Test','Postrunning'};
        for c2 = 1:length(condition)
            if sum(strcmp(condition{c2},condition_selected))>0
                spkposFR = cat(3,spkposFR,LapFRcorr.(RatID{rr}).(DateID{ii}).(condition{cc}).lapFR);
            end
        end
        rmv = max(squeeze(mean(spkposFR,3)),[],2) < 1;
        [~,PeakFRloc] = max(mean(spkposFR,3),[],2);
        [~,sortind] = sort(PeakFRloc);
        rmv = rmv(sortind);
        
        h = figure;
        set(h,'OuterPosition',[69,187,1806,823])
        for cc = 1:length(condition)
            spkposFR = LapFRcorr.(RatID{rr}).(DateID{ii}).(condition{cc}).lapFR;
            posaxis = LapFRcorr.(RatID{rr}).(DateID{ii}).(condition{cc}).posbinaxis;

            % sort and remove inactive cells
            spkposFR = spkposFR(sortind,:,:);
            spkposFR(rmv,:,:)=[];

            % Linearize the firing rate matrix of cell population
            Mat = zeros(size(spkposFR,1),size(spkposFR,2)*size(spkposFR,3));
            for ss = 1:size(spkposFR,1)
                temp = squeeze(spkposFR(ss,:,:));
                Mat(ss,:) = temp(:);
            end
            Mat = Mat./repmat(max(Mat,[],2),1,size(Mat,2));

            subplot(length(condition),1,cc)
            imagesc(1:size(Mat,2),1:size(Mat,1),Mat);axis tight;
            %surf(1:size(Mat,2),1:size(Mat,1),zeros(size(Mat)),Mat);view([0 90]);axis tight;shading interp;grid off
            colormap(rvgray)
            hold on; 
            for pp = 1:size(spkposFR,3)
                plot([length(posaxis)*(pp-1)+1 length(posaxis)*(pp-1)+1],ylim,'r--'); 
            end
            axis xy
            set(gca,'Xtick',length(posaxis):5*length(posaxis):length(posaxis)*size(spkposFR,3),'Xticklabel',1:5:size(spkposFR,3))
            if cc == length(condition)
                xlabel('Lap#')
            end
            ylabel('Cell ID')
            c=colorbar;
            set(c,'Ytick',[0 1],'Yticklabel',{'0','Max'})
            ylabel(c,'Firing rate')
            title(condition{cc},'Interpreter','none')
            
        end
        set(h,'Renderer','painters')
        saveas(h,strcat(outdir,'\PlaceCellAcrossLap_',RatID{rr},'_',DateID{ii}),'fig')
        saveas(h,strcat(outdir,'\PlaceCellAcrossLap_',RatID{rr},'_',DateID{ii}),'epsc')
        close(h)
    end
end

function [LapFR,LapSpkind,cellID,posbinaxis,active,CorrectInd,corrMat]...
    = PlaceCellPropertiesAcrossLap_v4(spkfile,spikeIndfile,eeg2posIndfile,posIndfile,Condition,trackInfo)

%% Example inputs
% spkfile = strcat(parentfd,'\TTList_dCA3.txt');
% spikeIndfile = strcat(parentfd,'\SpikeSpect\spike_index.mat');
% eeg2posIndfile = strcat(parentfd,'\SpikeSpect\eeg2pos.mat');
% posIndfile = strcat(parentfd,'\Data_angle_ontrack.mat');

%% default parameters
posbinaxis = 0:4:356;
run_thr = 5; % in cm/sec
sampfreq = 30;
kernel = gausskernel(5,2);
%% Load and read necessary inputs
temp = load(spikeIndfile);
SpikeInd = temp.spike_ind;

cellID = Readtextfile(spkfile);
for i2 = 1:length(cellID)
    [~,cellID{i2}] = fileparts(cellID{i2});
end

load(eeg2posIndfile);
load(posIndfile);

%% Make sure type of session match previously saved behavioral data
if length(Condition) ~= size(data_angle,2)
    error('Session type mismatch')
end

%% Compute spatial tuning for each cell in each lap
for rr = 1:size(data_angle,2)
    if strcmp(Condition{rr},'Sample')
        correctTrial = trackInfo.sign_correct_sample;
    elseif strcmp(Condition{rr},'Test')
        correctTrial = trackInfo.sign_correct_test;
    elseif strcmp(Condition{rr},'Postrunning')
        correctTrial = trackInfo.sign_correct_posttest;
    else
        correctTrial=[];
    end
    trial_data = data_angle{rr};
    lapN = size(trial_data,1);
    LapFR.(Condition{rr}) = zeros(length(cellID),length(posbinaxis),lapN);
    LapSpkind.(Condition{rr}) = cell(length(cellID),lapN);
    CorrectInd.(Condition{rr}) = NaN(lapN,1);
    for ll = 1:lapN
        lap_ts = trial_data{ll,1}(:,1);
        % position index
        ind_range = [find((data_angle_all(:,1)-lap_ts(1))== 0),...
                     find((data_angle_all(:,1)-lap_ts(end))== 0)];
        % translate into eeg index
        eegind_range = [find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1),...
                        find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1,'last')];
        eegind = eeg2pos(eegind_range(1):eegind_range(2));
        lap_pos = data_angle_all(ind_range(1):ind_range(2),2);
        lap_speed = data_angle_all(ind_range(1):ind_range(2),3);
        lap_pos_eeg = data_angle_all(eegind,2);

        run_ind = lap_speed>=run_thr;
        run_ind_eeg = data_angle_all(eegind,3)>=run_thr;
        posaxis = [posbinaxis 360];
        occupancy = histcounts(lap_pos(run_ind),posaxis/180*pi);
        occupancy = conv_cir(occupancy,kernel)';
        occupancy = occupancy+.0001; % Add offset to prevent zeros
        
        for ii = 1:length(cellID)
            spkind_log = false(length(eegind),1);
            if isfield(SpikeInd,cellID{ii})
                spkind = SpikeInd.(cellID{ii});
                spkind = spkind(spkind>=min(eegind_range) & spkind<=max(eegind_range))-min(eegind_range)+1;
                spkind_log(spkind) = true;
            end
            if sum(run_ind_eeg & spkind_log) > 5
                spkposcount = histcounts(lap_pos_eeg(run_ind_eeg & spkind_log),posaxis/180*pi);
                LapSpkind.(Condition{rr}){ii,ll} = lap_pos_eeg(run_ind_eeg & spkind_log);
            else
                spkposcount = zeros(size(occupancy));
            end
            LapFR.(Condition{rr})(ii,:,ll) = conv_cir(spkposcount./occupancy,kernel)*sampfreq;            
        end
        if ~isempty(correctTrial)
            CorrectInd.(Condition{rr})(ll) = correctTrial(ll);
        end
    end
    peakFR = nanmax(nanmean(LapFR.(Condition{rr}),3),[],2);
    active.(Condition{rr}) = peakFR >= 1;
    corrMat.(Condition{rr}) = PopVecCorr(LapFR.(Condition{rr})(active.(Condition{rr}),:,:));
end

function corrMat = PopVecCorr(LapFR)
LapFR = LapFR+1e-10; % add small offset to prevent zeros
% normalize the firing rate vector by maximal firing rate
LapFR = LapFR./repmat(max(max(LapFR,[],3),[],2),1,size(LapFR,2),size(LapFR,3));
% normalize the firing rate vector by total firing rate
%LapFR = LapFR./repmat(sum(LapFR,2),1,size(LapFR,2),1);
popVec = zeros(size(LapFR,1)*size(LapFR,2),size(LapFR,3));
for ii = 1:size(LapFR,3)
    temp = LapFR(:,:,ii);
    popVec(:,ii) = temp(:);
end
corrMat = corr(popVec);