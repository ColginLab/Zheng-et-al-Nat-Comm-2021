% Parentfd = {'C:\CA1Replayproject\Data\Rat139\circulartrack\2017-02-20-CT-2';...
%             'C:\CA1Replayproject\Data\Rat139\circulartrack\2017-02-21-CT-1'};
% Outdir = 'C:\CA1Replayproject\Analysis';
function RunConfusionMatrix_AcrossTrials_v2(Parentfd,trackdata,Outdir,cellList,spktime_fd)

outdir = strcat(Outdir,'\ConfusionMatrix');
if ~isdir(outdir)
    mkdir(outdir)
end
for ii = 1:length(Parentfd)
    
    % Extract ratID and dateID
    ind = strfind(Parentfd{ii},'Rat');
    RatID = Parentfd{ii}(ind:ind+5);
    ind = strfind(Parentfd{ii},'\');
    DateID = Parentfd{ii}(ind(end-1)+1:end-1);
    ind = strfind(DateID,'-');
    DateID(ind)=[];
    DateID = strcat('D',DateID);

    spikeIndfile = strcat(Parentfd{ii},'\',spktime_fd,'\spike_index.mat');
    eeg2posIndfile = strcat(Parentfd{ii},'\',spktime_fd,'\eeg2pos.mat');
    posIndfile = strcat(Parentfd{ii},'\Data_angle_ontrack.mat');
    spkfile = strcat(Parentfd{ii},'\',cellList);
    % Define encoding and decoding folders
    decodingfd = [1];
    encodingfd = [2 3 4];
    [confusionMatrix,error,cellID,posaxis] = ...
        ConfusionMatrix_AcrossTrials(Parentfd{ii},decodingfd,encodingfd,posIndfile,eeg2posIndfile,spikeIndfile,spkfile);
    
    % Obtain reward position
    load(trackdata{ii});
    sign_correct_sample(isnan(sign_correct_sample)) = false;
    rewardPos = mode(ang_sample_reward_ontrack(logical(sign_correct_sample)));
    
    ConfusionMat.(RatID).(DateID).confusionMatrix = confusionMatrix;
    ConfusionMat.(RatID).(DateID).error = error;
    ConfusionMat.(RatID).(DateID).cellID = cellID;
    ConfusionMat.(RatID).(DateID).posaxis = posaxis;
    ConfusionMat.(RatID).(DateID).rewardPos = rewardPos;
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
save(strcat(outdir,'\ConfusionMat.mat'),'ConfusionMat','FunctionName')

%% Plot results
Rat = fieldnames(ConfusionMat);
colorLabel = lines;
posbinsize = 4; % in degree
CumError = []; % output from plotting cumulative error graphs
for rr = 1:length(Rat)
    dateID = fieldnames(ConfusionMat.(Rat{rr}));
    posaxis = ConfusionMat.(Rat{rr}).(dateID{1}).posaxis*180/pi;
    nFigures = 0;
    nSubplot = 0;    
    for dd = 1:length(dateID)
        nSubplot = nSubplot+1;
        if nSubplot == 1
            nFigures = nFigures+1;
            h = figure;
            set(h,'OuterPosition',[110,269,1615,705])
            ha = zeros(5,2);
        end

        confusionMatrix = ConfusionMat.(Rat{rr}).(dateID{dd}).confusionMatrix;
        rewardPos = ConfusionMat.(Rat{rr}).(dateID{dd}).rewardPos*180/pi;
        error = ConfusionMat.(Rat{rr}).(dateID{dd}).error*180/pi;
        
        ha(nSubplot,1) = subplot(2,5,nSubplot); hold on
        imagesc(posaxis,posaxis,confusionMatrix,[0 .4]); axis xy square tight
        plot(xlim,[rewardPos, rewardPos],'r--')
        colormap(rvgray)
        xlabel('Actual position (deg)')
        if nSubplot == 1
            ylabel('Decoded position (deg)')
        elseif nSubplot == 5
            pos = get(ha(nSubplot,1),'Position');
            hb = colorbar;
            pos_bar = get(hb,'Position'); pos_bar(1) = pos_bar(1)+.05;
            ylabel(hb,'Probability density')
            set(hb,'Position',pos_bar);
            set(ha(nSubplot,1),'Position',pos);
        end
        title(dateID{dd},'Interpreter','none')
        
        ha(nSubplot,2) = subplot(2,5,nSubplot+5);hold on
        erroraxis = 0:posbinsize:180;
        count = hist(error,erroraxis);
        if ~isfield(CumError,'DecodeError')
            CumError.DecodeError = cumsum(count)/sum(count);
        else
            CumError.DecodeError = cat(1,CumError.DecodeError,cumsum(count)/sum(count));
        end
        plot(erroraxis,cumsum(count)/sum(count),'b'); axis square
        xlim([0 180])
        xlabel('Error (deg)')
        if nSubplot == 1
            ylabel('Cumulative fraction')
        end
        
       % cumulative graph for normalized confusion matrix
        confusionMatrix(isnan(confusionMatrix))=1; % place numerical values in null decoded position 
        confusionMatrix = confusionMatrix./repmat(sum(confusionMatrix),size(confusionMatrix,1),1);
        errorIndex = -ceil(size(confusionMatrix,1)/2)+1:floor(size(confusionMatrix,1)/2);
        errorMat = repmat(errorIndex',1,size(confusionMatrix,1));
        for cc = 1:length(errorIndex)
            errorMat(:,cc) = circshift(errorMat(:,cc),[errorIndex(cc),0]);
        end

        errorIndex = errorIndex(errorIndex>=0); % Remove the sign; only consider error distance
        errorDist = zeros(1,length(errorIndex));
        for ee = 1:length(errorIndex)
            errorDist(ee) = sum(sum(confusionMatrix(abs(errorMat)==errorIndex(ee))));
        end
        errorDist_cum = cumsum(errorDist)/sum(errorDist);
        plot(erroraxis,errorDist_cum,'k')
        if ~isfield(CumError,'ConfusionMatError') && ~isfield(CumError,'RatID') && ~isfield(CumError,'DateID')
            CumError.ConfusionMatError = errorDist_cum;
            CumError.RatID = Rat(rr);
            CumError.DateID = dateID(dd);
        else
            CumError.ConfusionMatError = cat(1,CumError.ConfusionMatError,errorDist_cum);
            CumError.RatID = cat(1,CumError.RatID,Rat(rr));
            CumError.DateID = cat(1,CumError.DateID,dateID(dd));
        end
        CumError.errorAxis = erroraxis;
        
        if nSubplot == 5 || dd == length(dateID)
            saveas(h,strcat(outdir,'\ConfusionMat_',Rat{rr},'_',num2str(nFigures)),'png')
            saveas(h,strcat(outdir,'\ConfusionMat_',Rat{rr},'_',num2str(nFigures)),'epsc')
            close(h)
            nSubplot = 0;
        end
    end
end
save(strcat(outdir,'\CumError.mat'),'CumError')

function [confusionMatrix,errorDis,cellID,posaxis] = ...
    ConfusionMatrix_AcrossTrials(parentfd,decodingfd,encodingfd,posIndfile,eeg2posIndfile,spikeIndfile,spklist)
% Example inputs
% parentfd = 'C:\CA1Replayproject\Data\Rat139\circulartrack\2017-02-20-CT-2';
% decodingfd = [2 3 4];
% encodingfd = [1]; (1 is Prerunning; 2 is sample; 3 is test; 4 is TrialsPost)
% posIndfile = strcat(parentfd,'\Data_angle_ontrack.mat');
% eeg2posIndfile = strcat(parentfd,'\SpikeTime\eeg2pos.mat');
% spikeIndfile = strcat(parentfd,'\SpikeTime\spike_index.mat');
% spklist = strcat(parentfd,'\TTList_dCA1_pyr.txt');

%parameters
run_thr = 5; % in cm/s
posbinsize = 4; % in degree
baysianTimeWindow = .5; %in sec
step = baysianTimeWindow/5;
FRkernel = gausskernel(5,2);
FRoffset = 0.0001; % in Hz
sampFreq = 2000; % eeg sampling frequency


%% Extract cells ID
cellID = Readtextfile(spklist);
[~,cellID_noExt] = cellfun(@fileparts,cellID,'UniformOutput',false);

%% Obtain spike index
temp = load(spikeIndfile);
spkind = temp.spike_ind;
id = fieldnames(spkind); % cell ID from spike index
if length(id) ~= length(cellID_noExt)
    error('# of cells does not match')
end
chk = strcmp(id,cellID_noExt); % check if cell IDs match
if sum(chk) ~= length(chk)
    error('cell IDs do not match')
end


%% Obtain firing rate from encoding trial
[~,spkFR] = GenerateSpkMatrix_v6(eeg2posIndfile,posIndfile,posbinsize,false,spikeIndfile,encodingfd);

%% Obtain spike raster from all trials
[SpikeMatrix,~,posaxis,~,Posbin_ind] = GenerateSpkMatrix_v6(eeg2posIndfile,posIndfile,posbinsize,false,spikeIndfile);
posaxis = circ_ang2rad(posaxis);

%Remove cells with peak firing rate less than 1 Hz
rmv = max(spkFR,[],2)<1;
spkFR = spkFR(~rmv,:);

smspkFR = zeros(size(spkFR));
for ss = 1:size(smspkFR,1)
    smspkFR(ss,:) = conv_cir(spkFR(ss,:),FRkernel);
end
smspkFR = smspkFR+FRoffset; %add offset to prevent 0 value

confusionMatrix = zeros(size(posaxis,2));
PosCount = zeros(1,size(posaxis,2));
errorDis = [];
load(eeg2posIndfile);
load(posIndfile);
for t = 1:length(decodingfd)
    trial_data = data_angle{decodingfd(t)};
    lapN = size(trial_data,1);
    truepos_all = data_angle_all(:,2);
    speed_all = data_angle_all(:,3);
    for ll = 1:lapN
        include = false(size(eeg2pos));
        lap_ts = trial_data{ll,1}(:,1);

        % position index
        ind_range = [find((data_angle_all(:,1)-lap_ts(1))== 0),...
                     find((data_angle_all(:,1)-lap_ts(end))== 0)];
        % translate into eeg index
        eegind_range = [find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1),...
                        find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1,'last')];
        include(min(eegind_range):max(eegind_range)) = true;  
       
        SpikeMatrix_lap = full(SpikeMatrix(~rmv,include));
        truepos = truepos_all(eeg2pos(include));
        speed = speed_all(eeg2pos(include));
        truepos_avg = zeros(1,ceil(sum(include)/sampFreq/step));
        for tt = 1:size(truepos_avg,2)-1
            range = round((tt-1)*step*sampFreq)+1:round((tt)*step*sampFreq);
            if mean(speed(range)) >= run_thr
                truepos_avg(1,tt) = circ_mean(truepos(range))+2*pi;
            else
                truepos_avg(1,tt)=NaN;
            end
        end
        truepos_avg(1,end) = circ_mean(truepos(max(range):end))+2*pi;
        truepos_avg(1,truepos_avg>2*pi) = truepos_avg(1,truepos_avg>2*pi)-2*pi;
        ind_nan = isnan(truepos_avg);
        truepos_avg = match(truepos_avg,posaxis)';
        truepos_avg(ind_nan)=NaN;
    
        p_x_n = BayesianDecoder(SpikeMatrix_lap,smspkFR,baysianTimeWindow,step,2000);    
        [~,dpos] = max(p_x_n); 
        dpos = posaxis(dpos);
        dpos(isnan(sum(p_x_n))) = NaN;
        if size(truepos_avg,2) > size(p_x_n,2)
            truepos_avg(end)=[];
        end
        errorDis = [errorDis abs(dpos(~isnan(dpos) & ~isnan(truepos_avg))-posaxis(truepos_avg(~isnan(dpos)& ~isnan(truepos_avg))))];    
        for ii = 1:size(confusionMatrix,1)
            confusionMatrix(:,ii) = confusionMatrix(:,ii)+nansum(p_x_n(:,truepos_avg==ii),2);
            PosCount(1,ii) = PosCount(1,ii)+nansum(truepos_avg==ii);
        end
    end
end
confusionMatrix = confusionMatrix./repmat(PosCount,size(confusionMatrix,1),1);
errorDis(errorDis>pi) = 2*pi-errorDis(errorDis>pi);