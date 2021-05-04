function CellPairInRipple(DataFD,AnalysisFD)

if exist(strcat(AnalysisFD,'\LapFRcorr\LapFRcorr.mat'),'file') ~= 2
    error('Please run RunPlaceCellPropertiesAcrossLap_v4.m first')
else
    load(strcat(AnalysisFD,'\LapFRcorr\LapFRcorr.mat'))
end

if exist(strcat(AnalysisFD,'\BayesianDecoding_ripple\BayesianDecodingResult.mat'),'file') ~= 2
    error('Please run RunBayesianDecoder_ripple.m first')
else
    load(strcat(AnalysisFD,'\BayesianDecoding_ripple\BayesianDecodingResult.mat'))
end

outdir = strcat(AnalysisFD,'\CellPairInRipple');
if ~isdir(outdir)
    mkdir(outdir)
end

Cellpair_label = [];
Xcorr_all = [];
PFS_all = [];
ratID = fieldnames(LapFRcorr);
for rr = 1:length(ratID)
    dateID = fieldnames(LapFRcorr.(ratID{rr}));
    DateFd = translateDateID(dateID);
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
        ind = match(rewardPos,trackdata.Ang_RewardLoc_ontrack);
        if ind == 1
            rewardPos_range = [rewardPos-(trackdata.Ang_RewardLoc_ontrack(ind+1)-rewardPos)/2 (trackdata.Ang_RewardLoc_ontrack(ind+1)+rewardPos)/2];
        elseif ind == length(trackdata.Ang_RewardLoc_ontrack)
            rewardPos_range = [(trackdata.Ang_RewardLoc_ontrack(ind-1)+rewardPos)/2 rewardPos+(rewardPos-trackdata.Ang_RewardLoc_ontrack(ind-1))/2];
        else
            rewardPos_range = [(trackdata.Ang_RewardLoc_ontrack(ind-1)+rewardPos)/2 (trackdata.Ang_RewardLoc_ontrack(ind+1)+rewardPos)/2];
        end
        % Obtain correct trials
        CorrectTest = find(trackdata.sign_correct_test==1);
        CorrectPostTest = find(trackdata.sign_correct_posttest==1);
        % Obtain position tunning
        FR_SampleTestRun = nanmean(cat(3,LapFRcorr.(ratID{rr}).(dateID{dd}).Sample.lapFR,...
                                         LapFRcorr.(ratID{rr}).(dateID{dd}).Test.lapFR,...
                                         LapFRcorr.(ratID{rr}).(dateID{dd}).Postrunning.lapFR),3);
        FR_PreRun = nanmean(LapFRcorr.(ratID{rr}).(dateID{dd}).Prerunning.lapFR,3);
        posaxis = LapFRcorr.(ratID{rr}).(dateID{dd}).Sample.posbinaxis;
        cellID = LapFRcorr.(ratID{rr}).(dateID{dd}).Sample.cellID;
        if ~isfield(BayesianDecodingResult.(ratID{rr}),dateID{dd})
            continue
        end
        [cellID,ind1,ind2] = intersect(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).cellID,cellID);
        % Select cells with place fields near reward position
        ind_reward = match(rewardPos_range,posaxis/180*pi);
        ActiveNearReward_preRun = mean(FR_PreRun(ind2,ind_reward(1):ind_reward(2)),2) >= 1;%mean(mean(FR_PreRun(ind2,ind(1):ind(2))));
        ActiveNearReward_SampleTestRun = mean(FR_SampleTestRun(ind2,ind_reward(1):ind_reward(2)),2) >= 1;%mean(mean(FR_SampleTestRun(ind2,ind(1):ind(2))));      
        % Select cells with place fields near rest box
        restPos_range = [trackdata.Ang_StartZone_depart_ontrack trackdata.Ang_StartZone_arrive_ontrack];
        ind_rest = match(restPos_range,posaxis/360*2*pi);
        logical_ind = false(size(posaxis));
        logical_ind(1:ind_rest(1)) = true; logical_ind(ind_rest(2):end) = true; % rest box located near 0 on circular track
        ActiveNearRest= mean(FR_PreRun(ind2,logical_ind),2) >= 1;%mean(mean(FR_PreRun(ind2,logical_ind)));
        ActiveOnTrack = (max(FR_PreRun(ind2,:),[],2) >= 1) | (max(FR_SampleTestRun(ind2,:),[],2) >= 1);
        
        sessionType = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sessionType;
        sType = unique(sessionType);
        %sType(strcmp(sType,'PreRun'))=[];
        sType = sType([4; 1; 2; 3]); %re-arrange the order
        Ncell = sum(ActiveNearReward_preRun)+...
                sum(~ActiveNearReward_preRun & ActiveNearReward_SampleTestRun)+...
                sum(~ActiveNearReward_preRun & ~ActiveNearReward_SampleTestRun & ActiveOnTrack);
        cellInd_in = [find(ActiveNearReward_preRun);... 
                      find(~ActiveNearReward_preRun & ActiveNearReward_SampleTestRun);...
                      find(~ActiveNearReward_preRun & ~ActiveNearReward_SampleTestRun & ActiveOnTrack)];
        for ss = 1:length(sType) % loop through each session type
            in = strcmp(sessionType,sType{ss});
            spkind = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).spkind(in);
            duration = (BayesianDecodingResult.(ratID{rr}).(dateID{dd}).RippleOffsetIndex(in)-...
                        BayesianDecodingResult.(ratID{rr}).(dateID{dd}).RippleOnsetIndex(in))/2000; % in sec
            XcorrMat.(sType{ss}) = NaN(length(spkind),Ncell);
            for s2 = 1:length(spkind) % loop through each ripple events
                cell_ind = ind1(spkind{s2}(:,1));
                temp = num2cell(1:Ncell)';
                id = cellfun(@num2str,temp,'UniformOutput',false);
                ind = categorical(cell_ind,cellInd_in,id);
                counts = histcounts(ind,id);
                XcorrMat.(sType{ss})(s2,:) = counts/duration(s2);
            end
            XcorrMat.(sType{ss}) = corr(XcorrMat.(sType{ss}),'Row','pairwise');
        end
        % Compute place field similarity
        positionTuning = FR_SampleTestRun(ind2,:);
        positionTuning = positionTuning(cellInd_in,:);
        pfs = corr(positionTuning','Row','pairwise');
        % Create cell pair labels
        temp = [ones(sum(ActiveNearReward_preRun),1);... 
                ones(sum(~ActiveNearReward_preRun & ActiveNearReward_SampleTestRun),1)+1;...
                ones(sum(~ActiveNearReward_preRun & ~ActiveNearReward_SampleTestRun & ActiveOnTrack),1)+3];
        cellpair_ind = repmat(temp,1,Ncell)+repmat(temp',Ncell,1);
        cellpair_label = cell(size(cellpair_ind));
        cellpair_label(cellpair_ind==2) = {'Pre-pre'};
        cellpair_label(cellpair_ind==3) = {'Pre-post'};
        cellpair_label(cellpair_ind==4) = {'Post-post'};
        cellpair_label(cellpair_ind>4) = {'Other'};
               
        % Select one cell pair using lower triangle of the matrix
        linearInd = find(tril(true(size(cellpair_label)),-1));
        cellpair_label = cellpair_label(linearInd);
        PFS = pfs(linearInd);
        Xcorr = [];
        for ss = 1:length(sType) % loop through each session type
            Xcorr = cat(2,Xcorr,XcorrMat.(sType{ss})(linearInd));
        end
        
        Cellpair_label = cat(1,Cellpair_label,cellpair_label);
        Xcorr_all = cat(1,Xcorr_all,Xcorr);
        PFS_all = cat(1,PFS_all,PFS);
    end
end

%% Plot results
h = figure;
set(h,'OuterPosition',[785,201,574.4,602.4])
label = unique(Cellpair_label);
color_label = {'k-','r-','g-','b-'};
for ii = 1:length(label)
    in = strcmp(Cellpair_label,label{ii});
    CI = bootci(5000,@nanmean,Xcorr_all(in,:));
    CI(1,:) = nanmean(Xcorr_all(in,:))-CI(1,:);
    CI(2,:) = CI(2,:)-nanmean(Xcorr_all(in,:));
    
    errorbar(1:size(Xcorr_all,2),nanmean(Xcorr_all(in,:)),CI(1,:),CI(2,:),color_label{ii});
    hold on
end
xlim([.5 4.5])
axis square
set(gca,'XTick',1:length(sType),'XTickLabel',sType,'Box','off')
legend(label)
ylabel('SWR cofiring (r)')
saveas(h,strcat(outdir,'\CellPairCorr'),'epsc')
saveas(h,strcat(outdir,'\CellPairCorr'),'png')
close(h)

%% Repeated-measured ANOVA
TimePoint = Xcorr_all;
t=[cell2table(Cellpair_label) array2table(PFS_all) array2table(TimePoint)];
TimePoint = table((1:size(TimePoint,2))','VariableNames',{'TimePoint'});
rm = fitrm(t,['TimePoint1-TimePoint',num2str(size(TimePoint,1)),'~Cellpair_label*PFS_all'],'WithinDesign',TimePoint);
Stats.ANOVA = ranova(rm);
Stats.posthoc = multcompare(rm,'TimePoint');

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
save(strcat(outdir,'\Results.mat'),'Xcorr_all','Cellpair_label','FunctionName')
save(strcat(outdir,'\Stats.mat'),'Stats')

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