function RunBayesianDecoder_ripple_v3(Path,CSClist,Outdir,cellList,spktime_fd,trackdata,plotOnly)
%% This version detect ripples outside of rest box
outdir = strcat(Outdir,'\BayesianDecoding_ripple_outsideOfRestBox_approachReward');
if ~isdir(outdir)
    mkdir(outdir)
end
if exist(strcat(outdir,'\BayesianDecodingResult.mat'),'file')==2 && plotOnly
    load(strcat(outdir,'\BayesianDecodingResult.mat'))
else
    for ii = 1:length(Path)

        % Extract ratID and dateID
        ind = strfind(Path{ii},'Rat');
        RatID = Path{ii}(ind:ind+5);
        ind = strfind(Path{ii},'\');
        DateID = Path{ii}(ind(end-1)+1:end-1);
        ind = strfind(DateID,'-');
        DateID(ind)=[];
        DateID = strcat('D',DateID);

        spikeIndfile = strcat(Path{ii},'\',spktime_fd,'\spike_index.mat');
        eeg2posIndfile = strcat(Path{ii},'\',spktime_fd,'\eeg2pos.mat');
        posIndfile = strcat(Path{ii},'\Data_angle_ontrack.mat');
        behavioralData = load(trackdata{ii});
        spkfile = strcat(Path{ii},'\',cellList);

        sampfreq = 2000;
        n_std = 5; %number of standard deviations for the threshold
        duration_threshold = .05; %in sec. duration of ripple event
        CA1eeg = {};
        for cc = 1:length(CSClist{ii})
            CSClist_CA1 = strcat(Path{ii},'\CSC',num2str(CSClist{ii}(cc)),'.ncs');
            [eeg,eegt,sampfreq] = LoadEEG(CSClist_CA1);
            CA1eeg = cat(1,CA1eeg,eeg);
        end
        if isempty(CA1eeg)
            continue
        end
        [RippleOnsetIndex, RippleOffsetIndex] = DetectRipples_v4(CA1eeg,sampfreq,n_std,duration_threshold);

        % Define encoding folders
        encodingfd = [2 3 4];
        %encodingfd = [1];
        [CurrentLoc,p_x_n,spkind,r_squared,calphase,slope,sessionType,NUniqueCell,NTrial,NTrial_post,NTrial_preRun,...
            RippleOnsetIndex,RippleOffsetIndex,beforeReward,RippleFR,ActiveNearReward,ActiveNearRest,RippleFR_xlabel,smspkFR,cellID,param] = ...    
        BayesianDecoding_v2(spkfile,encodingfd,posIndfile,eeg2posIndfile,spikeIndfile,behavioralData,RippleOnsetIndex,RippleOffsetIndex,eegt);

        BayesianDecodingResult.(RatID).(DateID).CurrentLoc = CurrentLoc;
        BayesianDecodingResult.(RatID).(DateID).p_x_n = p_x_n;
        BayesianDecodingResult.(RatID).(DateID).spkind = spkind;
        BayesianDecodingResult.(RatID).(DateID).r_squared = r_squared;
        BayesianDecodingResult.(RatID).(DateID).calphase = calphase;
        BayesianDecodingResult.(RatID).(DateID).slope = slope;
        BayesianDecodingResult.(RatID).(DateID).sessionType = sessionType;
        BayesianDecodingResult.(RatID).(DateID).NUniqueCell = NUniqueCell;
        BayesianDecodingResult.(RatID).(DateID).NTrial = NTrial; % [sample test]
        BayesianDecodingResult.(RatID).(DateID).NTrial_post = NTrial_post;
        BayesianDecodingResult.(RatID).(DateID).NTrial_preRun = NTrial_preRun;
        BayesianDecodingResult.(RatID).(DateID).RippleOnsetIndex = RippleOnsetIndex;
        BayesianDecodingResult.(RatID).(DateID).RippleOffsetIndex = RippleOffsetIndex;
        BayesianDecodingResult.(RatID).(DateID).beforeReward = beforeReward;
        BayesianDecodingResult.(RatID).(DateID).RippleFR = RippleFR;
        BayesianDecodingResult.(RatID).(DateID).PosTuning = smspkFR;
        BayesianDecodingResult.(RatID).(DateID).cellID = cellID;
        BayesianDecodingResult.(RatID).(DateID).ActiveNearReward = ActiveNearReward;
        BayesianDecodingResult.(RatID).(DateID).ActiveNearRest = ActiveNearRest;
        BayesianDecodingResult.(RatID).(DateID).RippleFR_xlabel = RippleFR_xlabel;
        BayesianDecodingResult.(RatID).(DateID).ang_test_reward_ontrack = behavioralData.ang_test_reward_ontrack;
        BayesianDecodingResult.(RatID).(DateID).ang_sample_reward_ontrack = behavioralData.ang_sample_reward_ontrack;
        BayesianDecodingResult.(RatID).(DateID).ang_posttest_reward_ontrack = behavioralData.ang_posttest_reward_ontrack;
        BayesianDecodingResult.(RatID).(DateID).sign_correct_test = behavioralData.sign_correct_test;
        BayesianDecodingResult.(RatID).(DateID).sign_correct_sample = behavioralData.sign_correct_sample;
        BayesianDecodingResult.(RatID).(DateID).sign_correct_posttest = behavioralData.sign_correct_posttest;
        BayesianDecodingResult.(RatID).(DateID).param = param;
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
    save(strcat(outdir,'\BayesianDecodingResult.mat'),'BayesianDecodingResult','FunctionName')
end

trial_ext = {'sample','test'};
for tt = 1:length(trial_ext)
%% Extract relevant data for plotting
p_x_n = {};
duration = [];
r_squared = [];
slope = [];
sessionType = {};
Correct_test = [];
stopLoc = [];
rewardLoc = [];
currentLoc = [];
beforeReward = [];
NUniqueCell_thr = 0;
NTrial_all = [];
ratID = fieldnames(BayesianDecodingResult);
for rr = 1:length(ratID)
    dateID = fieldnames(BayesianDecodingResult.(ratID{rr}));
    for dd = 1:length(dateID)
        correct_sample = ~isnan(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sign_correct_sample);
        correct_sample = find(correct_sample);
        NUniqueCell = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).NUniqueCell;
        if strcmp(trial_ext{tt},'sample')
            NTrial = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).NTrial(:,1); % 1 is sample; 2 is test
        elseif strcmp(trial_ext{tt},'test')
            NTrial = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).NTrial(:,2); % 1 is sample; 2 is test
        end
        NTrial_post = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).NTrial_post;
        NTrial_preRun = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).NTrial_preRun;
        rewardLoc_session = mode(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).ang_sample_reward_ontrack);
        
        % Remove leading and trailing time bins without spikes            
        pxn_all = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).p_x_n;
        for pp = 1:length(pxn_all)
            withspk = ~isnan(pxn_all{pp}(1,:));
            range = find(withspk,1):find(withspk,1,'last');
            pxn_all{pp} = pxn_all{pp}(:,range);
        end
        [~,dur] = cellfun(@size,pxn_all);
        dur = dur*BayesianDecodingResult.(ratID{rr}).(dateID{dd}).param.step; % in sec 
        
        if ~isempty(NTrial)
            include1 = NTrial~=0; % 0 corresponds to preRun and postRun trials 
            include2 = ~isnan(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared) & NUniqueCell >= NUniqueCell_thr;
            correct_test = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sign_correct_test;
            
            NTrial_all = cat(1,NTrial_all,NTrial(include1 & include2));
            p_x_n = cat(1,p_x_n,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).p_x_n(include1 & include2));
            duration = cat(1,duration,dur(include1 & include2));
            r_squared = cat(1,r_squared,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared(include1 & include2));
            slope = cat(1,slope,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).slope(include1 & include2));
            sessionType = cat(1,sessionType,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sessionType(include1 & include2));
            Correct_test = cat(1,Correct_test,correct_test(NTrial(include1 & include2)));
            stopLoc = cat(1,stopLoc,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).ang_test_reward_ontrack(NTrial(include1 & include2)));
            rewardLoc = cat(1,rewardLoc,repmat(rewardLoc_session,sum(include1 & include2),1));
            currentLoc = cat(1,currentLoc,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).CurrentLoc(include1 & include2));
            beforeReward = cat(1,beforeReward,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).beforeReward(include1 & include2));
        end

        if ~isempty(NTrial_post)
            include3 = NTrial_post~=0 & ~isnan(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared) & NUniqueCell >= NUniqueCell_thr;
            correct_posttest = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sign_correct_posttest;
            
            NTrial_all = cat(1,NTrial_all,NaN(sum(include3),1));
            p_x_n = cat(1,p_x_n,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).p_x_n(include3));
            duration = cat(1,duration,dur(include3));
            r_squared = cat(1,r_squared,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared(include3));
            slope = cat(1,slope,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).slope(include3));
            sessionType = cat(1,sessionType,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sessionType(include3));
            Correct_test = cat(1,Correct_test,correct_posttest(NTrial_post(include3)));
            stopLoc = cat(1,stopLoc,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).ang_posttest_reward_ontrack(NTrial_post(include3)));
            rewardLoc = cat(1,rewardLoc,repmat(rewardLoc_session,sum(include3),1));
            currentLoc = cat(1,currentLoc,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).CurrentLoc(include3));
            beforeReward = cat(1,beforeReward,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).beforeReward(include3));
        end

        if ~isempty(NTrial_preRun)
            include4 = NTrial_preRun~=0 & ~isnan(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared) & NUniqueCell >= NUniqueCell_thr;
            
            NTrial_all = cat(1,NTrial_all,NaN(sum(include4),1));
            p_x_n = cat(1,p_x_n,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).p_x_n(include4));
            duration = cat(1,duration,dur(include4));
            r_squared = cat(1,r_squared,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared(include4));
            slope = cat(1,slope,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).slope(include4));
            sessionType = cat(1,sessionType,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sessionType(include4));
            Correct_test = cat(1,Correct_test,NaN(sum(include4),1));
            stopLoc = cat(1,stopLoc,NaN(sum(include4),1));
            rewardLoc = cat(1,rewardLoc,repmat(rewardLoc_session,sum(include4),1));
            currentLoc = cat(1,currentLoc,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).CurrentLoc(include4));
            beforeReward = cat(1,beforeReward,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).beforeReward(include4));
        end
        
        if ~isempty(NTrial) && ~isempty(NTrial_post) && ~isempty(NTrial_preRun)
            %% Plot average posterior probability aligned by stop or reward locations for each recording session
            in = {(include1 & include2);include3;include4};
            LineColor = {'r';'k';'b'};
            h2 = figure; hold on
            for ii = 1:length(in)
                pxn_subset = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).p_x_n(in{ii});
                posbinsize = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).param.posbinsize;
                posaxis = 0:posbinsize:2*pi-posbinsize;
                shift = match(rewardLoc_session,posaxis)-length(posaxis)/2;
                Ntimebin = cellfun(@size,pxn_subset,repmat({2},size(pxn_subset)));
                pxn_aligned = zeros(length(posaxis),sum(Ntimebin));
                n = 1;
                for pp = 1:length(pxn_subset)
                    n2 = n+Ntimebin(pp);
                    pxn_aligned(:,n:n2-1) = circshift(pxn_subset{pp},[-shift,0]);
                    n = n2;
                end
                pxn_aligned = nanmedian(pxn_aligned,2);
                plot(posaxis-pi,pxn_aligned,LineColor{ii});
            end
            legend('SampleTest','PostTest','PreRun')
            set(gca,'YScale','log')
            ylabel('Median probability')
            xlabel('Aligned angular position (rad)')
            axis square tight
            title([ratID{rr},' Date ',dateID{dd}])
            %saveas(h2,strcat(outdir,'\',ratID{rr},dateID{dd}),'png')
            close(h2)
        end
    end
end
posbinsize = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).param.posbinsize;
posaxis = 0:posbinsize:2*pi-posbinsize;
timeStep = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).param.step; % in sec

%% Plot average posterior probability aligned by stop or reward locations and replay onset or offset
r2_thr = 0.5;
dur_thr = .05; % in sec
emptyBin_thr = .2; % percent number of empty bin in pxn allowed
LocDis = 14; % in deg
nLocAway = 8;
nShuffle = 5000; % shuffle position to create null distribution
nType = unique(sessionType);
Replay_alignment = {'onset','offset'};
Position_alignment = {'current position','goal position'};
tbinL = 6; % number of time bins
CrtErr = {'All','Correct','Error'};
for cc = 1:length(CrtErr)  
    for r1 = 1:length(Replay_alignment)
        for p1 = 1:length(Position_alignment)
            emptyBin_percent = zeros(length(pxn_aligned),1);
            if strcmp(Position_alignment{p1},'current position')
                % align current position to the middle
                ind_half = match(pi,posaxis);
                ind = match(currentLoc,posaxis);
                shift = ind-ind_half;
                pxn_aligned = p_x_n;
                for ii = 1:length(pxn_aligned)
                    pxn_aligned{ii} = circshift(pxn_aligned{ii},-shift(ii));
                    emptyBin_percent(ii) = sum(isnan(pxn_aligned{ii}(1,:)))/length(pxn_aligned{ii}(1,:));
                end

                fname2 = 'currPos';
            elseif strcmp(Position_alignment{p1},'goal position')
                % align reward position to the middle
                ind = match(rewardLoc,posaxis);
                ind_half = match(pi,posaxis);
                shift = ind-ind_half;
                pxn_aligned = p_x_n;
                for ii = 1:length(pxn_aligned)
                    pxn_aligned{ii} = circshift(pxn_aligned{ii},-shift(ii));
                    emptyBin_percent(ii) = sum(isnan(pxn_aligned{ii}(1,:)))/length(pxn_aligned{ii}(1,:));
                end

                fname2 = 'goalPos';
            end

            % combine replay events in each type of session
            sessions = unique(sessionType);
            sessions = sessions([4;1;2;3]);

            % Construct bar graph for possible goal location
            locationRangeInd = match(-LocDis*nLocAway:LocDis:LocDis*nLocAway,(posaxis-pi)*180/pi);
            ProbSum = zeros(length(locationRangeInd),length(sessions));

            h = figure;
            set(h,'OuterPosition',[13.8,309.8,1511.2,551.2])
            ha = zeros(length(sessions),1);
            for ss = 1:length(sessions)
                in = strcmp(sessionType,sessions{ss}) &...
                    r_squared >= r2_thr &...
                    duration >= dur_thr &...
                    emptyBin_percent <= emptyBin_thr;                
                if ~strcmp(sessions{ss},'PreRun')
                    if strcmp(CrtErr{cc},'Correct')
                        in_err = in & Correct_test == 0;
                        in = in &...
                             Correct_test == 1; % correct or error
                             %(strcmp(sessionType2,'BeforeTest') | strcmp(sessionType2,'PostTest'))
%                         if sum(in) > sum(in_err)
%                             % randomly downsample replay events
%                             count_ds = sum(in)-sum(in_err);
%                             ind = find(in);
%                             rmv = randperm(s,length(ind),count_ds);
%                             ind(rmv)=[];
%                             in = false(size(Correct_test));
%                             in(ind) = true;
%                         end
                    elseif strcmp(CrtErr{cc},'Error')
                        in = in &...
                             Correct_test == 0; % correct or error
                             %(strcmp(sessionType2,'BeforeTest') | strcmp(sessionType2,'PostTest'))                    
                    end
                end
                pxn = pxn_aligned(in & beforeReward);
                taxis = [];
                %[~,tbinL] = cellfun(@size,pxn);
                pxn_combined = NaN(length(posaxis),min(tbinL),length(pxn));
                for pp = 1:length(pxn)
                    if strcmp(Replay_alignment{r1},'onset')
                        range = 1:min(tbinL);% align onset of replay

                        taxis = 0:timeStep:(min(tbinL)-1)*timeStep;
                        fname1 = 'InitiationBias';
                    elseif strcmp(Replay_alignment{r1},'offset')
                        range = size(pxn{pp},2)-min(tbinL)+1:size(pxn{pp},2);% align offset of replay

                        taxis = 0-(min(tbinL)-1)*timeStep:timeStep:0;
                        fname1 = 'TerminationBias';
                    end
                    pxn_combined(:,:,pp) = pxn{pp}(:,range);
                end

                ha(ss) = subplot(2,length(sessions),ss);
                imagesc(taxis,(posaxis-pi)*180/pi,nanmean(pxn_combined,3))
                axis xy square
                hold on
                plot(xlim,[0 0],'w--','LineWidth',1)
                cb = colorbar;
                title([sessions{ss},' n=',num2str(sum(in & beforeReward))])
                if ss == 1
                    xlabel(['Time relative to replay ', Replay_alignment{r1}, ' (sec)'])
                    ylabel(['Position relative to ', Position_alignment{p1}, ' (deg)'])
                elseif ss == length(sessions)
                    ylabel(cb,'Mean posterior probability')
                end

                % Add up probability
                for pp = 1:length(locationRangeInd)
                    inrange = locationRangeInd(pp)-1:locationRangeInd(pp)+1;
                    ProbSum(pp,ss) = nansum(nansum(nansum(pxn_combined(inrange,:,:))));
                end
                % Find null distribution
                ProbSum_null = zeros(length(locationRangeInd),nShuffle);
                for ff = 1:nShuffle
                    s = RandStream('mt19937ar','Seed',ff);
                    ind_shu = randi(s,size(pxn_combined,1),size(pxn_combined,3),1);
                    pxn_combined_shu = pxn_combined;
                    for s1 = 1:size(pxn_combined,3)
                        pxn_combined_shu(:,:,s1) = circshift(pxn_combined(:,:,s1),ind_shu(s1));
                    end
                    for pp = 1:length(locationRangeInd)
                        inrange = locationRangeInd(pp)-1:locationRangeInd(pp)+1;
                        ProbSum_null(pp,ff) = nansum(nansum(nansum(pxn_combined_shu(inrange,:,:))));
                    end
                end
                [CI_u, CI_l] = CorrectionForMultipleComparsion(ProbSum_null');
                CI = [CI_l CI_u];

                subplot(2,length(sessions),ss+length(sessions))
                bar(-nLocAway:nLocAway,ProbSum(:,ss)); hold on
                plot(-nLocAway:nLocAway,CI,'k--')
                set(gca,'Box','off','XTick',[-nLocAway, 0, nLocAway])
                axis square
                colorbar
                if ss == 1
                    xlabel(['Location # relative to ', Position_alignment{p1}])
                    ylabel('sum of posterior probability')
                end
            end
            colormap hot
            set(ha,'CLim',[0 .04])

            if exist('outdir','var')
                saveas(h,strcat(outdir,'\',fname1,'_',fname2,'_',CrtErr{cc},'_',trial_ext{tt}),'epsc')
                saveas(h,strcat(outdir,'\',fname1,'_',fname2,'_',CrtErr{cc},'_',trial_ext{tt}),'png')
                close(h)
            end
        end
    end
end

%% Examine changes in r2 values as a function of trial numbers of Sample/Test
h = figure;
plot(NTrial_all+randn(size(NTrial_all))*.05,r_squared,'k.')
hold on;
xTickL = cell(1,7);
for ii = 1:7
    in = NTrial_all==ii;
    CI = bootci(5000,@nanmean,r_squared(in));
    CI(1) = nanmean(r_squared(in))-CI(1);
    CI(2) = CI(2)-nanmean(r_squared(in));
    errorbar(ii,nanmean(r_squared(in)),CI(1),CI(2),'k','LineWidth',1)
    xTickL{ii} = [num2str(ii),'-',num2str(ii+1)];
end
axis square
xlim([.5 7.5])
set(gca,'XTick',1:7,'XTickLabel',xTickL,'Box','off')
ylabel('Replay fidelity (r2)')
xlabel('Sample/Test trials')
saveas(h,strcat(outdir,'\r2overTrials_',trial_ext{tt}),'epsc')
saveas(h,strcat(outdir,'\r2overTrials_',trial_ext{tt}),'png')
close(h)

% Stats for r2 values
in = NTrial_all<8; % exclude the last trials since the last sleep time in rest box is longer than other sleep time
[~,tbl,stats]=kruskalwallis(r_squared(in),NTrial_all(in),'off');
posthoc = multcompare(stats,'Display','off');
Stats.r2.KKW = tbl;
Stats.r2.posthoc = posthoc;
[r,p] = corr(NTrial_all(in),r_squared(in),'Row','pairwise','Type','Spearman');
tbl2 = table(r,p);
Stats.r2.corr = tbl2;
save(strcat(outdir,'\Stats.mat'),'Stats')

%% Make separate plots for each type of behavioral session
% sort the events by r square values 
[~,sortind] = sort(r_squared,'descend');
p_x_n = p_x_n(sortind);
r_squared = r_squared(sortind);
slope = slope(sortind);
sessionType = sessionType(sortind);
Correct_test = Correct_test(sortind);
stopLoc = stopLoc(sortind);
rewardLoc = rewardLoc(sortind);
% remove events with too many empty bin
emptyBinThr = true(size(p_x_n));
emptyBinThrValue = .2;
for pp = 1:size(p_x_n,1)
    emptyBin = isnan(sum(p_x_n{pp}));
    range = find(~emptyBin,1):find(~emptyBin,1,'last');        
    p_x_n{pp} = p_x_n{pp}(:,range);
    emptyBinPercentage = sum(isnan(sum(p_x_n{pp})))/length(range);
    if emptyBinPercentage >= emptyBinThrValue
        emptyBinThr(pp) = false;
    end
end
for st = 1:length(nType)
    outdir_figs = strcat(outdir,'\',nType{st});
    if ~isdir(outdir_figs)
        mkdir(outdir_figs)
    end
    %% Plot Bayesian examples
    nExample = 5;
    in = strcmp(sessionType,nType{st}) & emptyBinThr;
    % correct before trials
    pxn_correct = p_x_n(Correct_test==1 & in);
    stopLoc_correct = stopLoc(Correct_test==1 & in);
    rewardLoc_correct = rewardLoc(Correct_test==1 & in);
    r2_correct = r_squared(Correct_test==1 & in);
    % incorrect before trials
    pxn_incorrect = p_x_n(Correct_test==0 & in);
    stopLoc_incorrect = stopLoc(Correct_test==0 & in);
    rewardLoc_incorrect = rewardLoc(Correct_test==0 & in);
    r2_incorrect = r_squared(Correct_test==0 & in);
    h = figure;
    set(h,'OuterPosition',[2007,6,601,1075]);
    ha = zeros(nExample*2,1);
    clim = NaN(nExample*2,2);
    for ii = 1:nExample
        % correct examples
        plotInd = 2*ii-1;
        ha(plotInd) = subplot(nExample,2,plotInd);
        if ~isempty(pxn_correct)            
            pxn = pxn_correct{ii};
            taxis = 0:timeStep:size(pxn,2)*timeStep;
            imagesc(taxis,posaxis,pxn);
            hold on; plot(xlim,[stopLoc_correct(ii) stopLoc_correct(ii)],'b--')
            plot(xlim,[rewardLoc_correct(ii) rewardLoc_correct(ii)],'g--')
            axis xy square
            colormap hot
            clim(plotInd,:) = [nanmin(pxn(:)) nanmax(pxn(:))];    
            if ii == 1
                title(['Correct trials; r2 = ', num2str(r2_correct(ii))])
            else
                title(['r2 = ',num2str(r2_correct(ii))])
            end
        else
            axis off
        end
        if ii == ceil(nExample/2)
            ylabel('Angular position (rad)')
        end
        if ii == nExample
            xlabel('Time (sec)')
        end
        % incorrect examples
        plotInd = 2*ii;
        ha(plotInd) = subplot(nExample,2,plotInd);
        if ~isempty(pxn_incorrect) 
            pxn = pxn_incorrect{ii};
            taxis = 0:timeStep:size(pxn,2)*timeStep;
            imagesc(taxis,posaxis,pxn);
            hold on; plot(xlim,[stopLoc_incorrect(ii) stopLoc_incorrect(ii)],'b--')
            plot(xlim,[rewardLoc_incorrect(ii) rewardLoc_incorrect(ii)],'g--')
            axis xy square
            colormap hot
            clim(plotInd,:) = [nanmin(pxn(:)) nanmax(pxn(:))];
            if ii == 1
                title(['Incorrect trials; r2 = ', num2str(r2_incorrect(ii))])
            else
                title(['r2 = ', num2str(r2_incorrect(ii))])
            end
            if plotInd == length(ha)
                hb = colorbar;
                pos_bar = get(hb,'Position'); pos_bar(1) = pos_bar(1)+.05;
                ylabel(hb,'probability')
                set(hb,'Position',pos_bar);
                set(ha,'CLim',[nanmin(clim(:)) prctile(clim(:),60)])
            end
        end
        if ii == nExample
            xlabel('Time (sec)')
        end    
    end
    saveas(h,strcat(outdir_figs,'\BayesianDecodingExamples_',trial_ext{tt}),'epsc')
    saveas(h,strcat(outdir_figs,'\BayesianDecodingExamples_',trial_ext{tt}),'png')
    close(h)
    
    %% Plot r2 distributions between correct and incorrect trials
    if ~strcmp(nType{st},'PreRun')
        in = strcmp(sessionType,nType{st});
        r2min = min(r_squared);
        r2max = max(r_squared);
        sigma = .1;
        [Wdis_correct,axisValue1] = WeightedProportion(r_squared(Correct_test==1 & in),r2min,r2max,sigma);
        sample = bootstrp(5000,@WeightedProportion,r_squared(Correct_test==1 & in),r2min,r2max,sigma);
        [CI_u, CI_l] = CorrectionForMultipleComparsion(sample);
        CI_correct = [CI_u'-Wdis_correct;Wdis_correct-CI_l'];

        [Wdis_incorrect,axisValue2] = WeightedProportion(r_squared(Correct_test==0 & in),r2min,r2max,sigma);
        sample = bootstrp(5000,@WeightedProportion,r_squared(Correct_test==0 & in),r2min,r2max,sigma);
        [CI_u, CI_l] = CorrectionForMultipleComparsion(sample);
        CI_incorrect = [CI_u'-Wdis_incorrect; Wdis_incorrect-CI_l'];

        h = figure; hold on
        ha1 = shadedErrorBar(axisValue1,Wdis_correct,CI_correct,'r-');
        ha2 = shadedErrorBar(axisValue2,Wdis_incorrect,CI_incorrect,'k-');
        ylabel('Proportion')
        xlabel('r^2')
        legend([ha1.mainLine; ha2.mainLine],'String',...
            {strcat('Correct(n=',num2str(sum(Correct_test==1 & in)),')'); strcat('Incorrect(n=',num2str(sum(Correct_test==0 & in)),')')})
        saveas(h,strcat(outdir_figs,'\R2distCorrectVSIncorrect_',trial_ext{tt}),'epsc')
        saveas(h,strcat(outdir_figs,'\R2distCorrectVSIncorrect_',trial_ext{tt}),'png')
        close(h)
    end
end
end
%% Extract ripple data for plotting
RippleFR=[];
ActiveNearReward=[];
ActiveNearRest=[];
ratID = fieldnames(BayesianDecodingResult);
for rr = 1:length(ratID)
    dateID = fieldnames(BayesianDecodingResult.(ratID{rr}));
    for dd = 1:length(dateID)
        RippleFR = cat(1,RippleFR,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).RippleFR);
        ActiveNearReward = cat(1,ActiveNearReward,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).ActiveNearReward);
        ActiveNearRest = cat(1,ActiveNearRest,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).ActiveNearRest);
    end
end
RippleFR_xlabel = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).RippleFR_xlabel;
celltype = cell(size(ActiveNearReward));
celltype(logical(ActiveNearReward)&~logical(ActiveNearRest))={'Reward'};
celltype(logical(ActiveNearRest)&~logical(ActiveNearReward))={'Rest'};
celltype(~logical(ActiveNearRest)&~logical(ActiveNearReward))={'NonReward'};
celltype(logical(ActiveNearRest)&logical(ActiveNearReward))={'Reward&Rest'};

h = figure;
set(h,'OuterPosition',[12,397,2476,587])
ha = zeros(length(RippleFR_xlabel),1);
yvalues = [];
for ii = 1:length(RippleFR_xlabel)
    ha(ii) = subplot(1,length(RippleFR_xlabel),ii);
    boxplot(RippleFR(:,ii),celltype,'notch','on')
    yvalues = cat(2,yvalues,ylim);
    axis square
    title(RippleFR_xlabel{ii})
    if ii == 1
        ylabel('Ripple firing rate (Hz)')
    end
end
set(ha,'YLim',[min(yvalues) max(yvalues)])
saveas(h,strcat(outdir,'\RippleFR_allCondition'),'png')
saveas(h,strcat(outdir,'\RippleFR_allCondition'),'epsc')
close(h)

%% Helper functions
function [CurrentLoc,p_x_n,spkInd,r_squared,calphase,slope,sessionType,NUniqueCell,NTrial,NTrial_post,NTrial_preRun,...
    RippleOnsetIndex,RippleOffsetIndex,beforeReward,RippleFR,ActiveNearReward,ActiveNearRest,RippleFR_xlabel,smspkFR,cellID,param] = ...
    BayesianDecoding_v2(spklist,encodingfd,posIndfile,eeg2posIndfile,spikeIndfile,trackdata,RippleOnsetIndex,RippleOffsetIndex,eegt)

%parameters
param.posbinsize = 4; % in degree
param.baysianTimeWindow = .02; %in sec
param.step = .01; % overlap 10 msec
param.FRkernel = gausskernel(5,2);
param.FRoffset = 0.0001; % in Hz
param.run_thr = 5; % in cm/s

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
[~,spkFR] = GenerateSpkMatrix_v6(eeg2posIndfile,posIndfile,param.posbinsize,false,spikeIndfile,encodingfd);

%% Obtain spike raster from all trials
[SpikeMatrix,~,posaxis] = GenerateSpkMatrix_v6(eeg2posIndfile,posIndfile,param.posbinsize,false,spikeIndfile);
posaxis = circ_ang2rad(posaxis);

%Remove cells with peak firing rate less than 1 Hz
rmv = max(spkFR,[],2)<1;
spkFR = spkFR(~rmv,:);
cellID = cellID_noExt(~rmv);

smspkFR = zeros(size(spkFR));
for ss = 1:size(smspkFR,1)
    smspkFR(ss,:) = conv_cir(spkFR(ss,:),param.FRkernel);
end
smspkFR = smspkFR+param.FRoffset; %add offset to prevent 0 value

%% Select cells with place fields near reward position
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
ind = match(rewardPos_range,posaxis);
ActiveNearReward = max(smspkFR(:,ind(1):ind(2)),[],2) >= 1;

%% Select cells with place fields near rest box position
restPos_range = [trackdata.Ang_StartZone_depart_ontrack trackdata.Ang_StartZone_arrive_ontrack];
ind = match(restPos_range,posaxis);
logical_ind = false(size(posaxis));
logical_ind(1:ind(1)) = true; logical_ind(ind(2):end) = true; % rest box located near 0 on circular track
ActiveNearRest= max(smspkFR(:,logical_ind),[],2) >= 1;

%% Extract time index outside of rest box
[Test,Sample,PostTest,PreRun]= ExtractTimeInd_OutOfRestBox(eeg2posIndfile,posIndfile);

%% Decode position from ripple events
load(eeg2posIndfile);
load(posIndfile);
pos = data_angle_all(eeg2pos,2);
CurrentLoc = NaN(length(RippleOnsetIndex),1);
p_x_n = cell(length(RippleOnsetIndex),1);
spkInd = cell(length(RippleOnsetIndex),1);
r_squared = zeros(length(RippleOnsetIndex),1);
calphase = cell(length(RippleOnsetIndex),1);
slope = zeros(length(RippleOnsetIndex),1);
sessionType = cell(length(RippleOnsetIndex),1);
NUniqueCell = zeros(length(RippleOnsetIndex),1);
NTrial = zeros(length(RippleOnsetIndex),2);
NTrial_post = zeros(length(RippleOnsetIndex),1);
NTrial_preRun = zeros(length(RippleOnsetIndex),1);
beforeReward = false(length(RippleOnsetIndex),1);
dur = (RippleOffsetIndex-RippleOnsetIndex)/sampFreq;
velocity = data_angle_all(eeg2pos,3);
ripple_in = false(size(RippleOnsetIndex));
ripple_time = zeros(sum(~rmv),4);
ripple_spkcount = zeros(sum(~rmv),4);
RippleFR_xlabel = cell(1,4);
for ii = 1:length(RippleOnsetIndex)
    InRange = RippleOnsetIndex(ii):RippleOffsetIndex(ii);
    if mean(velocity(InRange)) > param.run_thr % exclude ripple events while rats move
        continue
    end
    CurrentLoc(ii) = pos(RippleOnsetIndex(ii));
    timeaxis = 0:param.step:dur(ii);
    if sum(~isnan(Sample(InRange))) == length(InRange) || sum(~isnan(Test(InRange))) == length(InRange)...
        || sum(~isnan(PreRun(InRange))) == length(InRange) || sum(~isnan(PostTest(InRange))) == length(InRange)
        if full(sum(sum(SpikeMatrix(~rmv,InRange)))) > 0
            ripple_in(ii) = true;
            p_x_n{ii,1} = BayesianDecoder(full(SpikeMatrix(~rmv,InRange)),smspkFR,param.baysianTimeWindow,param.step,sampFreq);
            [cellind,tind] = find(full(SpikeMatrix(~rmv,InRange)));
            spkInd{ii} = [cellind, tind];
            pxn = p_x_n{ii,1};
            bin2use = find(~isnan(sum(pxn,1)));
            [r_squared(ii,1),calphase{ii,1},~,~,slope(ii,1)] = Cir_reg(pxn,posaxis',timeaxis,bin2use);
            NUniqueCell(ii) = sum(sum(full(SpikeMatrix(~rmv,InRange)),2)>0);
            spkcount = full(sum(SpikeMatrix(~rmv,InRange),2));
            if sum(~isnan(Test(InRange))) == length(InRange) || sum(~isnan(Sample(InRange))) == length(InRange)
                if mode(Test(InRange))<=4 || mode(Sample(InRange))<=4
                   sessionType{ii} = 'InitialFour';
                   if mode(Test(InRange))<=4
                        NTrial(ii,2) = mode(Test(InRange));
                   elseif mode(Sample(InRange))<=4
                        NTrial(ii,1) = mode(Sample(InRange));
                   end
                   ripple_time(:,1) = ripple_time(:,1) + length(InRange)/sampFreq;
                   ripple_spkcount(:,1) = ripple_spkcount(:,1) + spkcount;
                   RippleFR_xlabel{1} = sessionType{ii};
                elseif mode(Test(InRange))>4 || mode(Sample(InRange))>4
                   sessionType{ii} = 'LastFour';
                   if mode(Test(InRange))>4
                        NTrial(ii,2) = mode(Test(InRange));
                   elseif mode(Sample(InRange))>4
                        NTrial(ii,1) = mode(Sample(InRange));
                   end
                   ripple_time(:,2) = ripple_time(:,2) + length(InRange)/sampFreq;
                   ripple_spkcount(:,2) = ripple_spkcount(:,2) + spkcount;
                   RippleFR_xlabel{2} = sessionType{ii};
                end
                if sum(~isnan(Test(InRange))) == length(InRange)
                    ntrial = mode(Test(InRange));
                    onset = trackdata.ts_test_start(ntrial)/1e6;
                    offset = trackdata.ts_test_reward(ntrial)/1e6;
                    if eegt(RippleOnsetIndex(ii)) >= onset/1e6 && eegt(RippleOffsetIndex(ii)) <= offset/1e6
                        beforeReward(ii) = true;
                    end
                elseif sum(~isnan(Sample(InRange))) == length(InRange)
                    ntrial = mode(Sample(InRange));
                    onset = trackdata.ts_sample_start(ntrial)/1e6;
                    offset = trackdata.ts_sample_reward(ntrial)/1e6;
                    if eegt(RippleOnsetIndex(ii)) >= onset && eegt(RippleOffsetIndex(ii)) <= offset
                        beforeReward(ii) = true;
                    end                    
                end
            elseif sum(~isnan(PostTest(InRange))) == length(InRange)
               sessionType{ii} = 'PostTest';
               NTrial_post(ii) = mode(PostTest(InRange));
               ripple_time(:,3) = ripple_time(:,3) + length(InRange)/sampFreq;
               ripple_spkcount(:,3) = ripple_spkcount(:,3) + spkcount;
               RippleFR_xlabel{3} = sessionType{ii};
               
               onset = trackdata.ts_posttest_start(NTrial_post(ii))/1e6;
               offset = trackdata.ts_posttest_reward(NTrial_post(ii))/1e6;
               if eegt(RippleOnsetIndex(ii)) >= onset && eegt(RippleOffsetIndex(ii)) <= offset
                    beforeReward(ii) = true;
               end   
           elseif sum(~isnan(PreRun(InRange))) == length(InRange)
               sessionType{ii} = 'PreRun';
               NTrial_preRun(ii) = mode(PreRun(InRange));
               ripple_time(:,4) = ripple_time(:,4) + length(InRange)/sampFreq;
               ripple_spkcount(:,4) = ripple_spkcount(:,4) + spkcount;
               RippleFR_xlabel{4} = sessionType{ii};
            end
        end
    end
end

CurrentLoc = CurrentLoc(ripple_in);
p_x_n = p_x_n(ripple_in);
spkInd = spkInd(ripple_in);
r_squared = r_squared(ripple_in);
calphase = calphase(ripple_in);
slope = slope(ripple_in);
sessionType = sessionType(ripple_in);
NUniqueCell = NUniqueCell(ripple_in);
NTrial = NTrial(ripple_in,:);
NTrial_post = NTrial_post(ripple_in);
NTrial_preRun = NTrial_preRun(ripple_in);
RippleOnsetIndex = RippleOnsetIndex(ripple_in);
RippleOffsetIndex = RippleOffsetIndex(ripple_in);
beforeReward = beforeReward(ripple_in);
RippleFR = ripple_spkcount./ripple_time;
param.posbinsize = param.posbinsize/360*2*pi; % convert to radian

function [Test,Sample,PostTest,PreRun]= ExtractTimeInd_OutOfRestBox(eeg2posIndfile,posIndfile)
load(eeg2posIndfile);
load(posIndfile);

nTestTrial = length(data_angle{3});
Sample = NaN(size(eeg2pos));
Test = NaN(size(eeg2pos));
for t = 1:nTestTrial
    % sample trial
    SampleTrial_onset = data_angle{2}{t}(1,1);
    SampleTrial_offset = data_angle{2}{t}(end,1);
    % position index
    ind_range = [find((data_angle_all(:,1)-SampleTrial_onset)== 0),...
                 find((data_angle_all(:,1)-SampleTrial_offset)== 0)];
    % translate into eeg index
    eegind_range = [find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1),...
                    find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1,'last')];
    Sample(min(eegind_range):max(eegind_range)) = t;
    
    % test trial
    TestTrial_onset = data_angle{3}{t}(1,1);
    TestTrial_offset = data_angle{3}{t}(end,1);
     % position index
    ind_range = [find((data_angle_all(:,1)-TestTrial_onset)== 0),...
                 find((data_angle_all(:,1)-TestTrial_offset)== 0)];
    % translate into eeg index
    eegind_range = [find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1),...
                    find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1,'last')];
    Test(min(eegind_range):max(eegind_range)) = t;
end

nPrerunningTrial = length(data_angle{1});
PreRun = NaN(size(eeg2pos));
for t = 1:nPrerunningTrial
    % pre-running trial
    preRun_onset = data_angle{1}{t}(1,1);
    preRun_offset = data_angle{1}{t}(end,1);
    % position index
    ind_range = [find((data_angle_all(:,1)-preRun_onset)== 0),...
                 find((data_angle_all(:,1)-preRun_offset)== 0)];
    % translate into eeg index
    eegind_range = [find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1),...
                    find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1,'last')];
    PreRun(min(eegind_range):max(eegind_range)) = t;
end

nPostTestTrial = length(data_angle{4});
PostTest = NaN(size(eeg2pos));
for t = 1:nPostTestTrial
    % in between post test run in rest box
    PostTestTrial_onset = data_angle{4}{t}(1,1);
    PostTestTrial_offset = data_angle{4}{t}(end,1);
    % position index
    ind_range = [find((data_angle_all(:,1)-PostTestTrial_onset)== 0),...
                 find((data_angle_all(:,1)-PostTestTrial_offset)== 0)];
    % translate into eeg index
    eegind_range = [find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1),...
                    find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1,'last')];
    PostTest(min(eegind_range):max(eegind_range)) = t;
end

checkoverlap = ~isnan([Sample Test PreRun PostTest]);
checkoverlap = sum(checkoverlap,2);
if any(checkoverlap > 1)
    error('Event index overlaps')
end

function [CI_u, CI_l] = CorrectionForMultipleComparsion(sample)
%Simultaneous bounds are wider than separate bounds, because it is more stringent to require that the entire 
%curve be within the bounds than to require that the curve at a single predictor value be within the bounds.
    [sample_rank] = tiedrank(sample);
    sample_rank=sample_rank/size(sample_rank,1);

    sample_min=prctile(min(sample_rank,[],2),2.5);
    sample_max=prctile(max(sample_rank,[],2),97.5);

    CI_u=zeros(size(sample,2),1);
    for ii = 1:size(sample,2);
        if sum(sample_rank(:,ii)>=sample_max) == 0
            CI_u(ii) = max(sample(:,ii));
        else
            CI_u(ii) = min(sample(sample_rank(:,ii)>=sample_max,ii));
        end
    end

    CI_l=zeros(size(sample,2),1);
    for ii = 1:size(sample,2);
        if sum(sample_rank(:,ii)<=sample_min) == 0
            CI_l(ii) = min(sample(:,ii));
        else
            CI_l(ii) = max(sample(sample_rank(:,ii)<=sample_min,ii));
        end
    end
    
function [wDist,xaxis] = WeightedProportion(x,xmin,xmax,sigma)
    xaxis = xmin:sigma/5:xmax;
    delta = repmat(x,1,length(xaxis))-repmat(xaxis,length(x),1);
    W = exp(-0.5*delta.*delta/sigma^2);
    W = W./repmat(sum(W,2),1,size(W,2)); %normalize each kernel so the ones on the edge have same contribution
    wDist = sum(W)./sum(sum(W));