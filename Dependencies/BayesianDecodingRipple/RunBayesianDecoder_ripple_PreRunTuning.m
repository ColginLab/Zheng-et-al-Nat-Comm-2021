function RunBayesianDecoder_ripple_PreRunTuning(pathRats,CSClist,Outdir,cellList,spktime_fd,trackdata,plotOnly)

outdir = strcat(Outdir,'\BayesianDecoding_ripple_PreRunTuning');
if ~isdir(outdir)
    mkdir(outdir)
end
if exist(strcat(outdir,'\BayesianDecodingResult.mat'),'file')==2 && plotOnly
    load(strcat(outdir,'\BayesianDecodingResult.mat'))
else
    for ii = 1:length(pathRats)

        % Extract ratID and dateID
        ind = strfind(pathRats{ii},'Rat');
        RatID = pathRats{ii}(ind:ind+5);
        ind = strfind(pathRats{ii},'\');
        DateID = pathRats{ii}(ind(end-1)+1:end-1);
        ind = strfind(DateID,'-');
        DateID(ind)=[];
        DateID = strcat('D',DateID);

        spikeIndfile = strcat(pathRats{ii},'\',spktime_fd,'\spike_index.mat');
        eeg2posIndfile = strcat(pathRats{ii},'\',spktime_fd,'\eeg2pos.mat');
        posIndfile = strcat(pathRats{ii},'\Data_angle_ontrack.mat');
        behavioralData = load(trackdata{ii});
        spkfile = strcat(pathRats{ii},'\',cellList);

        sampfreq = 2000;
        n_std = 5; %number of standard deviations for the threshold
        duration_threshold = .05; %in sec. duration of ripple event
        CA1eeg = {};
        for cc = 1:length(CSClist{ii})
            CSClist_CA1 = strcat(pathRats{ii},'\CSC',num2str(CSClist{ii}(cc)),'.ncs');
            [eeg,~,sampfreq] = LoadEEG(CSClist_CA1);
            CA1eeg = cat(1,CA1eeg,eeg);
        end
        if isempty(CA1eeg)
            continue
        end
        [RippleOnsetIndex, RippleOffsetIndex] = DetectRipples_v4(CA1eeg,sampfreq,n_std,duration_threshold);

        % Define encoding folders
        %encodingfd = [2 3 4];
        encodingfd = [1];
        [CurrentLoc,p_x_n,spkind,r_squared,calphase,slope,sessionType,NUniqueCell,NTrial,NTrial_post,NTrial_preRun,...
            RippleOnsetIndex,RippleOffsetIndex,RippleFR,ActiveNearReward,ActiveNearRest,RippleFR_xlabel,smspkFR,cellID,param] = ...    
        BayesianDecoding(spkfile,encodingfd,posIndfile,eeg2posIndfile,spikeIndfile,behavioralData,RippleOnsetIndex,RippleOffsetIndex);

        BayesianDecodingResult.(RatID).(DateID).CurrentLoc = CurrentLoc;
        BayesianDecodingResult.(RatID).(DateID).p_x_n = p_x_n;
        BayesianDecodingResult.(RatID).(DateID).spkind = spkind;
        BayesianDecodingResult.(RatID).(DateID).r_squared = r_squared;
        BayesianDecodingResult.(RatID).(DateID).calphase = calphase;
        BayesianDecodingResult.(RatID).(DateID).slope = slope;
        BayesianDecodingResult.(RatID).(DateID).sessionType = sessionType;
        BayesianDecodingResult.(RatID).(DateID).NUniqueCell = NUniqueCell;
        BayesianDecodingResult.(RatID).(DateID).NTrial = NTrial;
        BayesianDecodingResult.(RatID).(DateID).NTrial_post = NTrial_post;
        BayesianDecodingResult.(RatID).(DateID).NTrial_preRun = NTrial_preRun;
        BayesianDecodingResult.(RatID).(DateID).RippleOnsetIndex = RippleOnsetIndex;
        BayesianDecodingResult.(RatID).(DateID).RippleOffsetIndex = RippleOffsetIndex;
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
%% Extract relevant data for plotting
p_x_n = {};
r_squared = [];
slope = [];
sessionType = {};
Correct_test = [];
stopLoc = [];
rewardLoc = [];
currentLoc = [];
NUniqueCell_thr = 8;
NTrial_all = [];
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
            %include1 = ismember(NTrial,correct_sample); % include only correct sample trial
            include1 = NTrial~=0; % 0 corresponds to preRun and postRun trials 
            include2 = ~isnan(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared) & NUniqueCell > NUniqueCell_thr;
            correct_test = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sign_correct_test;
            in = correct_test==1;
            correct_test = find(in);
            
            NTrial_all = cat(1,NTrial_all,NTrial(include1 & include2));
            p_x_n = cat(1,p_x_n,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).p_x_n(include1 & include2));
            r_squared = cat(1,r_squared,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared(include1 & include2));
            slope = cat(1,slope,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).slope(include1 & include2));
            sessionType = cat(1,sessionType,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sessionType(include1 & include2));
            Correct_test = cat(1,Correct_test,ismember(NTrial(include1 & include2),correct_test));
            stopLoc = cat(1,stopLoc,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).ang_test_reward_ontrack(NTrial(include1 & include2)));
            rewardLoc = cat(1,rewardLoc,repmat(rewardLoc_session,sum(include1 & include2),1));
            currentLoc = cat(1,currentLoc,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).CurrentLoc(include1 & include2));
        end

        if ~isempty(NTrial_post)
            include3 = NTrial_post~=0 & ~isnan(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared) & NUniqueCell > NUniqueCell_thr;
            correct_posttest = BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sign_correct_posttest;
            in = correct_posttest==1;
            correct_posttest = find(in);
            
            NTrial_all = cat(1,NTrial_all,NaN(sum(include3),1));
            p_x_n = cat(1,p_x_n,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).p_x_n(include3));
            r_squared = cat(1,r_squared,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared(include3));
            slope = cat(1,slope,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).slope(include3));
            sessionType = cat(1,sessionType,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sessionType(include3));
            Correct_test = cat(1,Correct_test,ismember(NTrial_post(include3),correct_posttest));
            stopLoc = cat(1,stopLoc,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).ang_posttest_reward_ontrack(NTrial_post(include3)));
            rewardLoc = cat(1,rewardLoc,repmat(rewardLoc_session,sum(include3),1));
            currentLoc = cat(1,currentLoc,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).CurrentLoc(include3));
        end

        if ~isempty(NTrial_preRun)
            include4 = NTrial_preRun~=0 & ~isnan(BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared) & NUniqueCell > NUniqueCell_thr;
            
            NTrial_all = cat(1,NTrial_all,NaN(sum(include4),1));
            p_x_n = cat(1,p_x_n,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).p_x_n(include4));
            r_squared = cat(1,r_squared,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).r_squared(include4));
            slope = cat(1,slope,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).slope(include4));
            sessionType = cat(1,sessionType,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).sessionType(include4));
            Correct_test = cat(1,Correct_test,false(sum(include4),1));
            stopLoc = cat(1,stopLoc,NaN(sum(include4),1));
            rewardLoc = cat(1,rewardLoc,repmat(rewardLoc_session,sum(include4),1));
            currentLoc = cat(1,currentLoc,BayesianDecodingResult.(ratID{rr}).(dateID{dd}).CurrentLoc(include4));
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
Correct_test = logical(Correct_test);

%% Plot average posterior probability aligned by stop or reward locations
nType = unique(sessionType);
% Type = {'Correct'; 'Incorrect stop'; 'Incorrect reward'};
% ha = zeros(length(Type),length(nType)-1);
% ha_cum = zeros(length(Type),length(nType)-1);
% yvalue = [];
% h = figure;
% h_cum = figure;
% set([h_cum h],'OuterPosition',[2018,3,1129,1078])
% ProbabilityThr = 0.3;
% for st = 1:length(nType)
%     in = strcmp(sessionType,nType{st});
%     if strcmp(nType{st},'PreRun')       
%         pxn_subset = p_x_n(in);
%         rewardLoc_subset = rewardLoc(in);
%         shift = match(rewardLoc_subset,posaxis)-length(posaxis)/2;
%         Ntimebin = cellfun(@size,pxn_subset,repmat({2},size(pxn_subset)));
%         pxn_aligned = zeros(length(posaxis),sum(Ntimebin));
%         n = 1;
%         for pp = 1:length(pxn_subset)
%             n2 = n+Ntimebin(pp);
%             pxn_aligned(:,n:n2-1) = circshift(pxn_subset{pp},[-shift(pp),0]);
%             n = n2;
%         end
%         
%         % Obtain cumulative proportion of time bins passing probability
%         % threshold
%         pxn_aligned_rmv = pxn_aligned(:,~isnan(sum(pxn_aligned)));
%         DistanceFromReward = 0:posbinsize:pi;
%         ProportionOfTimeBin = zeros(1,length(DistanceFromReward));
%         for ii = 1:length(DistanceFromReward)
%             in2 = abs(posaxis-pi)<=DistanceFromReward(ii);
%             CumProbability = sum(pxn_aligned_rmv(in2,:),1)';
%             ProportionOfTimeBin(ii) = sum(CumProbability>=ProbabilityThr,1)/size(pxn_aligned_rmv,2);
%         end
%         
%         pxn_aligned = nanmedian(pxn_aligned,2);        
%         
%         %Permutation test
%         nShuffle = 1000;
%         pxn_aligned_shu = zeros(length(posaxis),nShuffle);
%         PossibleShift = unique(shift);
%         ind = zeros(length(shift),1);
%         for ii = 1:length(PossibleShift)
%             ind(shift==PossibleShift(ii),1) = ii;
%         end
%         ProportionOfTimeBin_shu = zeros(nShuffle,length(DistanceFromReward));
%         for nn = 1:nShuffle
%             pxn_shu = zeros(length(posaxis),sum(Ntimebin));
%             n = 1;
%             %shuffle = shift(randperm(length(shift),length(shift)));
%             %shuffle = circshift(shift,randi(length(shift)));
%             shu = PossibleShift(randperm(length(PossibleShift),length(PossibleShift))); shuffle = shu(ind);
%             for pp = 1:length(pxn_subset)
%                 n2 = n+Ntimebin(pp);
%                 pxn_shu(:,n:n2-1) = circshift(pxn_subset{pp},[-shuffle(pp),0]);
%                 n = n2;
%             end 
%             
%             % Obtain cumulative proportion of time bins passing probability
%             % threshold
%             pxn_shu_rmv = pxn_shu(:,~isnan(sum(pxn_shu)));                
%             for ii = 1:length(DistanceFromReward)
%                 in2 = abs(posaxis-pi)<=DistanceFromReward(ii);
%                 CumProbability = sum(pxn_shu_rmv(in2,:),1)';
%                 ProportionOfTimeBin_shu(nn,ii) = sum(CumProbability>=ProbabilityThr,1)/size(pxn_shu_rmv,2);
%             end
%             pxn_aligned_shu(:,nn) = nanmedian(pxn_shu,2);
%         end
%           CI_u = prctile(pxn_aligned_shu,97.5,2);
%           CI_l = prctile(pxn_aligned_shu,2.5,2);
%         %[CI_u, CI_l] = CorrectionForMultipleComparsion(pxn_aligned_shu');
%         
%         h2 = figure;
%         set(h2,'OuterPosition',[2018,3,1129,1078])
%         subplot(length(Type),length(nType)-1,1);
%         hold on;
%         null_median = nanmedian(pxn_aligned_shu,2);
%         error = [CI_u-null_median null_median-CI_l];
%         temp = shadedErrorBar((posaxis-pi)/pi*180,null_median,error,'k');
%         set(temp.mainLine,'LineWidth',1)
%         ha1 = plot((posaxis-pi)/pi*180,pxn_aligned,'r');
%         set(ha1,'LineWidth',1)
%         axis tight square
%         set(gca,'YScale','log','YLim',[min(yvalue) max(yvalue)])
%         ylabel('Median probability')
%         xlabel('Aligned angular position (deg)')
%         legend([ha1 temp.mainLine],'Observed','Shuffle')
%         str = sprintf('%s (n=%d)',nType{st},length(pxn_subset));
%         title(str,'Interpreter','none')
%         
%         % plot cumulative proportion of time bins
%         h2_cum = figure;
%         set(h2_cum,'OuterPosition',[2018,3,1129,1078])
%         subplot(length(Type),length(nType)-1,1);        
%         hold on;
%         CI_u_cum = prctile(ProportionOfTimeBin_shu,97.5);
%         CI_l_cum = prctile(ProportionOfTimeBin_shu,2.5);
%         error = [CI_u_cum-mean(ProportionOfTimeBin_shu); mean(ProportionOfTimeBin_shu)-CI_l_cum];
%         temp = shadedErrorBar(DistanceFromReward/pi*180,mean(ProportionOfTimeBin_shu),error,'k');
%         ha1 = plot(DistanceFromReward/pi*180,ProportionOfTimeBin,'r');
%         axis square tight
%         set([temp.mainLine ha1],'LineWidth',1)
%         ylabel('Proportion of time bins')
%         xlabel('Angular distance from reward (deg)')
% 
%         saveas(h2,strcat(outdir,'\BayesianDecodingRipples_preRun'),'epsc')
%         saveas(h2,strcat(outdir,'\BayesianDecodingRipples_preRun'),'png')
%         saveas(h2_cum,strcat(outdir,'\BayesianDecodingRipples_preRun_cum'),'epsc')
%         saveas(h2_cum,strcat(outdir,'\BayesianDecodingRipples_preRun_cum'),'png')
%         close(h2,h2_cum)
%     else
%         for tt = 1:length(Type)
%             if ~isempty(strfind(Type{tt},'Correct'))
%                 pxn_subset = p_x_n(Correct_test & in);
%                 stopLoc_subset = stopLoc(Correct_test & in);
%             elseif ~isempty(strfind(Type{tt},'Incorrect'))
%                 pxn_subset = p_x_n(~Correct_test & in);
%                 stopLoc_subset = stopLoc(~Correct_test & in);
%                 rewardLoc_subset = rewardLoc(~Correct_test & in);
%             end 
% 
%             if ~isempty(strfind(Type{tt},'Incorrect reward'))
%                 shift = match(rewardLoc_subset,posaxis)-length(posaxis)/2;
%             else
%                 shift = match(stopLoc_subset,posaxis)-length(posaxis)/2;
%             end
% 
%             Ntimebin = cellfun(@size,pxn_subset,repmat({2},size(pxn_subset)));
%             pxn_aligned = zeros(length(posaxis),sum(Ntimebin));
%             n = 1;
%             for pp = 1:length(pxn_subset)
%                 n2 = n+Ntimebin(pp);
%                 pxn_aligned(:,n:n2-1) = circshift(pxn_subset{pp},[-shift(pp),0]);
%                 n = n2;
%             end
%             
%             % Obtain cumulative proportion of time bins passing probability
%             % threshold
%             pxn_aligned_rmv = pxn_aligned(:,~isnan(sum(pxn_aligned)));
%             DistanceFromReward = 0:posbinsize:pi;
%             ProportionOfTimeBin = zeros(1,length(DistanceFromReward));
%             for ii = 1:length(DistanceFromReward)
%                 in2 = abs(posaxis-pi)<=DistanceFromReward(ii);
%                 CumProbability = sum(pxn_aligned_rmv(in2,:),1)';
%                 ProportionOfTimeBin(ii) = sum(CumProbability>=ProbabilityThr,1)/size(pxn_aligned_rmv,2);
%             end
%                 
%             pxn_aligned = nanmedian(pxn_aligned,2);
% 
%             %Permutation test
%             nShuffle = 1000;
%             pxn_aligned_shu = zeros(length(posaxis),nShuffle);
%             PossibleShift = unique(shift);
%             ind = zeros(length(shift),1);
%             for ii = 1:length(PossibleShift)
%                 ind(shift==PossibleShift(ii),1) = ii;
%             end
%             ProportionOfTimeBin_shu = zeros(nShuffle,length(DistanceFromReward));
%             for nn = 1:nShuffle
%                 pxn_shu = zeros(length(posaxis),sum(Ntimebin));
%                 n = 1;
%                 %shuffle = shift(randperm(length(shift),length(shift)));
%                 %shuffle = circshift(shift,randi(length(shift)));
%                 shu = PossibleShift(randperm(length(PossibleShift),length(PossibleShift))); shuffle = shu(ind);
%                 for pp = 1:length(pxn_subset)
%                     n2 = n+Ntimebin(pp);
%                     pxn_shu(:,n:n2-1) = circshift(pxn_subset{pp},[-shuffle(pp),0]);
%                     n = n2;
%                 end
%                 
%                 % Obtain cumulative proportion of time bins passing probability
%                 % threshold
%                 pxn_shu_rmv = pxn_shu(:,~isnan(sum(pxn_shu)));                
%                 for ii = 1:length(DistanceFromReward)
%                     in2 = abs(posaxis-pi)<=DistanceFromReward(ii);
%                     CumProbability = sum(pxn_shu_rmv(in2,:),1)';
%                     ProportionOfTimeBin_shu(nn,ii) = sum(CumProbability>=ProbabilityThr,1)/size(pxn_shu_rmv,2);
%                 end
%                 pxn_aligned_shu(:,nn) = nanmedian(pxn_shu,2);
%             end
%               CI_u = prctile(pxn_aligned_shu,97.5,2);
%               CI_l = prctile(pxn_aligned_shu,2.5,2);
%             %[CI_u, CI_l] = CorrectionForMultipleComparsion(pxn_aligned_shu');
%             
%             figure(h);
%             ha(tt,st) = subplot(length(Type),length(nType)-1,st+(tt-1)*(length(nType)-1));
%             hold on;
%             null_median = nanmedian(pxn_aligned_shu,2);
%             error = [CI_u-null_median null_median-CI_l];
%             temp = shadedErrorBar((posaxis-pi)/pi*180,null_median,error,'k');
%             set(temp.mainLine,'LineWidth',1)
%             ha1 = plot((posaxis-pi)/pi*180,pxn_aligned,'r');
%             set(ha1,'LineWidth',1)
%             axis tight square
%             set(gca,'YScale','log')
%             yvalue = cat(2,yvalue,ylim);
%             if st == length(nType)-1 && tt == length(Type)            
%                 legend([ha1 temp.mainLine],'Observed','Shuffle')
%             end
%             if st == 1
%                 ylabel('Median probability')
%             end
%             if tt == length(Type)
%                 xlabel('Aligned angular position (deg)')
%             end
%             if tt == 1
%                 str = sprintf('%s\n%s (n=%d)',nType{st},Type{tt},length(pxn_subset));
%                 title(str,'Interpreter','none')
%             else
%                 str = sprintf('%s (n=%d)',Type{tt},length(pxn_subset));
%                 title(str,'Interpreter','none')
%             end
%             
%             % plot cumulative proportion of time bins
%             figure(h_cum)
%             ha_cum(tt,st) = subplot(length(Type),length(nType)-1,st+(tt-1)*(length(nType)-1));
%             hold on;
%             CI_u_cum = prctile(ProportionOfTimeBin_shu,97.5);
%             CI_l_cum = prctile(ProportionOfTimeBin_shu,2.5);
%             error = [CI_u_cum-mean(ProportionOfTimeBin_shu); mean(ProportionOfTimeBin_shu)-CI_l_cum];
%             temp = shadedErrorBar(DistanceFromReward/pi*180,mean(ProportionOfTimeBin_shu),error,'k');
%             ha1 = plot(DistanceFromReward/pi*180,ProportionOfTimeBin,'r');
%             axis square tight
%             set([temp.mainLine ha1],'LineWidth',1)
%             if st == 1
%                 ylabel('Proportion of time bins')
%             end
%             if tt == length(Type)
%                 xlabel('Angular distance from reward (deg)')
%             end
%         end
%     end
% end
% set(ha,'YLim',[min(yvalue) max(yvalue)]);
% set(ha_cum,'YLim',[0 1]);
% saveas(h,strcat(outdir,'\BayesianDecodingRipples'),'epsc')
% saveas(h,strcat(outdir,'\BayesianDecodingRipples'),'png')
% saveas(h_cum,strcat(outdir,'\BayesianDecodingRipples_cum'),'epsc')
% saveas(h_cum,strcat(outdir,'\BayesianDecodingRipples_cum'),'png')
% close(h,h_cum)

r2_thr = 0.6;
%% Examine whether there is bias from starting position at the beginning of replay
% align current position to the middle
ind_half = match(pi,posaxis);
ind = match(currentLoc,posaxis);
shift = ind-ind_half;
pxn_aligned = p_x_n;
for ii = 1:length(pxn_aligned)
    pxn_aligned{ii} = circshift(pxn_aligned{ii},-shift(ii));
end
% combine replay events in each type of session
sessions = unique(sessionType);
sessions = sessions([4;1;2;3]);

% Construct bar graph for possible goal location
LocDis = 16; % in deg
nLocAway = 9;
locationRangeInd = match(-LocDis*nLocAway:LocDis:LocDis*nLocAway,(posaxis-pi)*180/pi);
ProbSum = zeros(length(locationRangeInd),length(sessions));
pvalue = zeros(1,length(sessions));

h = figure;
set(h,'OuterPosition',[13.8,309.8,1511.2,551.2])
ha = zeros(length(sessions),1);
for ss = 1:length(sessions)
    in = strcmp(sessionType,sessions{ss}) & r_squared >= r2_thr;
    pxn = pxn_aligned(in);
    [~,tbinL] = cellfun(@size,pxn);
    pxn_combined = NaN(length(posaxis),min(tbinL),length(pxn));
    for pp = 1:length(pxn)
        % align onset of replay
        range = 1:min(tbinL);
        pxn_combined(:,:,pp) = pxn{pp}(:,range);
    end
    Out.(sessions{ss}) = pxn_combined;
    
    ha(ss) = subplot(2,length(sessions),ss);
    taxis = 0:timeStep:(min(tbinL)-1)*timeStep;
    imagesc(taxis,(posaxis-pi)*180/pi,nanmean(pxn_combined,3))
    axis xy square
    hold on
    plot(xlim,[0 0],'w--','LineWidth',1)
    cb = colorbar;
    title([sessions{ss},' n=',num2str(sum(in))])
    if ss == 1
        xlabel('Time relative to replay onset (sec)')
        ylabel('Position relative to current position (deg)')
    elseif ss == length(sessions)
        ylabel(cb,'Mean posterior probability')
    end
    
    % Add up probability
    for pp = 1:length(locationRangeInd)
        inrange = locationRangeInd(pp)-1:locationRangeInd(pp)+1;
        ProbSum(pp,ss) = nansum(nansum(nansum(pxn_combined(inrange,:,:))));
    end
    [~,pvalue(ss)] = lillietest(ProbSum(:,ss));
    subplot(2,length(sessions),ss+length(sessions))
    bar(-nLocAway:nLocAway,ProbSum(:,ss))
    axis square
    colorbar
    if ss == 1
        xlabel('Location # relative to current position')
        ylabel('sum of posterior probability')
    end
end
colormap hot
set(ha,'CLim',[0 .04])
saveas(h,strcat(outdir,'\InitiationBias'),'epsc')
saveas(h,strcat(outdir,'\InitiationBias'),'png')
close(h)

% test for outlier
tf = isoutlier(ProbSum,'gesd');
t = array2table(tf);
t2 = array2table(pvalue);
t.Properties.VariableNames=sessions;
t2.Properties.VariableNames=sessions;
Stats.ReplayOnset.OutlierTest = t;
Stats.ReplayOnset.Lillietest = t2;

%% Examine whether there is bias towards goal at the end of replay
% align reward position to the middle
ind = match(rewardLoc,posaxis);
ind_half = match(pi,posaxis);
shift = ind-ind_half;
pxn_aligned = p_x_n;
for ii = 1:length(pxn_aligned)
    pxn_aligned{ii} = circshift(pxn_aligned{ii},-shift(ii));
end
% Construct bar graph for possible goal location
LocDis = 16; % in deg
nLocAway = 9;
locationRangeInd = match(-LocDis*nLocAway:LocDis:LocDis*nLocAway,(posaxis-pi)*180/pi);
ProbSum = zeros(length(locationRangeInd),length(sessions));
pvalue = zeros(1,length(sessions));

% combine replay events in each type of session
h = figure;
set(h,'OuterPosition',[13.8,309.8,1511.2,551.2])
ha = zeros(length(sessions),1);
for ss = 1:length(sessions)
    in = strcmp(sessionType,sessions{ss}) & r_squared >= r2_thr;
    pxn = pxn_aligned(in);
    [~,tbinL] = cellfun(@size,pxn);
    pxn_combined = NaN(length(posaxis),min(tbinL),length(pxn));
    for pp = 1:length(pxn)
        % align offset of replay
        range = size(pxn{pp},2)-min(tbinL)+1:size(pxn{pp},2);
        pxn_combined(:,:,pp) = pxn{pp}(:,range);
    end
    Out.(sessions{ss}) = pxn_combined;
    
    ha(ss) = subplot(2,length(sessions),ss);
    taxis = 0-(min(tbinL)-1)*timeStep:timeStep:0;
    imagesc(taxis,(posaxis-pi)*180/pi,nanmean(pxn_combined,3))
    axis xy square
    hold on
    plot(xlim,[0 0],'w--','LineWidth',1)
    cb = colorbar;
    title([sessions{ss},' n=',num2str(sum(in))])
    if ss == 1
        xlabel('Time relative to replay offset (sec)')
        ylabel('Position relative to reward (deg)')
    elseif ss == length(sessions)
        ylabel(cb,'Mean posterior probability')
    end
    
    % Add up probability near possible goal location 
    for pp = 1:length(locationRangeInd)
        inrange = locationRangeInd(pp)-1:locationRangeInd(pp)+1;
        ProbSum(pp,ss) = nansum(nansum(nansum(pxn_combined(inrange,:,:))));
    end
    [~,pvalue(ss)] = lillietest(ProbSum(:,ss));
    subplot(2,length(sessions),ss+length(sessions))
    bar(-nLocAway:nLocAway,ProbSum(:,ss))
    axis square
    colorbar
    if ss == 1
        xlabel('Location # relative to goal')
        ylabel('sum of posterior probability')
    end
end
colormap hot
set(ha,'CLim',[0 .04])
saveas(h,strcat(outdir,'\GoalLoactionBias'),'epsc')
saveas(h,strcat(outdir,'\GoalLoactionBias'),'png')
close(h)
% test for outlier
[tf]= isoutlier(ProbSum,'gesd');
t = array2table(tf);
t2 = array2table(pvalue);
t.Properties.VariableNames=sessions;
t2.Properties.VariableNames=sessions;
Stats.ReplayOffset.OutlierTest = t;
Stats.ReplayOffset.Lillietest = t2;

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
saveas(h,strcat(outdir,'\r2overTrials'),'epsc')
saveas(h,strcat(outdir,'\r2overTrials'),'png')
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
    pxn_correct = p_x_n(Correct_test & in);
    stopLoc_correct = stopLoc(Correct_test & in);
    rewardLoc_correct = rewardLoc(Correct_test & in);
    r2_correct = r_squared(Correct_test & in);
    % incorrect before trials
    pxn_incorrect = p_x_n(~Correct_test & in);
    stopLoc_incorrect = stopLoc(~Correct_test & in);
    rewardLoc_incorrect = rewardLoc(~Correct_test & in);
    r2_incorrect = r_squared(~Correct_test & in);
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
    saveas(h,strcat(outdir_figs,'\BayesianDecodingExamples'),'epsc')
    saveas(h,strcat(outdir_figs,'\BayesianDecodingExamples'),'png')
    close(h)
    
    %% Plot r2 distributions between correct and incorrect trials
    if ~strcmp(nType{st},'PreRun')
        in = strcmp(sessionType,nType{st});
        r2min = min(r_squared);
        r2max = max(r_squared);
        sigma = .1;
        [Wdis_correct,axisValue1] = WeightedProportion(r_squared(Correct_test & in),r2min,r2max,sigma);
        sample = bootstrp(5000,@WeightedProportion,r_squared(Correct_test & in),r2min,r2max,sigma);
        [CI_u, CI_l] = CorrectionForMultipleComparsion(sample);
        CI_correct = [CI_u'-Wdis_correct;Wdis_correct-CI_l'];

        [Wdis_incorrect,axisValue2] = WeightedProportion(r_squared(~Correct_test & in),r2min,r2max,sigma);
        sample = bootstrp(5000,@WeightedProportion,r_squared(~Correct_test & in),r2min,r2max,sigma);
        [CI_u, CI_l] = CorrectionForMultipleComparsion(sample);
        CI_incorrect = [CI_u'-Wdis_incorrect; Wdis_incorrect-CI_l'];

        h = figure; hold on
        ha1 = shadedErrorBar(axisValue1,Wdis_correct,CI_correct,'r-');
        ha2 = shadedErrorBar(axisValue2,Wdis_incorrect,CI_incorrect,'k-');
        ylabel('Proportion')
        xlabel('r^2')
        legend([ha1.mainLine; ha2.mainLine],'String',...
            {strcat('Correct(n=',num2str(sum(Correct_test & in)),')'); strcat('Incorrect(n=',num2str(sum(~Correct_test & in)),')')})
        saveas(h,strcat(outdir_figs,'\R2distCorrectVSIncorrect'),'epsc')
        saveas(h,strcat(outdir_figs,'\R2distCorrectVSIncorrect'),'png')
        close(h)
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
    RippleOnsetIndex,RippleOffsetIndex,RippleFR,ActiveNearReward,ActiveNearRest,RippleFR_xlabel,smspkFR,cellID,param] = ...
    BayesianDecoding(spklist,encodingfd,posIndfile,eeg2posIndfile,spikeIndfile,trackdata,RippleOnsetIndex,RippleOffsetIndex)
% Example inputs
% parentfd = 'C:\CA1Replayproject\Data\Rat139\circulartrack\2017-02-20-CT-2';
% encodingfd = [1]; (1 is Prerunning; 2 is sample; 3 is test; 4 is TrialsPost)
% posIndfile = strcat(parentfd,'\Data_angle_ontrack.mat');
% eeg2posIndfile = strcat(parentfd,'\SpikeTime\eeg2pos.mat');
% spikeIndfile = strcat(parentfd,'\SpikeTime\spike_index.mat');

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

%% Extract time index in rest box
[AfterTest,BeforeTest,PostTest,PreRun]= ExtractTimeInd_RestBox(eeg2posIndfile,posIndfile);

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
NTrial = zeros(length(RippleOnsetIndex),1);
NTrial_post = zeros(length(RippleOnsetIndex),1);
NTrial_preRun = zeros(length(RippleOnsetIndex),1);
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
    if sum(~isnan(BeforeTest(InRange))) == length(InRange) || sum(~isnan(AfterTest(InRange))) == length(InRange)...
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
            if sum(~isnan(AfterTest(InRange))) == length(InRange) || sum(~isnan(BeforeTest(InRange))) == length(InRange)
                if mode(AfterTest(InRange))<=4 || mode(BeforeTest(InRange))<=4
                   sessionType{ii} = 'InitialFour';
                   ntrial = [mode(AfterTest(InRange)) mode(BeforeTest(InRange))];
                   NTrial(ii) = ntrial(~isnan(ntrial));
                   ripple_time(:,1) = ripple_time(:,1) + length(InRange)/sampFreq;
                   ripple_spkcount(:,1) = ripple_spkcount(:,1) + spkcount;
                   RippleFR_xlabel{1} = sessionType{ii};
                elseif mode(AfterTest(InRange))>4 || mode(BeforeTest(InRange))>4
                   sessionType{ii} = 'LastFour';
                   ntrial = [mode(AfterTest(InRange)) mode(BeforeTest(InRange))];
                   NTrial(ii) = ntrial(~isnan(ntrial));
                   ripple_time(:,2) = ripple_time(:,2) + length(InRange)/sampFreq;
                   ripple_spkcount(:,2) = ripple_spkcount(:,2) + spkcount;
                   RippleFR_xlabel{2} = sessionType{ii};
                end
            elseif sum(~isnan(PostTest(InRange))) == length(InRange)
               sessionType{ii} = 'PostTest';
               NTrial_post(ii) = mode(PostTest(InRange));
               ripple_time(:,3) = ripple_time(:,3) + length(InRange)/sampFreq;
               ripple_spkcount(:,3) = ripple_spkcount(:,3) + spkcount;
               RippleFR_xlabel{3} = sessionType{ii};
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
NTrial = NTrial(ripple_in);
NTrial_post = NTrial_post(ripple_in);
NTrial_preRun = NTrial_preRun(ripple_in);
RippleOnsetIndex = RippleOnsetIndex(ripple_in);
RippleOffsetIndex = RippleOffsetIndex(ripple_in);
RippleFR = ripple_spkcount./ripple_time;
param.posbinsize = param.posbinsize/360*2*pi; % convert to radian

function [AfterTest,BeforeTest,PostTest,PreRun]= ExtractTimeInd_RestBox(eeg2posIndfile,posIndfile)
load(eeg2posIndfile);
load(posIndfile);

nTestTrial = length(data_angle{3});
BeforeTest = NaN(size(eeg2pos));
AfterTest = NaN(size(eeg2pos));
for t = 1:nTestTrial
    % before test in rest box
    TestTrial_onset = data_angle{3}{t}(1,1);
    SampleTrial_offset = data_angle{2}{t}(end,1);
    % position index
    ind_range = [find((data_angle_all(:,1)-SampleTrial_offset)== 0),...
                 find((data_angle_all(:,1)-TestTrial_onset)== 0)];
    % translate into eeg index
    eegind_range = [find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1),...
                    find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1,'last')];
    BeforeTest(min(eegind_range):max(eegind_range)) = t;
    
    % after test in rest box
    TestTrial_offset = data_angle{3}{t}(end,1);
    if t < nTestTrial        
        SampleTrial_onset = data_angle{2}{t+1}(1,1);
    else
        SampleTrial_onset = data_angle{4}{1}(1,1); % onset from PostTrial
    end
     % position index
    ind_range = [find((data_angle_all(:,1)-SampleTrial_onset)== 0),...
                 find((data_angle_all(:,1)-TestTrial_offset)== 0)];
    % translate into eeg index
    eegind_range = [find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1),...
                    find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1,'last')];
    % exclude the period after last test trial due to long wait period
    % before post running initiates
    if t < nTestTrial
        AfterTest(min(eegind_range):max(eegind_range)) = t;
    end
end

nPrerunningTrial = length(data_angle{1});
PreRun = NaN(size(eeg2pos));
for t = 1:nPrerunningTrial
    % in between pre-run in rest box, except the last pre-run that include
    % time until the first sample trial onset
    PrerunningTrial_offset = data_angle{1}{t}(end,1);
    if t < nPrerunningTrial
        PrerunningTrial_onset = data_angle{1}{t+1}(1,1);
    else
        PrerunningTrial_onset = data_angle{2}{1}(1,1); % first sample trial onset
    end
    % position index
    ind_range = [find((data_angle_all(:,1)-PrerunningTrial_offset)== 0),...
                 find((data_angle_all(:,1)-PrerunningTrial_onset)== 0)];
    % translate into eeg index
    eegind_range = [find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1),...
                    find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1,'last')];
    PreRun(min(eegind_range):max(eegind_range)) = t;
end

nPostTestTrial = length(data_angle{4});
PostTest = NaN(size(eeg2pos));
for t = 1:nPostTestTrial-1
    % in between post test run in rest box
    PostTestTrial_onset = data_angle{4}{t+1}(1,1);
    PostTestTrial_offset = data_angle{4}{t}(end,1);
    % position index
    ind_range = [find((data_angle_all(:,1)-PostTestTrial_offset)== 0),...
                 find((data_angle_all(:,1)-PostTestTrial_onset)== 0)];
    % translate into eeg index
    eegind_range = [find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1),...
                    find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1,'last')];
    PostTest(min(eegind_range):max(eegind_range)) = t+1;
end

checkoverlap = ~isnan([BeforeTest AfterTest PreRun PostTest]);
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