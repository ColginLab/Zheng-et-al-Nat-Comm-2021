%% This version add position range to include spikes and add theta phase
%  and use only the first spike in a sequnece event
function group_spkphase_approachReward_v3(DownSample)
% DownSample: true or false; downsample correct trials
parentfd = fileparts(mfilename('fullpath'));

powthr = NaN; % in z score; Set to NaN for no threshold
SpkInSequnece = true; % only include spikes within sequenct events?
PosIncluded = [-10 -5]; % in location #
downsample_str = [];
if DownSample
    downsample_str = '_ds';
end
file_input1 = 'Data_gammaPhs_cyclenum.mat';
file_input2 = 'Data_gammaBandpass.mat';
%file_input3 = [parentfd,'\GroupData\group_gammaTFR_eachseq_20190529.mat'];
file_input3 = 'E:\ColginLab\Data Analysis\GroupData\group_gammaTFR_eachseq_20190529.mat';
file_input_cell = 'Cells_ds_ALLLaps_v2_vel_0.mat';
%file_output = [parentfd,'\GroupData Figures\group_spkphase_v3',downsample_str,'.mat'];
file_output = ['E:\ColginLab\Data Analysis\GroupData\group_spkphase_v3',downsample_str,'.mat'];
%fig_folder_out = [parentfd,'\GroupData Figures\group_spkphase_v3'];
fig_folder_out = 'E:\ColginLab\Figures\Figure5';
directories_allData_v1

numbins = 90; % Number of bins
bin_ang = 2*pi/numbins;
ang_bins = bin_ang/2:bin_ang:(2*pi-bin_ang/2);
ncell = 3;
nspk = 5;
nmove = 1;
% Load sequence data
seq = load(file_input3);

pos_prdc = nan(seq.Nseq,2);
pos_real = nan(seq.Nseq,2);
dist_prdc_real = nan(seq.Nseq,1);
vel = nan(seq.Nseq,1);
for nseq = 1:seq.Nseq
    min0 = min(seq.Pxn_all{nseq,4});
    max0 = max(seq.Pxn_all{nseq,4});
    pos_prdc(nseq,:) = [min0,max0];
    
    min0 = min(seq.Pxn_all{nseq,6});
    max0 = max(seq.Pxn_all{nseq,6});
    pos_real(nseq,:) = [min0,max0];
    
    dist_prdc_real(nseq) = min(abs(angle(exp((seq.Pxn_all{nseq,4}'-seq.Pxn_all{nseq,6})*1i))));
    vel(nseq,1) = mean(seq.Pxn_all{nseq,7});
end

ind_seq_log = seq.data_info(:,4) <= 8 & ...
    pos_real(:,2) < seq.data_info(:,8) &... % stop location
    seq.data_info(:,10) >= ncell & seq.data_info(:,11) >= nspk & ...  % limit active cells number
    seq.data_info(:,6) >= 0 &...  % only correct or error trials
    seq.data_info(:,12) > 0 &... % slope    
    seq.data_info(:,13) >= .6 &... % replay score    
    seq.para_all(:,1) <= 20 & seq.para_all(:,4) >= nmove &...
    dist_prdc_real <= 5*bin_ang;

Phasebin = 10; % in degree
PhaseAxis = 0:Phasebin:360;
numGroup = 3; % Early, middle, and late
ConditionType = {'SampleCorrect','SampleError','TestCorrect','TestError'};
SpkPhasePro_SG = zeros(length(PhaseAxis),numGroup,length(ConditionType));
SpkPhasePro_FG = zeros(length(PhaseAxis),numGroup,length(ConditionType));
SpkPhasePro_theta = zeros(length(PhaseAxis),numGroup,length(ConditionType));
CellGroup=[]; conditionType=[]; sessionID=[]; timestamp=[]; PhaseAngle=[];
SpkPhaseAngle_theta = table(CellGroup,conditionType,sessionID,timestamp,PhaseAngle);
for ns = 1:isession
    cd(pathRats{ns});
    disp(['Currentrly in ' pathRats{ns}])
    if exist(file_input1,'file') ~= 2
        continue
    end
    load(trackdata{ns})
    load(file_input1)
    load(file_input2)
    load(file_input_cell)    
    s = RandStream('mt19937ar','Seed',ns);
    
    % Use the tetrode with highest theta power to estimate gamma phase
    tet_gamma = CSClist_CA1_sort{ns}(1);
    
    % Obtain time for phase estimate
    phase_t = Phase{1};
    
    % Get time index within sequence events
    if SpkInSequnece
        in = seq.data_info(:,2) == ns;
        temp = [cellfun(@min,seq.Pxn_all(in & ind_seq_log,1)),... 
                       cellfun(@max,seq.Pxn_all(in & ind_seq_log,1))];
        [~,sortInd] = sort(temp(:,1));
        t_startstop = temp(sortInd,:);
        InSequence = false(size(phase_t));
        SequenceID = zeros(size(phase_t));
        for ss = 1:size(t_startstop,1)
            InRange = phase_t >= t_startstop(ss,1) & phase_t <= t_startstop(ss,2);
            InSequence(InRange) = true;
            SequenceID(InRange) = ss;
        end
        InSequence = find(InSequence);
        SequenceID(SequenceID==0)=[];
    else
        InSequence = 1:length(phase_t);
    end
    
    % Get time index when gamma power exceed threshold
    if ~isnan(powthr)
        SGpower_in = find(BP{2,4} >= powthr);
        FGpower_in = find(BP{2,3} >= powthr);
    else
        SGpower_in = 1:length(phase_t);
        FGpower_in = 1:length(phase_t);
    end
    
    % Get time index when animals are within position range
    locationRange = mean(diff(Ang_RewardLoc_ontrack));
    rewardpos = Ang_RewardLoc_ontrack(ind_rewardloc);
    posRange_spk = PosIncluded*locationRange+rewardpos;

    % Assign condition type index
    Condition_ind = zeros(size(phase_t));
    for cc = 1:length(ConditionType)
        if ~isempty(strfind(ConditionType{cc},'Sample'))
            tRange = [ts_sample_start ts_sample_reward]/1e6;            
        elseif ~isempty(strfind(ConditionType{cc},'Test'))
            tRange = [ts_test_start ts_test_reward]/1e6;
        end
        
        if ~isempty(strfind(ConditionType{cc},'Correct'))
            outcome = sign_correct_test==1;            
        elseif ~isempty(strfind(ConditionType{cc},'Error'))
            outcome = sign_correct_test==0;
        end
        
        if DownSample
            trialNumDiff = sum(sign_correct_test==1)-sum(sign_correct_test==0);
            if trialNumDiff > 0
                correctTrials = find(sign_correct_test==1);
                ind = randperm(s,length(correctTrials),trialNumDiff);
                outcome(correctTrials(ind)) = false;
            elseif trialNumDiff < 0
                errorTrials = find(sign_correct_test==0);
                ind = randperm(s,length(errorTrials),abs(trialNumDiff));
                outcome(errorTrials(ind)) = false;
            end
        end
        for ll = 1:length(outcome)
            if outcome(ll)
                InRange = phase_t>=tRange(ll,1) & phase_t<=tRange(ll,2);
                Condition_ind(InRange) = cc;
            end
        end
    end
    % Use peak FR to sort cells into Early, Middle, and Late position
    % while animals approach to a reward location
    posRange = [Ang_StartZone_depart_ontrack, Ang_RewardLoc_ontrack(ind_rewardloc)];
    posBoundary = posRange(1):diff(posRange)/numGroup:posRange(2);
    [~,ind] = max(Ratemap_AllLaps);
    peakFRpos = mapAxis(ind);
    [~,groupInd] = histc(peakFRpos,posBoundary);
    for ii = 1:numGroup
        in = groupInd==ii;
        spikes_in = spikes(in,1:2);
        spk = [cat(1,spikes_in{:,1}), cat(1,spikes_in{:,2})];
        if isempty(spk)
            continue
        end
        
        if isnan(posRange_spk)
            spkt = spk(:,1);
            spkt = sort(spkt);
        else
            in_spk = spk(:,2)>=min(posRange_spk) & spk(:,2)<=max(posRange_spk);
            spkt = spk(in_spk,1);
            [spkt,sortInd] = sort(spkt);
            in_spk = in_spk(sortInd);
        end
        
        spk2phase_ind = interp1(phase_t,1:length(phase_t),spkt,'nearest','extrap');
        tdiff = abs(phase_t(spk2phase_ind)-spkt');
        rmv = tdiff>1/30; % remove spikes that are too far away
        spk2phase_ind(rmv) = [];
        
        if ~isempty(spk2phase_ind)
            for cc = 1:length(ConditionType)
                [in] = ismember(spk2phase_ind,find(Condition_ind==cc));
                [in_gamma] = ismember(spk2phase_ind,FGpower_in);
                [in_sequence,ind] = ismember(spk2phase_ind,InSequence);
                spkind = find(in_sequence);
                spkSeqID = SequenceID(ind(ind~=0));
                [~,firstspkInd] = unique(spkSeqID);
                spkind = spkind(firstspkInd);
                in_sequence = false(size(in_sequence));
                in_sequence(spkind) = true;
                
                count = histc(Phase{1,3}(tet_gamma,spk2phase_ind(in & in_gamma & in_sequence)),PhaseAxis); % FG
                SpkPhasePro_FG(:,ii,cc) = SpkPhasePro_FG(:,ii,cc) + count';
                
                in_gamma = ismember(spk2phase_ind,SGpower_in);
                count = histc(Phase{1,4}(tet_gamma,spk2phase_ind(in & in_gamma & in_sequence)),PhaseAxis); % SG
                SpkPhasePro_SG(:,ii,cc) = SpkPhasePro_SG(:,ii,cc) + count';
                
                count = histc(Phase{1,2}(1,spk2phase_ind(in & in_sequence)),PhaseAxis); % Theta
                SpkPhasePro_theta(:,ii,cc) = SpkPhasePro_theta(:,ii,cc) + count';
                
                PhaseAngle = Phase{1,2}(1,spk2phase_ind(in & in_sequence))';
                CellGroup = repmat(ii,size(PhaseAngle));
                conditionType = repmat(ConditionType(cc),size(PhaseAngle));
                sessionID = repmat(ns,size(PhaseAngle));
                timestamp = Phase{1,1}(1,spk2phase_ind(in & in_sequence))';
                Tnew = table(CellGroup,conditionType,sessionID,timestamp,PhaseAngle);
                SpkPhaseAngle_theta = [SpkPhaseAngle_theta; Tnew];
            end
        end
    end
end
SpkPhasePro_SG = squeeze(sum(SpkPhasePro_SG,2));
SpkPhasePro_FG = squeeze(sum(SpkPhasePro_FG,2));
SpkPhasePro_theta = squeeze(sum(SpkPhasePro_theta,2));
SpkPhasePro_SG_norm = SpkPhasePro_SG./repmat( sum(SpkPhasePro_SG),length(PhaseAxis),1 );
SpkPhasePro_FG_norm = SpkPhasePro_FG./repmat( sum(SpkPhasePro_FG),length(PhaseAxis),1 );
SpkPhasePro_theta_norm = SpkPhasePro_theta./repmat( sum(SpkPhasePro_theta),length(PhaseAxis),1 );

% Smooth data
kernel = gausskernel(5,3);
for cc = 1:length(ConditionType)
    SpkPhasePro_SG_norm(:,cc) = conv_cir(SpkPhasePro_SG_norm(:,cc),kernel);
    SpkPhasePro_FG_norm(:,cc) = conv_cir(SpkPhasePro_FG_norm(:,cc),kernel);
    SpkPhasePro_theta_norm(:,cc) = conv_cir(SpkPhasePro_theta_norm(:,cc),kernel);
end

% SpkPhasePro_SG_norm = SpkPhasePro_SG./repmat( sum(SpkPhasePro_SG),length(PhaseAxis),1,1 );
% SpkPhasePro_FG_norm = SpkPhasePro_FG./repmat( sum(SpkPhasePro_FG),length(PhaseAxis),1,1 );
% SpkPhasePro_theta_norm = SpkPhasePro_theta./repmat( sum(SpkPhasePro_theta),length(PhaseAxis),1,1 );
% 
% % Smooth data
% kernel = gausskernel(5,3);
% for cc = 1:length(ConditionType)
%     for c2 = 1:size(SpkPhasePro_SG_norm,2)
%         SpkPhasePro_SG_norm(:,c2,cc) = conv_cir(SpkPhasePro_SG_norm(:,c2,cc),kernel);
%         SpkPhasePro_FG_norm(:,c2,cc) = conv_cir(SpkPhasePro_FG_norm(:,c2,cc),kernel);
%         SpkPhasePro_theta_norm(:,c2,cc) = conv_cir(SpkPhasePro_theta_norm(:,c2,cc),kernel);
%     end
% end

%% Do stats
numGroup = 1;
pvalue = zeros(numGroup,1);
pvalue_samp = zeros(numGroup,1); % between correct and error in sample trials
pvalue_test = zeros(numGroup,1); % between correct and error in test trials
pvalue_combined = zeros(numGroup,1); % between correct and error in all trials
medianAng_combined = zeros(numGroup,2);
SpkCount_combined = zeros(numGroup,2);
for ii = 1:numGroup
    in = true(size(SpkPhaseAngle_theta.PhaseAngle));
    %in=SpkPhaseAngle_theta.CellGroup==ii;
    alpha=SpkPhaseAngle_theta.PhaseAngle(in);
    idx = zeros(sum(in),1);
    ctype = SpkPhaseAngle_theta.conditionType(in);
    for cc = 1:length(ConditionType)
        in2 = strcmp(ctype,ConditionType{cc});
        idx(in2)=cc;
    end
    pvalue(ii) = circ_cmtest(alpha/180*pi,idx);
    in3 = idx==1 | idx==2;
    pvalue_samp(ii) = circ_cmtest(alpha(in3)/180*pi,idx(in3));
    in4 = idx==3 | idx==4;
    pvalue_test(ii) = circ_cmtest(alpha(in4)/180*pi,idx(in4));
    idx_crtErr = zeros(size(idx));
    idx_crtErr(idx==1 | idx==3) = 1; idx_crtErr(idx==2 | idx==4) = 2;
    pvalue_combined(ii) = circ_cmtest(alpha/180*pi,idx_crtErr);
    medianAng_combined(ii,:) = [circ_median(alpha(idx_crtErr==1)/180*pi), circ_median(alpha(idx_crtErr==2)/180*pi)];
    SpkCount_combined(ii,:) = [sum(idx_crtErr==1), sum(idx_crtErr==2)];
end
StatsOut = table(pvalue_combined,medianAng_combined,SpkCount_combined);
save(file_output,'SpkPhasePro_SG','SpkPhasePro_FG','SpkPhasePro_theta','SpkPhaseAngle_theta',...
                 'SpkPhasePro_SG_norm','SpkPhasePro_FG_norm','SpkPhasePro_theta_norm',...
                 'StatsOut',...
                 'ConditionType','PhaseAxis','numGroup');

%% Plot results
% Shift phase
shift = floor(length(PhaseAxis)/2);
PhaseAxis = PhaseAxis-PhaseAxis(shift);
% SpkPhasePro_SG_norm = circshift(SpkPhasePro_SG_norm,shift);
% SpkPhasePro_FG_norm = circshift(SpkPhasePro_FG_norm,shift);
SpkPhasePro_theta_norm = circshift(SpkPhasePro_theta_norm,shift);
h1 = figure('Units','normalized','Position',[0 0 1/2 1/2]);
imagesc(1:length(ConditionType),PhaseAxis(1:end-1),squeeze(SpkPhasePro_theta_norm(1:end-1,:,:)))
ylabel(sprintf('Theta spike\nphase (degrees)'))
axis xy square
set(gca,'XTick',1:length(ConditionType),'XTickLabel',ConditionType)
set(gca,'CLim',[0 .05])
colorbar

if ~isdir(fig_folder_out)
    mkdir(fig_folder_out)
end

saveas(h1,[fig_folder_out,'\SpkPhaseLock_AcrossTrack',downsample_str],'fig')
saveas(h1,[fig_folder_out,'\SpkPhaseLock_AcrossTrack',downsample_str],'epsc')
close(h1)