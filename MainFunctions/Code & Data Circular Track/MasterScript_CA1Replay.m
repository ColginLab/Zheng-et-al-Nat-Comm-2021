%% Master script %%
parentfd = fileparts(mfilename('fullpath'));
ind = strfind(parentfd,'MATLAB');
%Datafd = [parentfd(1:2),'Data\Circular Track'];     % this should point to data directory
Datafd = 'E:\ColginLab\Data\Circular Track\';

FigDir{1} = 'E:\ColginLab\Figures\Figure6';
FigDir{2} = 'E:\ColginLab\Figures\Figure7';
FigDir{3} = 'E:\ColginLab\Figures\SupFigure5';
FigDir{4} = 'E:\ColginLab\Figures\SupFigure9';

directories_allData_v1;
if verLessThan('matlab','8.4')
    error('Please use Matlab version 2014b or above')
end
%% Analysis output directory
%Outdir = strcat(parentfd,'\Analysis_EH_ds');
Outdir = 'E:\ColginLab\Data Analysis\Analysis_EH_ds';
if ~isdir(Outdir)
    mkdir(Outdir)
end
%% Cell list to use
cellList = 'TTList_dCA1_pyr_ds.txt';
spktime_fd = 'SpikeTime_ds';
% cellList = 'TTList_dCA1_pyr.txt';
% spktime_fd = 'SpikeTime';

%% Preprocessing: Match spkie timestamp with EEG timestamp
Preprocess_done = true;
if ~Preprocess_done
    doall_data_pos2ang;
    for ff = 1:length(pathRats)
        MatchSpikeTimeToEEGTime_v3(pathRats{ff},spktime_fd,cellList);
    end
end

%% Behavior performance
%BehaviorPerformance(trackdata,Outdir);

%% Bayesian decoding analysis
plotOnly = false;
%RunConfusionMatrix_AcrossTrials_v2(pathRats,trackdata,Outdir,cellList,spktime_fd);          % % WORKS
%RunConfusionMatrix_AcrossTrials_v3(pathRats,trackdata,Outdir,cellList,spktime_fd);          % % WORKS
% RunBayesianDecoder_v3(path,CSClist_CA1,Outdir,plotOnly); % use population firing rate to find trajectory events
RunBayesianDecoder_ripple(pathRats,CSClist_CA1,Outdir,cellList,spktime_fd,trackdata,plotOnly,FigDir); %use position tunning after reward presentation    % % WORKS
% RunBayesianDecoder_ripple_v2(path,CSClist_CA1,Outdir,cellList,spktime_fd,trackdata,plotOnly); %use position tunning before and after reward presentation
%RunBayesianDecoder_ripple_v3(pathRats,CSClist_CA1,Outdir,cellList,spktime_fd,trackdata,plotOnly); % Look at ripples outside of rest box          % % WORKS
%RunBayesianDecoder_ripple_v3_afterReward(pathRats,CSClist_CA1,Outdir,cellList,spktime_fd,trackdata,plotOnly);           % % WORKS

%RunBayesianDecoder_ripple_PreRunTuning(pathRats,CSClist_CA1,Outdir,cellList,spktime_fd,trackdata,plotOnly);         % % WORKS

%RunBayesianDecoder_ripple_BeforeTest(pathRats,CSClist_CA1,Outdir,cellList,spktime_fd,trackdata,plotOnly); % include only the time between sample and test during learning phase  % % WORKS
% RunBayesianDecoder_ripple_ShiftedPlaceTuning(path,CSClist_CA1,Outdir,plotOnly)
%ReplayQualityAcrossConditions(Outdir);                  % % WORKS
% plotSequenceEvents(path,trackdata,Outdir);

%% Place map analysis
%RunPlaceCellPropertiesAcrossLap_v4(path,trackdata,Outdir); % Estimate firing rate with kernel density estimate
%RunPlaceCellPropertiesAcrossLap_v5(pathRats,trackdata,Outdir,cellList,spktime_fd);          % % WORKS
%RunPlaceCellPropertiesAcrossLap_v6(pathRats,trackdata,Outdir,cellList,spktime_fd); % Estimate firing rate with histogram % % WORKS
%CompareFRrelativeToReward(pathRats,trackdata,Outdir,spktime_fd);                            % % WORKS
%PlaceFieldsDist(trackdata,Outdir,Outdir);                                                   % % WORKS
%PlaceFieldsDist_alignedToRestBox(Outdir)                                                    % % WORKS
%PlaceTuningAcrossLap(Datafd,Outdir);                                                        % % WORKS
%RunSpeedDist(pathRats,Outdir);                                                              % % WORKS
% CorrErrPopVecCorr(Outdir);
% CorrErrPopVecCorr_v2(trackdata,Outdir); % align to stop location
%CorrErrPopVecCorr_v9(trackdata,Outdir); % separate cell types                                % % WORKS
%CorrErrPopVecCorr_v10(trackdata,Outdir); % separate cell types and align to stop location    % % WORKS
%CorrErrPopVecCorr_v7(trackdata,Outdir); % include all cells                                  % % WORKS
%CorrErrPopVecCorr_v11(trackdata,Outdir); % include all cells and align to stop location      % % WORKS


%% Ripple Analysis
RippleTrigFR_reward(Datafd,Outdir,spktime_fd);                              % % WORKS
RippleTrigFR_reward_v2(Datafd,Outdir,spktime_fd);                           % % WORKS
RippleTrigFR_reward_PreRunPosTuning(Datafd,Outdir,spktime_fd);              % % WORKS
CellPairInRipple(Datafd,Outdir);                                            % % WORKS

%% LFP analysis
% PVCorr_VFR(path,trackdata,CSClist_CA1,Outdir);
% PFRAcrossConditions(path,trackdata,CSClist_CA1,Outdir);
% SpectDiff_CortErr(path,trackdata,CSClist_CA1,Outdir);