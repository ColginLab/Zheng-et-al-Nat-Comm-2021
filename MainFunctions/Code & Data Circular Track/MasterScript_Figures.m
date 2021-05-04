%% Master script to produce figures
% Figures are stored in the following directory: 
% ~\MATLAB\MainFunctions\Code & Data Circular Track\GroupData Figures
% Figures made for publication are stored in the folowing directory:
% ~\MATLAB\MainFunctions\Code & Data Circular Track\Figures_manuscript

%% Set up dependencies and path
clearvars
close all
clc

parentfd = pwd;
depDir = [parentfd(1:end-41) '\Dependencies'];
genFDir = [parentfd(1:end-41) '\GeneralFunctions'];
codeDataSeqDir = [parentfd(1:end-41) '\MainFunctions\Code & Data Sequence'];
addpath([depDir '\matbugs'])
addpath(genFDir)
addpath(codeDataSeqDir)

%Code & Data Sequence

%%-----------------------------------Main figures------------------------------------------%%

%% Figure 1: Behavioral performance
clear
%group_behavior_statespace_all  %%
%code_behavior_v1               %%
% Relevant figures: Behavior_barGraph, 
%                   Behavior_learningCurve

%% Figure 2: Predictive firing developed with learning
clear
doall_group_Pxn
Stats_PxnSum
% Relevant figures: PxnSum_line_ahead,
%                   PxnSum_line_behind                  
% Relevant MAT files: Stats_PxnSum.mat

%% Figure 3: Difference in the slope of place cell sequence between correct and error test trials
clear
doall_get_gammaTFR_eachseq_EH
Stats_gammaTFR_eachseq_forMultiCompare
BayesianDecodingExample
% Relevant figures: SeqSlope_test, 
%                   SeqExample_test_Rat149_2017-08-08-CT-2
% Relevant MAT files: Stats_SeqSlope.mat

%% Figure 4: Difference in the slope of place cell sequence between correct and error sample trials
% figures are already produced by the function
%   "Stats_gammaTFR_eachseq_forMultiCompare" and "BayesianDecodingExample"
% Relevant figures: SeqSlope_samp,
%                   SeqExample_sample_Rat149_2017-08-08-CT-2

%% Figure 5: Theta phase locking
clear
doall_get_gammaTFR_eachseq_EH
group_spkphase_approachReward_v2(false)
% Relevant figures: SpkPhaseLock_AcrossTrack
% Relevant MAT files: group_spkphase_v2.mat

%% Figure 6: No difference in replay fidelity between correct and error trials
clear
parentfd = ReturnParentPath;
directories_allData_v1;
cellList = 'TTList_dCA1_pyr_ds.txt';
spktime_fd = 'SpikeTime_ds';
% Analysis output directory
Outdir = strcat(parentfd,'\Analysis_EH_ds\');
if ~isdir(Outdir)
    MasterScript_CA1Replay % Run Ernie's script
end
% Plot analysis result again in the same folder as other figures
Outdir2 = [parentfd,'\GroupData Figures'];
RunBayesianDecoder_ripple(path,CSClist_CA1,Outdir2,cellList,spktime_fd,trackdata,true);
% Relevant figures: BayesianDecoding_ripple\R2distCorrectVSIncorrect, 
%                   BayesianDecoding_ripple\BayesianDecodingExamples

%% Figure 7: Termination bias at goal location develops with learning
% figures are already produced by the function "RunBayesianDecoder_ripple"
% Relevant figures: BayesianDecoding_ripple\TerminationBias_goalPos_Error, 
%                   BayesianDecoding_ripple\TerminationBias_goalPos_Correct

%% Figure 8: Gamma phase coding is intact in different types of trials
clear
group_gammaphase_eachseq
Stats_spkphaseLock
% Relevant figures: SlowG_spkPhaseLock, 
%                   FastG_spkPhaseLock,
%                   Stats_spkphaseLock
% Relevant MAT files: Stats_spkphaseLock.mat

%%
%%--------------------------------Supplementary figures-------------------------------------%%

%% Supplementary figure 1: Over-representation of reward location by place cells
clear
parentfd = ReturnParentPath;
directories_allData_v1;
cellList = 'TTList_dCA1_pyr.txt';
spktime_fd = 'SpikeTime';
Outdir = strcat(parentfd,'\Analysis_EH\');

Preprocess_done = true;
if ~Preprocess_done
    for ff = 1:length(path)
        MatchSpikeTimeToEEGTime_v3(path{ff},spktime_fd,cellList);
    end
end

if ~isdir(Outdir)
    if verLessThan('matlab','8.4')
        error('Please use Matlab version 2014b or above')
    else
        RunPlaceCellPropertiesAcrossLap_v6(path,trackdata,Outdir,cellList,spktime_fd); 
    end
end
% Plot analysis result in the same folder as other figures
Outdir2 = strcat(parentfd,'\GroupData Figures\');
PlaceFieldsDist(trackdata,Outdir,Outdir2);

% Do cross validation
InputFile = strcat(parentfd,'\Analysis_EH_ds\ConfusionMatrix\ConfusionMat.mat');
if exist(InputFile,'file')~=2
    MasterScript_CA1Replay % Run Ernie's script
end
Outdir = strcat(parentfd,'\GroupData Figures\');
PlotCrossValidation(InputFile,Outdir);
% Relevant figures: PlaceFieldsDist_alignedToReward\PlaceFieldDis_after,
%                   CrossValidation

%% Supplementary figure 2: Place cell sequence with negative slope
clear
doall_get_gammaTFR_eachseq_EH_v2
Stats_gammaTFR_eachseq_forMultiCompare_rev
% Relevant figures: SeqSlopeNeg_average_samp, 
%                   SeqSlopeNeg_average_test
% Relevant MAT files: Stats_SeqSlopeNeg.mat

%% Supplementary figure 3: More Bayesian decoding examples
% figures are already produced by the function "BayesianDecodingExample"
% Relevant figures: SeqExample_test_Rat149_2017-08-08-CT-1, 
%                   SeqExample_test_Rat148_2017-12-02-CT

%% Supplementary figure 4: Termination bias at goal location during Pre-run and Post-test
% figures are already produced by the function "RunBayesianDecoder_ripple"
% Relevant figures: BayesianDecoding_ripple\TerminationBias_goalPos_Error, 
%                   BayesianDecoding_ripple\TerminationBias_goalPos_Correct

%% Supplementary figure 5: No difference in gamma power during different types of trials
clear
doall_PFR
% one figure is already produced by the function "Stats_gammaTFR_eachseq_forMultiCompare"
% Relevant figures: SampTest_CrtErr_CI_SR, 
%                   SampTest_CrtErr_SR, 
%                   SeqSlopePower
% Relevant MAT files: PFR_output_SR.mat

%%
%%--------------------------------Stats reported in Main text with no figure-------------------------------------%%

%% No difference in running speed between correct and error trials
clear
code_behavior_speed
% Relevant MAT files: Stats_Behavior_speed.mat

%% No difference in time spent at reward location or in the beginning of the rest box
clear
code_behavior_measure_time
% Relevant MAT files: Stats_Behavior_time.mat

%% No difference in duration of sequences between correct and error trials
clear
Stats_gammaTFR_eachseq_forMultiCompare_tspan
% Relevant MAT files: Stats_SeqDuration.mat

%% Difference in decoded path length of sequences between correct and error trials
clear
Stats_gammaTFR_eachseq_forMultiCompare_xspan
% Relevant MAT files: Stats_SeqXspan.mat

%%
%%---------------------------------New figures added after--------------%%

%% Separate replay events into before and after test trials
clear
parentfd = ReturnParentPath;
directories_allData_v1;
cellList = 'TTList_dCA1_pyr_ds.txt';
spktime_fd = 'SpikeTime_ds';
% Plot analysis result again in the same folder as other figures
Outdir2 = [parentfd,'\GroupData Figures'];
RunBayesianDecoder_ripple_BeforeAfterTest(path,CSClist_CA1,Outdir2,cellList,spktime_fd,trackdata,true);

%% Use bandpass filter instead of wavelet
clear
doall_BandpassFilter
doall_get_gammaTFR_v6_SR_bp
doall_get_gammaTFR_eachseq_EH_v2_bp
Stats_gammaTFR_eachseq_forMultiCompare_bp

%% Normalize duration of replay events
clear
parentfd = ReturnParentPath;
directories_allData_v1;
cellList = 'TTList_dCA1_pyr_ds.txt';
spktime_fd = 'SpikeTime_ds';
RunBayesianDecoder_ripple_NormTime(path,CSClist_CA1,Outdir2,cellList,spktime_fd,trackdata,true);