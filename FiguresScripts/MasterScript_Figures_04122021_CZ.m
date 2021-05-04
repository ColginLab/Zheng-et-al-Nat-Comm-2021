%% Master script to produce figures
% Figures are stored in the following directory: 
% ~\MATLAB\MainFunctions\Code & Data Circular Track\GroupData Figures
% Figures made for publication are stored in the folowing directory:
% ~\MATLAB\MainFunctions\Code & Data Circular Track\Figures_manuscript
% Run all figures using Matlab2014a except for Figure 6, use Matlab2018b


%%-----------------------------------Main figures------------------------------------------%%

%% Figure 1: Behavioral performance
clear
group_behavior_statespace_all
code_behavior_v1
% Relevant figures: Behavior_barGraph, 
%                   Behavior_learningCurve

%% Figure 2: Predictive firing in place cell ensembles developed with learning
% 
doall_group_Pxn
Stats_PxnSum_CI_Fig2

%% Figure 3: Place cell sequences exhibited steeper slopes in correct test trials than in error test trials
clear
doall_get_gammaTFR_eachseq_EH
Stats_gammaTFR_eachseq_forMultiCompare
BayesianDecodingExample
% Relevant figures: SeqSlope_test, 
%                   SeqExample_test_Rat149_2017-08-08-CT-2
% Relevant MAT files: Stats_SeqSlope.mat

%% Figure 4: Place cell sequences exhibited steeper slopes in correct sample trials than in error sample trials
% figures are already produced by the function
%   "Stats_gammaTFR_eachseq_forMultiCompare" and "BayesianDecodingExample"
% Relevant figures: SeqSlope_samp,
%                   SeqExample_sample_Rat149_2017-08-08-CT-2

%% Figure 5: Place cells active at the beginning of the trajectory toward the reward fired at earlier theta phases... 
%            in correct trials than in error trials
clear
doall_get_gammaTFR_eachseq_EH
group_spkphase_approachReward_v2(false)
% Relevant figures: SpkPhaseLock_AcrossTrack
% Relevant MAT files: group_spkphase_v2.mat

%% Figure 6: Replay quality in the rest box was similar between correct and error trials
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
RunBayesianDecoder_ripple_NormTime(path,CSClist_CA1,Outdir2,cellList,spktime_fd,trackdata,true);
% Relevant figures: BayesianDecoding_ripple\R2distCorrectVSIncorrect, 
%                   BayesianDecoding_ripple\BayesianDecodingExamples

%% Figure 7: A bias for replay events to terminate at the correct goal location emerged after learning during correct trials
% figures are already produced by the function "RunBayesianDecoder_ripple_NormTime"
% Relevant figures: BayesianDecoding_ripple_NormTime_v2\TerminationBias_goalPos_Error, 
%                   BayesianDecoding_ripple_NormTime_v2\TerminationBias_goalPos_Correct

%% Figure 8: Locations represented at the end of replay sequences plotted against animals?actual stop locations.
% figures are already produced by the function "RunBayesianDecoder_ripple"
% Relevant figures: BayesianDecoding_ripple\StopVSPredictLocCorr

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

%% Supplementary figure 2: Examples of posterior probability distributions from test trials with different types of errors
% figures are already produced by the function "BayesianDecodingExample"
% Relevant figures: SeqExample_test_Rat149_2017-08-08-CT-1, 
%                   SeqExample_test_Rat148_2017-12-02-CT

%% Supplementary figure 3: Running speeds during the approach to the stop location
clear
Stats_speedAcrossPos
% Relevant figures: SpeedThetaAcrossPos_test, 
%                   SpeedThetaAcrossPos_samp

%% Supplementary figure 4: Distances (x-span? and durations (t-span? of sequences during correct and error trials
% 
clear
Stats_gammaTFR_eachseq_forMultiCompare_pos_neg
% need to comment/uncomment line259-269 to do x-span/x-span diff/t-span for supp fig4ad/be/cf respectively

%% Supplementary Figure 5: Gamma phase modulation of place cell spikes was similar between different trial types
clear
group_gammaphase_eachseq
Stats_spkphaseLock
% Relevant figures: SlowG_spkPhaseLock, 
%                   FastG_spkPhaseLock,
%                   Stats_spkphaseLock
% Relevant MAT files: Stats_spkphaseLock.mat

%% Supplementary figure 6: No difference in gamma power was observed between trial types as rats approached a stop location
clear
doall_PFR
% one figure is already produced by the function "Stats_gammaTFR_eachseq_forMultiCompare"
% Relevant figures: SampTest_CrtErr_CI_SR, 
%                   SampTest_CrtErr_SR, 
%                   SeqSlopePower
% Relevant MAT files: PFR_output_SR.mat

%% Supplementary figure 7: Place cell sequences with negative slopes detected during correct and error trials
clear
doall_get_gammaTFR_eachseq_EH_v2
Stats_gammaTFR_eachseq_forMultiCompare_rev
% Relevant figures: SeqSlopeNeg_average_samp, 
%                   SeqSlopeNeg_average_test
% Relevant MAT files: Stats_SeqSlopeNeg.mat

%% Supplementary Figure 8: Forward and reverse replay fidelity was similar between rest periods from correct and error trials
% figures are already produced by the function "RunBayesianDecoder_ripple"
% Relevant figures: BayesianDecoding_ripple\R2distCorrectVSIncorrect_forVSrev

%% Supplementary Figure 9: A bias for replay events to terminate at the correct reward location occurs in rest periods before... 
%                          and after the test phase during correct trials but not error trials
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
RunBayesianDecoder_ripple_BeforeAfterTest(path,CSClist_CA1,Outdir2,cellList,spktime_fd,trackdata,true);
% Relevant figures: BayesianDecoding_ripple_BeforeAfterTest\TerminationBias_goalPos_Error, 
%                   BayesianDecoding_ripple_BeforeAfterTest\TerminationBias_goalPos_Correct

%% Supplementary Figure 10: A significant bias for replay events to terminate at the correct reward location did not occur in...
%                           post-test trials and pre-running trials
% figures are already produced by the function "RunBayesianDecoder_ripple"
% Relevant figures: BayesianDecoding_ripple\TerminationBias_goalPos_Error, 
%                   BayesianDecoding_ripple\TerminationBias_goalPos_Correct
%% Supplementary Fig. 11:  No bias for replay events to start at the correct reward location was observed
% figures are already produced by the function "RunBayesianDecoder_ripple_NormTime"
% Relevant figures: BayesianDecoding_ripple_NormTime_v2\InitiationBias_goalPos_Correct
%                   BayesianDecoding_ripple_NormTime_v2\InitiationBias_goalPos_Error
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

