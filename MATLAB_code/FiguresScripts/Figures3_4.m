%% Figures 3 and 4 - DONE
% For figure 3: Difference in the slope of place cell sequence between correct and error test trials
% Relevant figures: SeqSlope_test, 
%                   SeqExample_test_Rat149_2017-08-08-CT-2
% Relevant MAT files: Stats_SeqSlope.mat
% For Figure 4: figures are already produced by the function
% "Stats_gammaTFR_eachseq_forMultiCompare" and "BayesianDecodingExample"
% Relevant figures: SeqSlope_samp,
%                   SeqExample_sample_Rat149_2017-08-08-CT-2
% NOTE: These scripts take a while to run

function Figures3_4
clearvars
clc
close all
FiguresDir = pwd;

%%
ind = strfind(FiguresDir,'MATLAB');
CodeDatCirDir = [FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'];
addpath(CodeDatCirDir)
BayesDecodDir = [FiguresDir(1:ind+5) '\Dependencies\BayesianDecodingRipple'];
addpath(BayesDecodDir)
GenDir = [FiguresDir(1:ind+5) '\GeneralFunctions'];
addpath(GenDir)

cd([FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'])

%%
doall_get_gammaTFR_eachseq_EH
Stats_gammaTFR_eachseq_forMultiCompare
BayesianDecodingExample

cd(FiguresDir)
fprintf('Done! \n')

end