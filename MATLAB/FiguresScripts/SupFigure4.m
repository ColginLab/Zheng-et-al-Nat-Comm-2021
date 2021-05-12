%% Suplementary Figure 4 - PROBLEM: DOES NOT REPRODUCE THE SAME FIGURE IN THE MANUSCRIPT
% Relevant figures: SeqSlopeNeg_average_samp, 
%                   SeqSlopeNeg_average_test
% Relevant MAT files: Stats_SeqSlopeNeg.mat

function SupFigure2
clearvars
close all
clc
FiguresDir = pwd;

%%
ind = strfind(FiguresDir,'MATLAB');
CodeDatCirDir = [FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'];
addpath(CodeDatCirDir)
BayesDecodDir = [FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track\Analysis_EH\BayesianDecoding_ripple'];
addpath(BayesDecodDir)
GenDir = [FiguresDir(1:ind+5) '\GeneralFunctions'];
addpath(GenDir)

cd([FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'])

%%
%doall_get_gammaTFR_eachseq_EH_v2
Stats_gammaTFR_eachseq_forMultiCompare_rev

cd(FiguresDir)
fprintf('Done! \n')

end