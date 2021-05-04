function Figure8
clc
close all
clear all
FiguresDir = pwd;

%%
ind = strfind(FiguresDir,'MATLAB');
CodeDatCirDir = [FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'];
addpath(CodeDatCirDir)
BayesDecodDir = [FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track\Analysis_EH\BayesianDecoding_ripple'];
addpath(BayesDecodDir)
genFDir = '..\GeneralFunctions';
addpath(genFDir)
circpackDir = [FiguresDir(1:ind+5) '\Dependencies\circpackage'];
addpath(circpackDir)
%%

%group_gammaphase_eachseq        % WORKS
Stats_spkphaseLock
disp('done')

end