%% Figure 5: Theta phase locking - DONE
% Relevant figures: SpkPhaseLock_AcrossTrack
% Relevant MAT files: group_spkphase_v2.mat
% This script takes a while to run

function Figure5
clearvars
close all
clc

FiguresDir = pwd;

%%
ind = strfind(FiguresDir,'MATLAB');
CodeDatCirDir = [FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'];
addpath(CodeDatCirDir)
BayesDecodDir = [FiguresDir(1:ind+5) '\Dependencies\BayesianDecodingRipple'];
addpath(BayesDecodDir)
circpackDir = '..\Dependencies\circpackage';
addpath(circpackDir)
genFDir = '..\GeneralFunctions';
addpath(genFDir)

%%
%doall_get_gammaTFR_eachseq_EH
group_spkphase_approachReward_v3(false)

%%
cd(FiguresDir)
fprintf('Done! \n')
end