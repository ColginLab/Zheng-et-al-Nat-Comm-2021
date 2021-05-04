%% Figure 1: Behavioral performance - DONE
% Relevant figures: Behavior_barGraph, 
%                   Behavior_learningCurve
function Figure1
clearvars
close all
clc
FiguresDir = pwd;
matbugsDir = '..\Dependencies\matbugs';
addpath(matbugsDir)
genFDir = '..\GeneralFunctions';
addpath(genFDir)

%%
ind = strfind(FiguresDir,'MATLAB');
cd([FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'])
group_behavior_statespace_all('BUGS')                                        
cd([FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'])
code_behavior_v1   

%%
cd(FiguresDir)
fprintf('Done! \n')
end