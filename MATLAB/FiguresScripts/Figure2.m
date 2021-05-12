function Figure2(DataDir, AnalysisDir, FiguresDir)
%FIGURE2 Function that plots Figure 2
%   ARGUMENTS:
%       DataDir:        Raw data directory
%       AnalysisDir:    Directory to store analysis results
%       FiguresDir:     Directory to save figure(s). eps and fig format
close all
clc
FuncDir = pwd;
%%
genFDir = '..\GeneralFunctions';
addpath(genFDir)
ind = strfind(FuncDir,'MATLAB');
CoDatCirDir = [FuncDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'];
addpath(CoDatCirDir)
cd([FuncDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'])
%%
doall_group_Pxn(DataDir, AnalysisDir)
close all
Stats_PxnSum_CI_Fig2(AnalysisDir, FiguresDir)
%%
cd(FuncDir)
fprintf('Done! \n')
end