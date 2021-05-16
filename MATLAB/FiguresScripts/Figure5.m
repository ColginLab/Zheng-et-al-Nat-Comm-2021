function Figure5(DataDir, AnalysisDir, DownSampleFlag, file_analysis_name_ext, FiguresDir)
%FIGURE5 Function that plots Figure 5
%   ARGUMENTS:
%       DataDir:                Raw data directory
%       AnalysisDir:            Directory to store analysis results
%       DownSampleFlag:         Boolean flag for downsampling cells setting
%       file_analysis_name_ext: Full path of meta-analysis file
%       FiguresDir:             Directory to save figure(s). eps and fig format
%   RELEVANT FIGURES:
%       Fig5
%   NOTE: These scripts take a while to run

close all
clc
FuncDir = pwd;
%%
ind = strfind(FuncDir,'MATLAB');
CodeDatCirDir = [FuncDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'];
addpath(CodeDatCirDir)
BayesDecodDir = [FuncDir(1:ind+5) '\Dependencies\BayesianDecodingRipple'];
addpath(BayesDecodDir)
circpackDir = [FuncDir(1:ind+5) '\Dependencies\circpackage'];
addpath(circpackDir)
GenDir = [FuncDir(1:ind+5) '\GeneralFunctions'];
addpath(GenDir)
cd([FuncDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'])
%%
doall_get_gammaTFR_eachseq_EH(file_analysis_name_ext, DataDir);
group_spkphase_approachReward_v3(DataDir, AnalysisDir, DownSampleFlag, file_analysis_name_ext, FiguresDir)
%%
cd(FuncDir)
fprintf('Done! \n')
end