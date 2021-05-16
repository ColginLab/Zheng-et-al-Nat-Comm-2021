function SupFigure5(DataDir, AnalysisDir, file_analysis_name_ext, file_analysis_out_ext, FiguresDir)
%SUPFIGURE5 Function that plots Supplementary Figure 5
%   ARGUMENTS:
%       DataDir:                Raw data directory
%       AnalysisDir:            Directory to store analysis results
%       file_analysis_name_ext: Full path of meta-analysis file
%       file_analysis_out_ext:  Full path of meta-analysis file (gamma)
%       FiguresDir:             Directory to save figure(s). eps and fig format
%   RELEVANT FIGURES:
%       SupFigure5a, SupFigure5b, SupFigure_cf
%   NOTE: These scripts take a while to run
clc
close all
FuncDir = pwd;
%%
ind = strfind(FuncDir,'MATLAB');
CodeDatCirDir = [FuncDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'];
addpath(CodeDatCirDir)
BayesDecodDir = [FuncDir(1:ind+5) '\Dependencies\BayesianDecodingRipple'];
addpath(BayesDecodDir)
GenDir = [FuncDir(1:ind+5) '\GeneralFunctions'];
addpath(GenDir)
circpackDir = [FuncDir(1:ind+5) '\Dependencies\circpackage'];
addpath(circpackDir)
cd([FuncDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'])
%%
group_gammaphase_eachseq(DataDir, AnalysisDir, file_analysis_name_ext, file_analysis_out_ext, FiguresDir);
Stats_spkphaseLock(file_analysis_out_ext, AnalysisDir, FiguresDir);
%%
cd(FuncDir)
fprintf('Done! \n')
end