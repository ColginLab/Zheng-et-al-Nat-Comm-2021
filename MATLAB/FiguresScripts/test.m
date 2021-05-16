%%
close all
clear all
clc

DataDir = 'E:\ColginLab\Data\Circular Track';
AnalysisDir = 'E:\ColginLab\Data Analysis\GroupData';
DownSampleFlag = false;
file_analysis_name_ext = fullfile(AnalysisDir, 'group_gammaTFR_eachseq_20190529.mat');  % 4 rats % remove jumping-out points, ind_approach = 3;
FiguresDir = 'E:\ColginLab\github-repo\Figures\Figure5';

Figure5(DataDir, AnalysisDir, DownSampleFlag, file_analysis_name_ext, FiguresDir)

%Figures3_4(DataDir, AnalysisDir, file_analysis_name, FiguresDir)