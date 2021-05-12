%% Figure 6: No difference in replay fidelity between correct and error trials
% Relevant figures: BayesianDecoding_ripple\R2distCorrectVSIncorrect, 
%                   BayesianDecoding_ripple\BayesianDecodingExamples

function Figure6

clc
close all
clearvars
FiguresDir = pwd;

%%
ind = strfind(FiguresDir,'MATLAB');
CodeDatCirDir = [FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'];
addpath(CodeDatCirDir)
%BayesDecodDir = [FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track\Analysis_EH\BayesianDecoding_ripple'];
BayesDecodDir = [FiguresDir(1:ind+5) '\Dependencies\BayesianDecodingRipple'];
addpath(BayesDecodDir)

directories_allData_v1;
cellList = 'TTList_dCA1_pyr_ds.txt';
spktime_fd = 'SpikeTime_ds';
% Analysis output directory
%Outdir = strcat(parentfd,'\Analysis_EH_ds\');
Outdir = 'E:\ColginLab\Data Analysis\Analysis_EH_ds';
%if ~isdir(Outdir)
   %MasterScript_CA1Replay % Run Ernie's script
%end
% Plot analysis result again in the same folder as other figures
%Outdir2 = [parentfd,'\GroupData Figures'];
%Outdir2 = 'E:\ColginLab\Figures\Figure6';

FigDir{1} = 'E:\ColginLab\Figures\Figure6';
FigDir{2} = 'E:\ColginLab\Figures\Figure7';
FigDir{3} = 'E:\ColginLab\Figures\SupFigure5';
FigDir{4} = 'E:\ColginLab\Figures\SupFigure9';

RunBayesianDecoder_ripple(pathRats,CSClist_CA1,Outdir,cellList,spktime_fd,trackdata,false,FigDir);
RunBayesianDecoder_ripple_NormTime(pathRats,CSClist_CA1,Outdir,cellList,spktime_fd,trackdata,true);
% Relevant figures: BayesianDecoding_ripple\R2distCorrectVSIncorrect, 
%                   BayesianDecoding_ripple\BayesianDecodingExamples
end