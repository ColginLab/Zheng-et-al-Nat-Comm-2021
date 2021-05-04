%% Suplementary Figure 1 - DONE
% Relevant figures: PlaceFieldsDist_alignedToReward\PlaceFieldDis_after,
%                   CrossValidation

function SupFigure1
clearvars
close all
clc
FiguresDir = pwd;
ind = strfind(FiguresDir,'MATLAB');
CoDatCirDir = [FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track'];
addpath(CoDatCirDir)
BayDecDir = [FiguresDir(1:ind+5) '\MainFunctions\Code & Data Circular Track\Analysis_EH\BayesianDecoding_ripple'];
addpath(BayDecDir)

%%
%%parentfd = ReturnParentPath;
directories_allData_v1;
cellList = 'TTList_dCA1_pyr.txt';
spktime_fd = 'SpikeTime';
%Outdir = strcat(parentfd,'\Analysis_EH\');
Outdir = 'E:\ColginLab\Figures\SupFigure1';

Preprocess_done = true;
%Preprocess_done = false;
if ~Preprocess_done
    for ff = 1:length(pathRats)
        MatchSpikeTimeToEEGTime_v3(pathRats{ff},spktime_fd,cellList);           % IMPORTANT: FUNCTION NOT FOUND
    end
end

%if ~isdir(Outdir)
if verLessThan('matlab','8.4')
    error('Please use Matlab version 2014b or above')
else
    RunPlaceCellPropertiesAcrossLap_v6(pathRats(1:end-1),trackdata,Outdir,cellList,spktime_fd);
end
%end
% Plot analysis result in the same folder as other figures
%Outdir2 = strcat(parentfd,'\GroupData Figures\');
Outdir2 = Outdir;
PlaceFieldsDist(trackdata,Outdir,Outdir2);

% Do cross validation
%InputFile = strcat(parentfd,'\Analysis_EH_ds\ConfusionMatrix\ConfusionMat.mat');
InputFile = 'E:\ColginLab\Data Analysis\Analysis_EH_ds\ConfusionMatrix\ConfusionMat.mat';
if exist(InputFile,'file')~=2
    MasterScript_CA1Replay % Run Ernie's script
end
%Outdir = strcat(parentfd,'\GroupData Figures\');
PlotCrossValidation(InputFile,Outdir);

%%
cd(FiguresDir)
fprintf('Done! \n')
end

