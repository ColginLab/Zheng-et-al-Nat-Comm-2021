%% us State-Space model to calculate the learning curve
% edited on 4/24/2017

%clearvars
%close all
%clc
function group_behavior_statespace_all(alg)
parentfd = fileparts(mfilename('fullpath'));
indMATLAB = strfind(parentfd,'MATLAB');
parentfd = parentfd(1:indMATLAB-2);

%alg = 'BUGS';                 % Algorithm to estimate learning curve, options: BUGS, EM

if strcmpi(alg, 'BUGS')
    %dirfolder_alg = [parentfd '\Behavior\State-space Model\Learning\IndividualWinBUGS'];
    %dirfolder_alg = [parentfd '\Behavior'];
    
    %dirfolder_alg = [parentfd,'\State-space Model\Learning Analysis\IndividualAnalysisWinBUGS'];
    dirfolder_alg = [parentfd,'\Behavior\State-space Model\Learning\IndividualAnalysisWinBUGS'];   
    %dirfolder_figs = [parentfd '\Behavior\State-space Model\Figures\WinBUGS'];
    dirfolder_figs = 'E:\ColginLab\Figures\Figure1';
    dirfolder_output = 'E:\ColginLab\Behavior\State-space Model\Data Analysis\WinBUGS';
    %dirfolder_output = [parentfd '\Behavior\State-space Model\WinBUGS'];
    
    %dirfolder_figs = [parentfd,'\MainFunctions\Code & Data Circular Track\GroupData Figures\WinBUGS'];
    %dirfolder_output = [parentfd,'\MainFunctions\Code & Data Circular Track\GroupData\WinBUGS'];
elseif strcmpi(alg, 'EM')           % TODO
    dirfolder_alg = [parentfd,'\State-space Model\Learning Analysis\IndividualAnalysisEM'];
    dirfolder_figs = [parentfd,'\MainFunctions\Code & Data Circular Track\GroupData Figures\EM'];
    dirfolder_output = [parentfd,'\MainFunctions\Code & Data Circular Track\GroupData\EM'];
    MaxResponse = 1;
end
file_output = 'group_behavior_statespace_all.mat';

if ~isdir(dirfolder_figs)
    mkdir(dirfolder_figs)
end
if ~isdir(dirfolder_output)
    mkdir(dirfolder_output)
end


%dirfolder_BUGS = [parentfd,'\MATLAB\State-space Model\Learning Analysis\IndividualAnalysisWinBUGS\'];
%dirfolder_EM = [parentfd,'\MATLAB\State-space Model\Learning Analysis\IndividualAnalysisEM\'];
%dirfolder_figs = [parentfd,'\MATLAB\MainFunctions\Code & Data Circular Track\GroupData Figures'];

%file_input = [parentfd(1:end-7) '\Behavior\data_behavior_statespace_all.mat'];
file_input = [parentfd '\Behavior\data_behavior_statespace_all.mat'];
load(file_input)

%% ========= use all trials across days to make a learning curve =========
%% learning curve for stop around reward location vs. non-stop
% correct = stop within +-3 around the correct location
% incorrect = non-stop


% % calculate the learning curve for test trials only
% ind = find(data_behavior(:,3) == 3);
% data_behavior_test = data_behavior(ind,:);
% ind_crt = find(~isnan(data_behavior_test(:,4)));
% ind_incrt = find(isnan(data_behavior_test(:,4)));
% Responses = [];
% Responses(1,ind_crt) = 1;
% Responses(1,ind_incrt) = 0;
% 
% % using BUGS method
% cd(dirfolder_BUGS)
% MaxResponse = ones(size(Responses));
% BackgroundProb = 5/19;
% runanalysis(Responses,MaxResponse, BackgroundProb);
% plotresults
% 
% % using EM method
% cd(dirfolder_EM)
% MaxResponse = 1;
% BackgroundProb = 5/19;
% runanalysis(Responses,MaxResponse, BackgroundProb);
% plotresults


%% ========= use trials in each session to make a learning curve =========

%% calculate the learning curve 

ind = find(data_behavior(:,3) == 3);
data_behavior_test = data_behavior(ind,:);
ns_all = max(data_behavior_test(:,2));
pdata_behavior0 = nan(20,4,1000); % max 20 trials X 4 col X ns sessions
Ind_session = [];
for ns = 1:ns_all
    ind = find(data_behavior(:,2) == ns);
    data_behavior_ns = data_behavior(ind,:);
    
    % Case1: 
    % ** correct = stop within +-3 around the correct location
    % ** incorrect = non-stop
%     ind_crt = find(~isnan(data_behavior_ns(:,4)));
%     ind_incrt = find(isnan(data_behavior_ns(:,4)));
%     BackgroundProb = 5/19;
    
    % Case2: 
    % ** correct = stop at the correct location
    % ** incorrect = non-stop and/or stop at the incorrect location
    ind_crt = find(data_behavior_ns(:,4) == 1);
    ind_incrt = find(data_behavior_ns(:,4) ~= 1);
    BackgroundProb = 1/5;
    
    if ~isempty(ind_crt) & ~isempty(ind_incrt)
        Responses = [];
        Responses(1,ind_crt) = 1;
        Responses(1,ind_incrt) = 0;
        
        % using BUGS method
        %cd(dirfolder_BUGS)
        %
        
        if strcmpi(alg, 'BUGS')
            MaxResponse = ones(size(Responses));
        end
        
        % using EM method
        cd(dirfolder_alg)
        %MaxResponse = 1;
        
        
        runanalysis(Responses,MaxResponse, BackgroundProb);
        plotresults
        pause(0.1)
        [m,n] = size(pdata);
        pdata_behavior0(1:m,1:n,ns) = pdata;
        Ind_session = [Ind_session;ns];
    end
end
close all

cd(dirfolder_output)
save(file_output);

%% choose data from different rats to plot
Rat_plot = [139,150,149,148];  %% NEED TO CHANGE

[~,ia,~] = unique(data_behavior(:,2));
data_behavior_unique = data_behavior(ia,1:2);
Lia = ismember(data_behavior_unique(:,1),Rat_plot);
data_behavior_unique = data_behavior_unique(Lia,:);
Lia = ismember(Ind_session,data_behavior_unique(:,2));

pdata_behavior = pdata_behavior0(:,:,Ind_session);
pdata_behavior = pdata_behavior(:,:,Lia);
pdata_behavior_mean = nanmean(pdata_behavior,3);
x_trial = pdata_behavior_mean(:,1);
ffa = figure('Units','normalized','Position',[0 0.2 .7 0.5]);
color_CI = [0 0 0]./255;
color_test = [240,199,148]./255;
color_posttest = [231,229,178]./255;
hold on;
fill([0.5,8.5, fliplr([0.5,8.5])],[0,0,1.2,1.2],color_test,'linestyle','none');
fill([8.5,14.5, fliplr([8.5,14.5])],[0,0,1.2,1.2],color_posttest,'linestyle','none');
plot(pdata_behavior_mean(:,3),'s-','color', 'k','MarkerFaceColor','k','MarkerEdgeColor', 'k','linewidth',2); 
plot(pdata_behavior_mean(:,2),'color', color_CI,'linewidth',2); 
plot(pdata_behavior_mean(:,4),'color', color_CI,'linewidth',2); 
line([1,14],[BackgroundProb BackgroundProb],'color','b')
text(4.5,1.1,'Test trials','HorizontalAlignment','center','FontSize',20)
text(10.5,1.1,'Post-Test trials','HorizontalAlignment','center','FontSize',20)
xlim([0.5, 14.5])
ylim([0,1.15])
set(gca,'XTick',1:length(x_trial));
set(gca,'XTickLabel',{'t1','t2','t3','t4','t5','t6','t7','t8','pt1','pt2','pt3','pt4','pt5','pt6'});
set(gca,'YTick',0:0.2:1);
set(gca,'fontsize',22);
xlabel('Trial Number');
ylabel('Prob (Correct Response)')
hold off
title('Learning curve across all recording days')

saveas(ffa,[dirfolder_figs,'\Behavior_learningCurve'],'epsc')
saveas(ffa,[dirfolder_figs,'\Behavior_learningCurve'],'fig')
close all
end