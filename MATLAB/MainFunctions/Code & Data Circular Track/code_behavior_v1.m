% edited on 4/24/2017
% Behavior data
% check the spatial preference and memory interference in Pre-running session
% calculate the stops count for correct and error trials in Test and Post-test trials

%% data in pre-running session
% check the spatial preference and memory interference in Pre-running session

function code_behavior_v1

pwdir = pwd;

directories_allData_v1
%parentfd = fileparts(mfilename('fullpath'));
%dirfolder = [parentfd,'\GroupData\'];
dirfolder = 'E:\ColginLab\Data Analysis\GroupData';
file_input = 'Data_angle_ontrack.mat';
file_output = 'data_behavior_v1.mat';

id_rewardloc = 6:15;  % this is very important! Keep it the same for each rat in the experiment

time_rewardloc_all = []; % time spent at each reward location, save data for each session in columns
perc_rewardloc_all = []; % percentage of time spent at each reward location, save data for each session in columns
vel_all = []; % running speed
vel_ang_all = []; % angular velocity
perc_rewardloc_3 = []; % percentage of time spent at each reward location categories
Ratid = [];

for ns = 1:isession
    path_ns = pathRats{ns};
    cd(path_ns);
    
    % load angle data in pre-running trials
    load(file_input)
    nlap = size(data_angle{1,1},1);
    data0 = [];
    for nl = 1:nlap
        data0 = [data0;data_angle{1,1}{nl}];
    end
    
    % read reward location information
    trackdata_ns = trackdata{ns};
    load(trackdata_ns,'Ang_RewardLoc_ontrack','N_loc','Ind_rewardloc','vfs')
    
    % left and right border of each reward locations
    para=1/2;
    Ang_RewardLoc_ontrack_zone=nan(N_loc,2);
    for i=2:N_loc
        Ang_RewardLoc_ontrack_zone(i,1)=Ang_RewardLoc_ontrack(i)-(Ang_RewardLoc_ontrack(i)-Ang_RewardLoc_ontrack(i-1))*para;
    end
    Ang_RewardLoc_ontrack_zone(1,1)=Ang_RewardLoc_ontrack(1);
    for i=1:N_loc-1
        Ang_RewardLoc_ontrack_zone(i,2)=Ang_RewardLoc_ontrack(i)+(Ang_RewardLoc_ontrack(i+1)-Ang_RewardLoc_ontrack(i))*para;
    end
    Ang_RewardLoc_ontrack_zone(N_loc,2)=Ang_RewardLoc_ontrack(N_loc);
    
    % 1. Calculate: spatial preference (time spent) at each reward locations
    time_rewardloc_ns = nan(1,N_loc);
    perc_rewardloc_ns = nan(1,N_loc);
    vel_ns = nan(1,N_loc);
    vel_ang_ns = nan(1,N_loc);
    for iloc = id_rewardloc
        ind = find(data0(:,2) >= Ang_RewardLoc_ontrack_zone(iloc,1) & ...
            data0(:,2) <= Ang_RewardLoc_ontrack_zone(iloc,2));
        time_rewardloc_ns(iloc) = length(ind)./vfs;
        perc_rewardloc_ns(iloc) = length(ind)./size(data0,1);
        vel_ns(iloc) = nanmean(abs(data0(ind,3)));
        vel_ang_ns(iloc) = nanmean(abs(data0(ind,4)));
    end
    time_rewardloc_all(ns,:) = time_rewardloc_ns;
    perc_rewardloc_ns = perc_rewardloc_ns./nansum(perc_rewardloc_ns); % relative value for this percentage
    perc_rewardloc_all(ns,:) = perc_rewardloc_ns;
    vel_all(ns,:) = vel_ns;
    vel_ang_all(ns,:) = vel_ang_ns;
    
    % 2. Calculate time spent at 3 categories of reward locations:
    %(1) current reward location (in the following sample-test trials)
    
    id_loc_now = Ind_rewardloc;
    if id_loc_now ~= rewardcurrent(ns)
        warning('The current reward location is not matching up with notes.')
    end
    time_rewardloc_now = time_rewardloc_ns(id_loc_now);
    perc_rewardloc_now = perc_rewardloc_ns(id_loc_now);
    %(2) last reward location from the preceding session (memory interference)
    id_loc_last = rewardlast(ns);
    time_rewardloc_last = time_rewardloc_ns(id_loc_last);
    perc_rewardloc_last = perc_rewardloc_ns(id_loc_last);
    
    %(3) other locations
    id_loc_other = setdiff(id_rewardloc,id_loc_now);
    id_loc_other = setdiff(id_loc_other,id_loc_last);
    time_rewardloc_other = nanmean(time_rewardloc_ns(id_loc_other));
    perc_rewardloc_other = nanmean(perc_rewardloc_ns(id_loc_other));
    
    time_rewardloc_3(ns,:) = [time_rewardloc_last,time_rewardloc_now,time_rewardloc_other];
    perc_rewardloc_3(ns,:) = [perc_rewardloc_last,perc_rewardloc_now,perc_rewardloc_other];
    Ratid(ns,1) = Ind_Rat(ns);
end

cd(dirfolder)

save(file_output,'Ratid','perc_rewardloc_all','time_rewardloc_all','vel_all','vel_ang_all','perc_rewardloc_3','time_rewardloc_3','id_rewardloc');

%% Plots in pre-running session
file_input = 'data_behavior_v1.mat';
%%load([dirfolder,file_input])

Plot_rat = [139,150,149,148];
ind = find(ismember(Ratid,Plot_rat) > 0);
perc_rewardloc_all = perc_rewardloc_all(ind,:);
vel_all = vel_all(ind,:);
vel_ang_all = vel_ang_all(ind,:);
perc_rewardloc_3 = perc_rewardloc_3(ind,:);

% Plot the percentage of time spent at each location, and corresponding
% running speed
num_rewardloc = length(id_rewardloc);
perc_rewardloc = perc_rewardloc_all(:,id_rewardloc);
perc_rewardloc_mean = 100*mean(perc_rewardloc);
perc_rewardloc_sem = 100*std(perc_rewardloc)./sqrt(size(perc_rewardloc,1));
vel0 = vel_all(:,id_rewardloc);
vel0_mean = mean(vel0);
vel0_sem = std(vel0)./sqrt(size(vel0,1));
vel_ang0 = vel_ang_all(:,id_rewardloc);
vel_ang0_mean = mean(vel_ang0);
vel_ang0_sem = std(vel_ang0)./sqrt(size(vel_ang0,1));

ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
bar_color = [100,180,220]./255;
hold on
bar(1:num_rewardloc,perc_rewardloc_mean,'FaceColor',bar_color,'EdgeColor', bar_color);
h11=errorbar(1:num_rewardloc,perc_rewardloc_mean,perc_rewardloc_sem,'color',[0 0 0]./255,'LineStyle','none','LineWidth',2);
%errorbar_tick(h11, 0);                 BUG
hold off
set(gca, 'FontSize',20);
xlim([0, num_rewardloc+1])
ylim([0,20])
set(gca, 'XTick',1:num_rewardloc);
xlabel('Reward location #');
ylabel('Percentage of time spent at each location (%)');
title('Pre-running Probe')

ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
bar_color = [100,180,220]./255;
hold on
h11=errorbar(1:num_rewardloc,vel0_mean,vel0_sem,'s','color',[0 0 0]./255,'LineStyle','none','LineWidth',2);
%errorbar_tick(h11, 0);
hold off
set(gca, 'FontSize',20);
xlim([0, num_rewardloc+1])
ylim([0,80])
set(gca, 'XTick',1:num_rewardloc);
xlabel('Reward location #');
ylabel('Running speed (cm/s)');
title('Pre-running Probe')

% Plot the time spent at each location, and corresponding running speed
num_rewardloc = length(id_rewardloc);
time_rewardloc = time_rewardloc_all(:,id_rewardloc);
time_rewardloc_mean = mean(time_rewardloc);
time_rewardloc_sem = std(time_rewardloc)./sqrt(size(time_rewardloc,1));
vel0 = vel_all(:,id_rewardloc);
vel0_mean = mean(vel0);
vel0_sem = std(vel0)./sqrt(size(vel0,1));
vel_ang0 = vel_ang_all(:,id_rewardloc);
vel_ang0_mean = mean(vel_ang0);
vel_ang0_sem = std(vel_ang0)./sqrt(size(vel_ang0,1));

ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
bar_color = [100,180,220]./255;
hold on
bar(1:num_rewardloc,time_rewardloc_mean,'FaceColor',bar_color,'EdgeColor', bar_color);
h11=errorbar(1:num_rewardloc,time_rewardloc_mean,time_rewardloc_sem,'color',[0 0 0]./255,'LineStyle','none','LineWidth',2);
%errorbar_tick(h11, 0);
hold off
set(gca, 'FontSize',20);
xlim([0, num_rewardloc+1])
ylim([0,6])
set(gca, 'XTick',1:num_rewardloc);
xlabel('Reward location #');
ylabel('Time spent at each location (%)');
title('Pre-running Probe')

ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
bar_color = [100,180,220]./255;
hold on
h11=errorbar(1:num_rewardloc,vel0_mean,vel0_sem,'s','color',[0 0 0]./255,'LineStyle','none','LineWidth',2);
%errorbar_tick(h11, 0);
hold off
set(gca, 'FontSize',20);
xlim([0, num_rewardloc+1])
ylim([0,80])
set(gca, 'XTick',1:num_rewardloc);
xlabel('Reward location #');
ylabel('Running speed (cm/s)');
title('Pre-running Probe')

% Plot the percentage of time spent at (1) last reward location; (2)
% current reward location; (3) Other locations
ind = find(~isnan(sum(perc_rewardloc_3,2)));
perc_rewardloc_3 = perc_rewardloc_3(ind,:);
perc_rewardloc_3_mean = 100*mean(perc_rewardloc_3);
perc_rewardloc_3_sem = 100*std(perc_rewardloc_3)./sqrt(size(perc_rewardloc_3,1));

ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
bar_color = [100,180,220]./255;
hold on
bar(1:3,perc_rewardloc_3_mean,'FaceColor',bar_color,'EdgeColor', bar_color);
h11=errorbar(1:3,perc_rewardloc_3_mean,perc_rewardloc_3_sem,'color',[0 0 0]./255,'LineStyle','none','LineWidth',2);
%errorbar_tick(h11, 0);
hold off
set(gca, 'FontSize',20);
xlim([0, 4])
ylim([0,15])
set(gca, 'XTick',1:3);
set(gca,'XTickLabel',{'Last Loc','Current Loc','Other Locs'});
xlabel('Reward location categories');
ylabel('Percentage of time spent at each location (%)');
title('Pre-running Probe')


% Plot the time spent at (1) last reward location; (2)
% current reward location; (3) Other locations
ind = find(~isnan(sum(time_rewardloc_3,2)));
time_rewardloc_3 = time_rewardloc_3(ind,:);
time_rewardloc_3_mean = mean(time_rewardloc_3);
time_rewardloc_3_sem = std(time_rewardloc_3)./sqrt(size(time_rewardloc_3,1));

ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
bar_color = [100,180,220]./255;
hold on
bar(1:3,time_rewardloc_3_mean,'FaceColor',bar_color,'EdgeColor', bar_color);
h11=errorbar(1:3,time_rewardloc_3_mean,time_rewardloc_3_sem,'color',[0 0 0]./255,'LineStyle','none','LineWidth',2);
errorbar_tick(h11, 0);
hold off
set(gca, 'FontSize',20);
xlim([0,4])
ylim([0,5])
set(gca, 'XTick',1:3);
set(gca, 'YTick',0:5);
set(gca,'XTickLabel',{'Last Loc','Current Loc','Other Locs'});
xlabel('Reward location categories');
ylabel('Time spent at each location (%)');
title('Pre-running Probe')

%% data in Test session and Post-test session
% calculate the stops count for correct and error trials in Test and Post-test trials

%clear

cd(pwdir)

directories_allData_v1

%parentfd = fileparts(mfilename('fullpath'));
%dirfolder = [parentfd,'\GroupData\'];
dirfolder = 'E:\ColginLab\Data Analysis\GroupData';
file_input = 'Data_angle_ontrack.mat';
file_output = 'data_behavior_errorcount.mat';

id_rewardloc = 6:15;  % this is very important! Keep it the same for each rat in the experiment

Loc_rlt_test = [];
Loc_rlt_posttest = [];
Ratid_test = [];
Ratid_posttest = [];
for ns = 1:isession
    path_ns = pathRats{ns};
    cd(path_ns);
    
    % read reward location information
    trackdata_ns = trackdata{ns};
    load(trackdata_ns,'Sign_correct_test','Ind_rewardloc_test',...
        'Sign_correct_posttest','Ind_rewardloc_posttest',...
        'Ang_RewardLoc_ontrack','N_loc','Ind_rewardloc')
    
    % 1. Calculate: counts of correct/error stops in Test trials
    Loc_rlt_test_ns = Ind_rewardloc_test{1}-Ind_rewardloc;
    Loc_rlt_test = [Loc_rlt_test;Loc_rlt_test_ns];
    Ratid_test = [Ratid_test;ones(length(Loc_rlt_test_ns),1)*Ind_Rat(ns)];
    
    % 2. Calculate: counts of correct/error stops in PostTest trials
    Loc_rlt_posttest_ns = Ind_rewardloc_posttest{1}-Ind_rewardloc;
    Loc_rlt_posttest = [Loc_rlt_posttest;Loc_rlt_posttest_ns];
    Ratid_posttest = [Ratid_posttest;ones(length(Loc_rlt_posttest_ns),1)*Ind_Rat(ns)];
end

cd(dirfolder)
save(file_output,'Ratid_test','Ratid_posttest','Loc_rlt_test','Loc_rlt_posttest');

%% Plots in Test session and Post-test session
file_input = 'data_behavior_errorcount.mat';
%%load([dirfolder,file_input])

Plot_rat = [139,150,149,148];
ind = find(ismember(Ratid_test,Plot_rat) > 0);
Loc_rlt_test = Loc_rlt_test(ind,:);
ind = find(ismember(Ratid_posttest,Plot_rat) > 0);
Loc_rlt_posttest = Loc_rlt_posttest(ind,:);

% Test trials
Loc_rlt = Loc_rlt_test;
title0 = 'Test Probe';
% Post-test trials
Loc_rlt = Loc_rlt_posttest;
title0 = 'Post-Test Probe';
% Post and Post-test trials
Loc_rlt = [Loc_rlt_test;Loc_rlt_posttest];
title0 = 'Test and Post-Test Probe';

% plot errors for -2,-1,0,1,2,non-stop
count_correct = [];
for i = -2:1:2
    ind = find(Loc_rlt == i);
    count_correct = [count_correct,length(ind)];
end
ind = find(isnan(Loc_rlt));
count_correct = [count_correct,length(ind)];
sum(count_correct)

ffa = figure('Units','normalized','Position',[0 0 0.6 0.6]);
bar_color = [100,180,220]./255;
hold on
bar(-2:1:2,count_correct(1:5),'FaceColor',bar_color,'EdgeColor', bar_color);
bar(4,count_correct(6),'FaceColor',bar_color,'EdgeColor', bar_color);
hold off
set(gca, 'FontSize',20);
xlim([-3,5])
set(gca, 'XTick',-2:4);
set(gca,'XTickLabel',{'Err-2','Err-1','Correct','Err+1','Err+2','','Nonstop'});
xlabel('Reward location #');
ylabel('Stop counts');
title(title0)

ffa = figure('Units','normalized','Position',[0 0 0.6 0.6]);
bar_color = [100,180,220]./255;
hold on
bar(-2:1:2,100*count_correct(1:5)./sum(count_correct),'FaceColor',bar_color,'EdgeColor', bar_color);
bar(4,100*count_correct(6)./sum(count_correct),'FaceColor',bar_color,'EdgeColor', bar_color);
hold off
set(gca, 'FontSize',20);
xlim([-3,5])
set(gca, 'XTick',-2:4);
set(gca,'XTickLabel',{'Err-2','Err-1','Correct','Err+1','Err+2','','Nonstop'});
xlabel('Reward location #');
ylabel('Percentage of stop counts (%)');
title(title0)

% plot each error locations
min_loc = min(Loc_rlt);
max_loc = max(Loc_rlt);

count_correct = [];
for i = min_loc:1:max_loc
    ind = find(Loc_rlt == i);
    count_correct = [count_correct,length(ind)];
end
ind = find(isnan(Loc_rlt));
count_correct = [count_correct,length(ind)];
sum(count_correct)

ffa = figure('Units','normalized','Position',[0 0 0.6 0.6]);
bar_color = [100,180,220]./255;
hold on
bar(min_loc:max_loc,count_correct(1:7),'FaceColor',bar_color,'EdgeColor', bar_color);
bar(5,count_correct(8),'FaceColor',bar_color,'EdgeColor', bar_color);
hold off
set(gca, 'FontSize',20);
xlim([-4,7])
set(gca, 'XTick',-3:5);
set(gca,'XTickLabel',{'Err-3','Err-2','Err-1','Correct','Err+1','Err+2','Err+3','','Nonstop'});
xlabel('Reward location #');
ylabel('Stop counts');
title(title0)

ffa = figure('Units','normalized','Position',[0 0 0.6 0.6]);
bar_color = [100,180,220]./255;
hold on
bar(min_loc:max_loc,100*count_correct(1:7)./sum(count_correct),'FaceColor',bar_color,'EdgeColor', bar_color);
bar(5,100*count_correct(8)./sum(count_correct),'FaceColor',bar_color,'EdgeColor', bar_color);
hold off
set(gca, 'FontSize',20);
xlim([-4,7])
set(gca, 'XTick',-3:5);
set(gca,'XTickLabel',{'Err-3','Err-2','Err-1','Correct','Err+1','Err+2','Err+3','','Nonstop'});
xlabel('Reward location #');
ylabel('Percentage of stop counts (%)');
title(title0)

FiguresDir = 'E:\ColginLab\Figures\Figure1';

saveas(ffa,[FiguresDir, '\Behavior_barGraph'],'epsc')
saveas(ffa,[FiguresDir, '\Behavior_barGraph'],'fig')
close all

cd(pwdir)

end