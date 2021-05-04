% edited on 12/3/2018
% Behavioral data
% running speed across positions in sample and test trials

%% data in sample-test session
% check the running speed in samples and tests

clear

directories_allData_v1
Parentfd = fileparts(mfilename('fullpath'));
dirfolder = [Parentfd,'\GroupData\'];
file_input = 'Data_angle_ontrack.mat';
file_output = 'data_behavior_speed.mat';
file_output_Stats = [Parentfd,'\GroupData Figures\Stats_Behavior_speed.mat'];

id_rewardloc = 6:15;  % this is very important! Keep it the same for each rat in the experiment

vel_all = []; % running speed
vel_ang_all = []; % angular velocity
Ratid = [];
n_seg = 10;  % use script_meanspeed_v1
% n_seg = 4;  % use script_meanspeed_v2

speed_samp = nan(isession,n_seg);
speed_test = nan(isession,n_seg);
speed_samp_corr = nan(isession,n_seg);
speed_test_corr = nan(isession,n_seg);
speed_samp_err = nan(isession,n_seg);
speed_test_err = nan(isession,n_seg);
for ns = 1:isession
    path_ns = pathRats{ns};
    cd(path_ns);
    
    % load angle data in sample trials
    load(file_input)
    trackdata_ns = trackdata{ns};
    load(trackdata_ns,'Ang_RewardLoc_ontrack','N_loc','Ind_rewardloc')
    load(trackdata_ns,'Sign_correct_test',...
        'Ang_StartZone_depart_ontrack',...
        'Ang_StartZone_arrive_ontrack',...
        'Ang_RewardLoc_ontrack','Diam_inner',...
        'Ts_sample','Ts_test','Ind_rewardloc','Ind_rewardloc_sample','Ind_rewardloc_test');
    
    nlap = size(data_angle{1,2},1);
    speed_seg_samp = nan(nlap,n_seg);
    speed_seg_test = nan(nlap,n_seg);
    
    for nl = 1:nlap
        if isnan(Ind_rewardloc_sample{1}(nl)) | isnan(Ind_rewardloc_test{1}(nl))
            speed_seg_samp(nl,:) = nan(1,n_seg);
            speed_seg_test(nl,:) = nan(1,n_seg);
        else
            % sample
            Ind_rewardloc0 = Ind_rewardloc_sample;
            Ts0 = Ts_sample;
            data_angle0 = data_angle{1,2};
            script_meanspeed_v1;  % v1: make segments evenly from Loc1 to stop location
%             script_meanspeed_v2;  % v2: from stop location -3 to stop location
            speed_seg_samp(nl,:) = speed_seg;
            
            % test
            Ind_rewardloc0 = Ind_rewardloc_test;
            Ts0 = Ts_test;
            data_angle0 = data_angle{1,3};
            script_meanspeed_v1;  % v1: make segments evenly from Loc1 to stop location
%             script_meanspeed_v2;  % v2: from stop location -3 to stop location
            speed_seg_test(nl,:) = speed_seg;
        end
    end
    speed_samp(ns,:) = nanmean(speed_seg_samp);
    speed_test(ns,:) = nanmean(speed_seg_test);
    
    ind_corr = find(Sign_correct_test{1} == 1);
    if ~isempty(ind_corr)
        speed_samp_corr(ns,:) = nanmean(speed_seg_samp(ind_corr,:));
        speed_test_corr(ns,:) = nanmean(speed_seg_test(ind_corr,:));
    end
    ind_err = find(Sign_correct_test{1} == 0);
    if ~isempty(ind_err)
        speed_samp_err(ns,:) = nanmean(speed_seg_samp(ind_err,:));
        speed_test_err(ns,:) = nanmean(speed_seg_test(ind_err,:));
    end
    
    Ratid(ns,1) = Ind_Rat(ns);
end

cd(dirfolder)

% save(file_output,'Ratid','speed_samp','speed_test','speed_samp_corr','speed_test_corr',...
%     'speed_samp_err','speed_test_err');


%% PLOTS
data1 = speed_samp; data2 = speed_test;
data1 = speed_samp_corr; data2 = speed_test_corr;
data1 = speed_samp_err; data2 = speed_test_err;

ind = find(~isnan(mean(data1-data2,2)));
data1 = data1(ind,:); data2 = data2(ind,:);

sample_color = [255,165,0]./255;
test_color = [111,57,214]./255;
ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
data1_mean = mean(data1);
data1_sem = std(data1)./sqrt(size(data1,1));
data2_mean = mean(data2);
data2_sem = std(data2)./sqrt(size(data2,1));
hold on
h1 = errorbar(1:length(data1_mean),data1_mean,data1_sem,'s-','color',sample_color,'LineWidth',2);
h2 = errorbar(1:length(data2_mean),data2_mean,data2_sem,'s-','color',test_color,'LineWidth',2);
%errorbar_tick(h1, 0);
%errorbar_tick(h2, 0);
legend('Samples','Tests','Location','SouthEast')
hold off
set(gca,'XTick',[1:length(data1_mean)]);
xlim([0,length(data1_mean)+1])
ylim([15,50])
set(gca,'fontsize',20);
xlabel('Position segments');
ylabel('Running speed (cm/s)');

set(gca,'XTickLabel',{'Loc1','','','','','','','','','Stop Loc'});
% set(gca,'XTickLabel',{'-3','-2','-1','Stop Loc'});

title('Correct trials')
title('Error trials')
close all
%% Matlab Stats
StatsData_LineGraph = array2table([speed_test_corr, speed_test_err],'VariableNames',{'Loc1_crt','Loc2_crt','Loc3_crt','Loc4_crt','Loc5_crt',...
                                                                                     'Loc6_crt','Loc7_crt','Loc8_crt','Loc9_crt','Loc10_crt',...
                                                                                     'Loc1_err','Loc2_err','Loc3_err','Loc4_err','Loc5_err',...
                                                                                     'Loc6_err','Loc7_err','Loc8_err','Loc9_err','Loc10_err'});
within_fact = fullfact([10,2]);
within_tbl2 = array2table(within_fact,'VariableNames',{'Location','CrtErr'});
within_tbl2.Location = categorical(within_tbl2.Location);
within_tbl2.CrtErr = categorical(within_tbl2.CrtErr);
rm1 = fitrm(StatsData_LineGraph,'Loc1_crt-Loc10_err ~ 1','WithinDesign',within_tbl2);
StatsOut = ranova(rm1,'WithinModel','CrtErr');

save(file_output_Stats,'StatsOut','StatsData_LineGraph')
%% SPSS STATS
% 2 way ANOVA, with locations being repeated factor
data1 = data1(:,1:10); data2 = data2(:,1:10);
N = size(data1,1);
data0 = [ [ones(N,1)*1;ones(N,1)*2],[data1;data2] ]; % locations

% 2 way ANOVA, with sample and test being repeated factor
data1 = data1(:,1:10); data2 = data2(:,1:10);
data1 = reshape(data1',size(data1,1)*size(data1,2),1); % running speed
data2 = reshape(data2',size(data2,1)*size(data2,2),1);
data0 = [repmat([1:n_seg]',N,1),data1,data2]; % locations

% mixed model in SPSS
data1 = reshape(speed_samp',size(speed_samp,1)*size(speed_samp,2),1); % running speed
data2 = reshape(speed_test',size(speed_test,1)*size(speed_test,2),1);

data1 = [repmat([1:size(speed_samp,2)]',size(speed_samp,1),1),data1]; % locations
data2 = [repmat([1:size(speed_test,2)]',size(speed_test,1),1),data2];

data00 = repmat([1:size(speed_samp,1)],size(speed_samp,2),1); % sessions
data00 = reshape(data00,size(data00,1)*size(data00,2),1);
data1 = [data00, data1];
data2 = [data00, data2];

data1 = [ones(size(data1,1),1), data1];% sample/test
data2 = [ones(size(data2,1),1), data2];

data0 = [data1;data2];