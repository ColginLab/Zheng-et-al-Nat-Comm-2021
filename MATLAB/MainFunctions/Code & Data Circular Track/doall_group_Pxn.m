% Edited on 07/28/2017
% Modified from 'doall_group_decodingErr.m'
% Group the P(x|n) regardless of theta cycles
% 
% Since we need continuous decoding for a running lap,
% the slopes of theta sequences can not be precisely measured
% due to missing bins and incorrect theta cycle cutoff.
%
% Modified on 05/07/2019
% add measurements group_Pxn_approach_vel_prereward_2bin and group_Pxn_approach_vel_postreward_2bin

%% settings
function doall_group_Pxn(DataDir, AnalysisDir)

%fig_dir = 'F:\MATLAB\MainFunctions\Code & Data Circular Track\';
ind_approach = 5;  % approaching reward location (3 locations away from reward location):
nbin_sum = 5;  % CHANGE: the number to sum up probabilities

% (1-2) use the overall place field excluding pre-ruuning trials as a decoder
% downsample place cells so that the field distribution equally on the track
file_input = 'BayesData_CircMap_v1_dt40ms_10ms_1cells_1spk_v2_ds.mat';
TTList0 = 'TTList_dCA1_pyr_ds.txt';
file_input_speed = 'Data_angle_ontrack.mat';  % input the speed information
fig_folder ='GroupData\';

% file_output1 = 'group_PxnSum_5LocAway_25cells_20200617_5bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
file_output1 = 'group_PxnSum_5LocAway_25cells_20201217_5bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
file_output2 = 'group_Pxn_forErr_5LocAway_25cells_20200617_5bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
file_output3 = 'group_Pxn_stopzone_5LocAway_25cells_20200623_5bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells

code_version = 'v3';

file_output1 = strcat(AnalysisDir,file_output1);
file_output2 = strcat(AnalysisDir,file_output2);
file_output3 = strcat(AnalysisDir,file_output3);

ncell_threshold = 25; % Only use data with more than 30 cells

dt = 0.04; % 40 ms window

%% Group data
directories_allData_v1

Nlaps = 0;
Ns = 0;

group_Pxn_prereward = [];  % Pxn using all location bins before getting reward
group_Pxn_vel_prereward = [];  % Pxn using all location bins with vel>5cm/s
group_Pxn_approach_vel_prereward = [];%  Pxn using all location bins with vel>5cm/s, and the rat is approaching the reward location
group_Pxn_approach_vel_prereward_2bin = [];%  Pxn using all location bins with vel>5cm/s, and the rat is approaching the reward location, from current location+2bins
group_Pxn_postreward = [];  % Pxn using all location bins after getting reward
group_Pxn_vel_postreward = [];  % Pxn using all location bins with vel>5cm/s
group_Pxn_approach_vel_postreward = [];%  Pxn using all location bins with vel>5cm/s, and the rat is approaching the reward location
group_Pxn_approach_vel_postreward_2bin = [];%  Pxn using all location bins with vel>5cm/s, and the rat is approaching the reward location, from current location+2bins

group_Pxn_prestop = [];
group_Pxn_vel_prestop = [];
group_Pxn_approach_vel_prestop = [];
group_Pxn_approach_vel_prestop_2bin = [];
group_Pxn_poststop = [];
group_Pxn_vel_poststop = [];
group_Pxn_approach_vel_poststop = [];
group_Pxn_approach_vel_poststop_2bin = [];

group_Pxn_stopzone_info = [];
group_Pxn_stopzone_vel_info = [];
group_Pxn_stopzone_approach_vel_info = [];
group_Pxn_stopzone_samp = cell(1,8);
group_Pxn_stopzone_samp_vel = cell(1,8);
group_Pxn_stopzone_samp_approach_vel = cell(1,8);
group_Pxn_stopzone_test = cell(1,8);
group_Pxn_stopzone_test_vel = cell(1,8);
group_Pxn_stopzone_test_approach_vel = cell(1,8);

path_group = cell(1,1);
for ns = 1:isession
    path_ns = pathRats{ns};
    cd(path_ns);
    
    Ncell = getnumberofcells_cz_v1(TTList0);
    
    if exist(file_input,'file')>0 && Ncell >= ncell_threshold
        disp(['Currently in ' pwd])
        Ns = Ns+1;
        path_group{Ns} = pathRats{ns};
        
        load (file_input);
        load(file_input_speed,'data_angle');
        
        trackdata_ns = trackdata{ns};
        load(trackdata_ns);
        
        % Only focus on sample-test laps for now
        scores_sample = scores{1,2};
        scores_test = scores{1,3};
        nlaps = size(scores_sample,1);
        
        % Limit the time period between RewardLoc(1) and stop location
        ts_start_stop_sample = Ts_sample{1}./ 1000000;
        ts_start_stop_test = Ts_test{1}./ 1000000;
        
        ang_vel_limit = 5/(Diam_inner/2); % 5cm/s
        
        for nl = 1:nlaps
            Nlaps = Nlaps+1;
            ts_start_stop_sample_nl_0 = ts_start_stop_sample(nl,:);
            ts_start_stop_test_nl_0 = ts_start_stop_test(nl,:);
            ts_start_stop_sample_nl = nan(1,4);
            ts_start_stop_test_nl = nan(1,4);
            
            %% start time: when the rat passed reward location #1
            ang_reward1 = Ang_RewardLoc_ontrack(1);
            
            % Sample trials
            ind = find(data_angle{1,2}{nl}(1:end-1,2) <= ang_reward1 &...
                data_angle{1,2}{nl}(2:end,2) >= ang_reward1 & diff(data_angle{1,2}{nl}(:,2))<pi);
            ind = ind(end)+1;
            ts_start_stop_sample_nl(1) = data_angle{1,2}{nl}(ind,1);
            
            % Test trials
            ind = find(data_angle{1,3}{nl}(1:end-1,2) <= ang_reward1 &...
                data_angle{1,3}{nl}(2:end,2) >= ang_reward1 & diff(data_angle{1,3}{nl}(:,2))<pi);
            ind = ind(end)+1;
            ts_start_stop_test_nl(1) = data_angle{1,3}{nl}(ind,1);
            
            %% reward time: when the rat passed correct/incorrect reward location
            
            % Sample trials
            if ~isnan(Sign_correct_sample{1}(nl))
                ind = find(data_angle{1,2}{nl}(:,1) <= ts_start_stop_sample_nl_0(2) &...
                    data_angle{1,2}{nl}(:,4) >= ang_vel_limit);
                ind = ind(end);
                ts_start_stop_sample_nl(2) = data_angle{1,2}{nl}(ind,1);
            else
                ang_reward_sample_border = Ang_RewardLoc_ontrack(Ind_rewardloc,1);
                ind = find(data_angle{1,2}{nl}(1:end-1,2) <= ang_reward_sample_border &...
                    data_angle{1,2}{nl}(2:end,2) >= ang_reward_sample_border & diff(data_angle{1,2}{nl}(:,2))<pi);
                ind = ind(1);
                ts_start_stop_sample_nl(2) = min(data_angle{1,2}{nl}(ind,1),ts_start_stop_sample_nl_0(2));
            end
            
            % Test trials
            if ~isnan(Sign_correct_test{1}(nl))
                ind = find(data_angle{1,3}{nl}(:,1) <= ts_start_stop_test_nl_0(2) &...
                    data_angle{1,3}{nl}(:,4) >= ang_vel_limit);
                ind = ind(end);
                ts_start_stop_test_nl(2) = data_angle{1,3}{nl}(ind,1);
            else
                ang_reward_test_border = Ang_RewardLoc_ontrack(Ind_rewardloc,1);
                ind = find(data_angle{1,3}{nl}(1:end-1,2) <= ang_reward_test_border &...
                    data_angle{1,3}{nl}(2:end,2) >= ang_reward_test_border & diff(data_angle{1,3}{nl}(:,2))<pi);
                ind = ind(1);
                ts_start_stop_test_nl(2) = min(data_angle{1,3}{nl}(ind,1),ts_start_stop_test_nl_0(2));
            end
            
            %% Re-start time: when the rat restarts after getting rewards
            % Starting from reward location +1
            
            % Sample trials
            ang_reward_sample_border = Ang_RewardLoc_ontrack(Ind_rewardloc+1,1);
            ind = find(data_angle{1,2}{nl}(:,1) > ts_start_stop_sample_nl(2) &...
                data_angle{1,2}{nl}(:,2) >= ang_reward_sample_border &...
                data_angle{1,2}{nl}(:,4) >= ang_vel_limit);
            ind = ind(1);
            ts_start_stop_sample_nl(3) = data_angle{1,2}{nl}(ind,1);
            
            % Test trials
            if ~isnan(Sign_correct_test{1}(nl))
                ind_rewardloc = Ind_rewardloc_test{1}(nl);
            else
                ind_rewardloc = Ind_rewardloc;
            end
            ang_reward_test_border = Ang_RewardLoc_ontrack(ind_rewardloc+1,1);
            ind = find(data_angle{1,3}{nl}(:,1) > ts_start_stop_test_nl(2) &...
                data_angle{1,3}{nl}(:,2) >= ang_reward_test_border &...
                data_angle{1,3}{nl}(:,4) >= ang_vel_limit);
            ind = ind(1);
            ts_start_stop_test_nl(3) = data_angle{1,3}{nl}(ind,1);
            
            
            %% end time: when the rat reaches the last location
            ang_rewardlast = Ang_RewardLoc_ontrack(end);
            
            % Sample trials
            ind = find(data_angle{1,2}{nl}(1:end-1,2) <= ang_rewardlast &...
                data_angle{1,2}{nl}(2:end,2) >= ang_rewardlast &...
                data_angle{1,2}{nl}(1:end-1,1) > ts_start_stop_sample_nl(3) &...
                data_angle{1,2}{nl}(1:end-1,1) <= ts_start_stop_sample_nl_0(3));
            ind = ind(1);
            ts_start_stop_sample_nl(4) = data_angle{1,2}{nl}(ind,1);
            
            % Test trials
            ind = find(data_angle{1,3}{nl}(1:end-1,2) <= ang_rewardlast &...
                data_angle{1,3}{nl}(2:end,2) >= ang_rewardlast &...
                data_angle{1,3}{nl}(1:end-1,1) > ts_start_stop_test_nl(3) &...
                data_angle{1,3}{nl}(1:end-1,1) <= ts_start_stop_test_nl_0(3));
            ind = ind(1);
            ts_start_stop_test_nl(4) = data_angle{1,3}{nl}(ind,1);
            
          %% Group the sum Pxn before and after getting reward
            group_PxnSum;
            group_Pxn_prereward(Nlaps,:) = data0_prereward;
            group_Pxn_vel_prereward(Nlaps,:) = data0_vel_prereward;
            group_Pxn_approach_vel_prereward(Nlaps,:) = data0_approach_vel_prereward;
            group_Pxn_approach_vel_prereward_2bin(Nlaps,:) = data0_approach_vel_prereward_2bin;
            group_Pxn_postreward(Nlaps,:) = data0_postreward;
            group_Pxn_vel_postreward(Nlaps,:) = data0_vel_postreward;
            group_Pxn_approach_vel_postreward(Nlaps,:) = data0_approach_vel_postreward;
            group_Pxn_approach_vel_postreward_2bin(Nlaps,:) = data0_approach_vel_postreward_2bin;
            
          %% Group the sum Pxn before and after stop
            group_PxnSum_v2;
            group_Pxn_prestop(Nlaps,:) = data0_prestop;
            group_Pxn_vel_prestop(Nlaps,:) = data0_vel_prestop;
            group_Pxn_approach_vel_prestop(Nlaps,:) = data0_approach_vel_prestop;
            group_Pxn_approach_vel_prestop_2bin(Nlaps,:) = data0_approach_vel_prestop_2bin;
            group_Pxn_poststop(Nlaps,:) = data0_poststop;
            group_Pxn_vel_poststop(Nlaps,:) = data0_vel_poststop;
            group_Pxn_approach_vel_poststop(Nlaps,:) = data0_approach_vel_poststop;
            group_Pxn_approach_vel_poststop_2bin(Nlaps,:) = data0_approach_vel_poststop_2bin;
            
          %% Group the averaged Pxn, check if it can predict correct vs. error trials
            group_Pxn_forErr;
            group_Pxn_forErr_prereward_info(Nlaps,:) = data0_prereward;
            group_Pxn_forErr_vel_prereward_info(Nlaps,:) = data0_vel_prereward;
            group_Pxn_forErr_approach_vel_prereward_info(Nlaps,:) = data0_approach_vel_prereward;
            group_Pxn_forErr_prereward_mean{ns,1}(:,nl) = Pxn_sample_prereward_mean;
            group_Pxn_forErr_prereward_mean{ns,2}(:,nl) = Pxn_test_prereward_mean;
            group_Pxn_forErr_vel_prereward_mean{ns,1}(:,nl) = Pxn_sample_vel_prereward_mean;
            group_Pxn_forErr_vel_prereward_mean{ns,2}(:,nl) = Pxn_test_vel_prereward_mean;
            group_Pxn_forErr_approach_vel_prereward_mean{ns,1}(:,nl) = Pxn_sample_approach_vel_prereward_mean;
            group_Pxn_forErr_approach_vel_prereward_mean{ns,2}(:,nl) = Pxn_test_approach_vel_prereward_mean;
            group_Pxn_forErr_prereward_binsum{ns,1}(:,nl) = Pxn_sample_prereward_binsum;
            group_Pxn_forErr_prereward_binsum{ns,2}(:,nl) = Pxn_test_prereward_binsum;
            group_Pxn_forErr_vel_prereward_binsum{ns,1}(:,nl) = Pxn_sample_vel_prereward_binsum;
            group_Pxn_forErr_vel_prereward_binsum{ns,2}(:,nl) = Pxn_test_vel_prereward_binsum;
            group_Pxn_forErr_approach_vel_prereward_binsum{ns,1}(:,nl) = Pxn_sample_approach_vel_prereward_binsum;
            group_Pxn_forErr_approach_vel_prereward_binsum{ns,2}(:,nl) = Pxn_test_approach_vel_prereward_binsum;
            
            %% Group the summed Pxn in the stop zone, check if it can predict correct vs. error trials
            group_Pxn_stopzone;
            group_Pxn_stopzone_info(Nlaps,:) = data0_prereward;
            group_Pxn_stopzone_vel_info(Nlaps,:) = data0_vel_prereward;
            group_Pxn_stopzone_approach_vel_info(Nlaps,:) = data0_approach_vel_prereward;
            group_Pxn_stopzone_samp{ns,nl} = [Pxn_sample_stop;angvel_sample_stop];
            group_Pxn_stopzone_samp_vel{ns,nl} = [Pxn_sample_stop_vel;angvel_sample_stop_vel];
            group_Pxn_stopzone_samp_approach_vel{ns,nl} = [Pxn_sample_stop_approach_vel;angvel_sample_stop_approach_vel];
            group_Pxn_stopzone_test{ns,nl} = [Pxn_test_stop;angvel_test_stop];
            group_Pxn_stopzone_test_vel{ns,nl} = [Pxn_test_stop_vel;angvel_test_stop_vel];
            group_Pxn_stopzone_test_approach_vel{ns,nl} = [Pxn_test_stop_approach_vel;angvel_test_stop_approach_vel];
            
        end
    end
    
    cd ../
end

close all

ns_RatID = nan(size(group_Pxn_prereward,1),2);
ns_id = unique(group_Pxn_prereward(:,1));
for i = ns_id'
    ind = find(group_Pxn_prereward(:,1)==i);
    ns_RatID(ind,1) = group_Pxn_prereward(ind,1);
    ns_RatID(ind,2) = Ind_Rat(i)*ones(length(ind),1);
end

save(file_output1,...
    'group_Pxn_prereward','group_Pxn_vel_prereward','group_Pxn_approach_vel_prereward',...
    'group_Pxn_postreward','group_Pxn_vel_postreward','group_Pxn_approach_vel_postreward',...
    'group_Pxn_approach_vel_prereward_2bin','group_Pxn_approach_vel_postreward_2bin',...
    'path_group','ns_RatID',...
    'group_Pxn_prestop','group_Pxn_vel_prestop','group_Pxn_approach_vel_prestop',...
    'group_Pxn_poststop','group_Pxn_vel_poststop','group_Pxn_approach_vel_poststop',...
    'group_Pxn_approach_vel_prestop_2bin','group_Pxn_approach_vel_poststop_2bin',...
    'path_group','ns_RatID');

% Data structure:
% col1:ns; col2:nl; 
% col3:1(sample); col4:sample correct sign; col5:reward loc id;
% col6:Pxn; col7:angvel_sample
% col8:2(test); col9:test correct sign; col10:reward loc id;
% col11:Pxn; col12:angvel_test

% save(file_output2,...
%     'group_Pxn_forErr_prereward_info',...
%     'group_Pxn_forErr_vel_prereward_info',...
%     'group_Pxn_forErr_approach_vel_prereward_info',...
%     'group_Pxn_forErr_prereward_mean',...
%     'group_Pxn_forErr_vel_prereward_mean',...
%     'group_Pxn_forErr_approach_vel_prereward_mean',...
%     'group_Pxn_forErr_prereward_binsum',...
%     'group_Pxn_forErr_vel_prereward_binsum',...
%     'group_Pxn_forErr_approach_vel_prereward_binsum',...
%     'path_group','ns_RatID');
% 

save(file_output3,...
    'group_Pxn_stopzone_info',...
    'group_Pxn_stopzone_vel_info',...
    'group_Pxn_stopzone_approach_vel_info',...
    'group_Pxn_stopzone_samp',...
    'group_Pxn_stopzone_samp_vel',...
    'group_Pxn_stopzone_samp_approach_vel',...
    'group_Pxn_stopzone_test',...
    'group_Pxn_stopzone_test_vel',...
    'group_Pxn_stopzone_test_approach_vel',...
    'path_group','ns_RatID');
%% do Stats
% file_output1 = 'group_PxnSum_5LocAway_35cells_20170729.mat';
% file_output1 = 'group_PxnSum_5LocAway_35cells_20171006_v2.mat';
% file_output1 = 'group_PxnSum_5LocAway_35cells_20171008_v3.mat'; %3 rats
% file_output1 = 'group_PxnSum_5LocAway_35cells_20190118_v3.mat'; %4 rats
% file_output1 = 'group_PxnSum_5LocAway_35cells_20190507_5bins_v3.mat'; %4 rats, change nbin_sum = 5 in group_PxnSum.m and group_Pxn_forErr.m
% file_output1 = 'group_PxnSum_5LocAway_40cells_20190118_v3.mat'; %4 rats, >40cells
% file_output1 = 'group_PxnSum_5LocAway_35cells_20190507_v3.mat'; %4 rats, >35cells, method(2)
% file_output1 = 'group_PxnSum_5LocAway_35cells_20190516_5bins_v3.mat'; %4 rats, >35cells, method(2)
% file_output1 = 'group_PxnSum_5LocAway_20cells_20190508_10bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
% file_output1 = 'group_PxnSum_5LocAway_20cells_20190508_5bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
% file_output1 = 'group_PxnSum_5LocAway_30cells_20190518_5bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
% file_output1 = 'group_PxnSum_5LocAway_25cells_20190518_5bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
% file_output1 = 'group_PxnSum_3LocAway_25cells_20190526_5bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
file_output1 = 'group_PxnSum_5LocAway_25cells_20201217_5bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
file_output1 = strcat(AnalysisDir,file_output1);

load(file_output1)

% Stats_PxnSum

Stats_PxnSum_CI
%Stats_PxnSum_CI_Fig2  % directly plot final Fig2
%% do Stats
% file_output2 = 'group_Pxn_forErr_5LocAway_35cells_20170731.mat';
% file_output2 = 'group_Pxn_forErr_5LocAway_35cells_20171006_v2.mat';
% file_output2 = 'group_Pxn_forErr_5LocAway_35cells_20171008_v3.mat'; %3 rats
% file_output2 = 'group_Pxn_forErr_5LocAway_35cells_20190118_v3.mat'; %4 rats
% file_output2 = 'group_Pxn_forErr_5LocAway_35cells_20190507_5bins_v3.mat'; %4 rats, change nbin_sum = 5 in group_PxnSum.m and group_Pxn_forErr.m
% file_output2 = 'group_Pxn_forErr_5LocAway_40cells_20190118_v3.mat'; %4 rats, >40cells
% file_output2 = 'group_Pxn_forErr_5LocAway_35cells_20190507_v3.mat'; %4 rats, >35cells, method(2)
% file_output2 = 'group_Pxn_forErr_5LocAway_35cells_20190516_5bins_v3.mat'; %4 rats, >35cells, method(2)
% file_output2 = 'group_Pxn_forErr_5LocAway_20cells_20190508_10bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
% file_output2 = 'group_Pxn_forErr_5LocAway_20cells_20190508_5bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
% file_output2 = 'group_Pxn_forErr_5LocAway_30cells_20190518_10bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
% file_output2 = 'group_Pxn_forErr_5LocAway_25cells_20190518_10bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells

% load(file_output2)
% Stats_Pxn_forErr

%% do Stats
file_output3 = 'group_Pxn_stopzone_5LocAway_25cells_20200623_5bins_v3_ds.mat';
%file_output3 = 'group_Pxn_stopzone_5LocAway_25cells_20200617_5bins_v3_ds.mat';
file_output3 = strcat(AnalysisDir,file_output3);

load(file_output3)
%Stats_Pxn_stopzone

end