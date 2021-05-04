%% Group summed Pxn in the stop zone, check if it can predict correct vs. error trials
% modified on 6/17/2020
%% 
%disp(strcat('ntrial = ',num2str(nl)));
ang_reward = Ang_RewardLoc_ontrack(Ind_rewardloc);
ang_vel_limit = 5/(Diam_inner/2);

% find each reward location zone
loc_bins = scores_sample{nl,5};
df_reward = mean(diff(Ang_RewardLoc_ontrack));

% Method1: get reward location borders [ half zones around the reward location]
% ang_reward_zone_border = mean([Ang_RewardLoc_ontrack(1:end-1),Ang_RewardLoc_ontrack(2:end)],2);
% ang_reward_zone_border = [Ang_RewardLoc_ontrack(1)-df_reward;ang_reward_zone_border];
% ang_reward_zone_border = [ang_reward_zone_border;Ang_RewardLoc_ontrack(end)+df_reward];
% for i = 1:length(ang_reward_zone_border)
%     [~,ang_reward_zone_border_bin(i,1)] = min(abs(ang_reward_zone_border(i)-loc_bins));
% end
% ang_reward_zone = [ang_reward_zone_border(1:end-1),ang_reward_zone_border(2:end)];
% ang_reward_zone_bin = [ang_reward_zone_border_bin(1:end-1),ang_reward_zone_border_bin(2:end)];

% Method2: get reward location borders [one zone before the reward location]
ang_reward_zone_border = [Ang_RewardLoc_ontrack(1)-df_reward;Ang_RewardLoc_ontrack];
for i = 1:length(ang_reward_zone_border)
    [~,ang_reward_zone_border_bin(i,1)] = min(abs(ang_reward_zone_border(i)-loc_bins));
end
ang_reward_zone = [ang_reward_zone_border(1:end-1),ang_reward_zone_border(2:end)];
ang_reward_zone_bin = [ang_reward_zone_border_bin(1:end-1),ang_reward_zone_border_bin(2:end)];

%% Sample Lap data

stop_loc_id = Ind_rewardloc_sample{1}(nl);
if ~isnan(stop_loc_id)
    stop_minus1_loc = Ang_RewardLoc_ontrack(stop_loc_id-1);
else
    stop_minus1_loc = Ang_RewardLoc_ontrack(Ind_rewardloc-1);
    stop_loc_id = Ind_rewardloc;
end

vel_sample = scores_sample{nl,8}(1,:);
angvel_sample = scores_sample{nl,8}(4,:);

% find actual locations bins
bin_actual_sample = scores_sample{nl,4};
loc_actual_sample_binned = loc_bins(bin_actual_sample);

% in each time bin, calculate averaged Pxn in the reward zone
% (exclude the +1 bins around the rat's current bin)
t_start_sample = scores_sample{nl,6};
Pxn_sample = scores_sample{nl,3};
Pxn_sample_stop_all = nan(1,size(Pxn_sample,2));
ind1 = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2));
ind2 = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4));
ind = unique([ind1,ind2]);
reward_zone_range = ang_reward_zone_bin(stop_loc_id,:);
Pxn_sample_stop_all(1,ind) = max(Pxn_sample(reward_zone_range(1):reward_zone_range(2),ind));

%% group sample lap data

% Before getting stop-1 location
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2) & ...
    loc_actual_sample_binned' < stop_minus1_loc);
Pxn_sample_stop = nan(1,size(Pxn_sample,2));
angvel_sample_stop = nan(1,size(Pxn_sample,2));
Pxn_sample_stop(1,ind) = Pxn_sample_stop_all(1,ind);
angvel_sample_stop(1,ind) = angvel_sample(ind);
nwin_sample_stop = length(ind);

% Before getting stop-1 location
% + running speed > 5cm/s
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2) & ...
    loc_actual_sample_binned' < stop_minus1_loc &...
    angvel_sample >= ang_vel_limit);
Pxn_sample_stop_vel = nan(1,size(Pxn_sample,2));
angvel_sample_stop_vel = nan(1,size(Pxn_sample,2));
Pxn_sample_stop_vel(1,ind) = Pxn_sample_stop_all(1,ind);
angvel_sample_stop_vel(1,ind) = angvel_sample(ind);
nwin_sample_stop_vel = length(ind);

% Before getting stop-1 location
% + running speed > 5cm/s
% + approaching reward location (5 locations away from reward location):
if ~isnan(Sign_correct_sample{1}(nl))
    ang_reward_approach_sample = Ang_RewardLoc_ontrack(max(Ind_rewardloc_sample{1}(nl)-ind_approach,1));
else
    ang_reward_approach_sample = Ang_RewardLoc_ontrack(max(Ind_rewardloc-ind_approach,1));
end
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2) & ...
    loc_actual_sample_binned' < stop_minus1_loc &...
    angvel_sample >= ang_vel_limit &...
    loc_actual_sample_binned' >= ang_reward_approach_sample);
ind0 = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2) & ...
    loc_actual_sample_binned' <= ang_reward_approach_sample);  % find the last point around -5 location
if ~isempty(ind0)
    ind0 = ind0(end);
    ind = ind(ind>=ind0);
end
Pxn_sample_stop_approach_vel = nan(1,size(Pxn_sample,2));
angvel_sample_stop_approach_vel = nan(1,size(Pxn_sample,2));
Pxn_sample_stop_approach_vel(1,ind) = Pxn_sample_stop_all(1,ind);
angvel_sample_stop_approach_vel(1,ind) = angvel_sample(ind);
nwin_sample_stop_approach_vel = length(ind);

%% Test Lap data

stop_loc_id = Ind_rewardloc_test{1}(nl);
if ~isnan(stop_loc_id)
    stop_minus1_loc = Ang_RewardLoc_ontrack(stop_loc_id-1);
else
    stop_minus1_loc = Ang_RewardLoc_ontrack(Ind_rewardloc-1);
    stop_loc_id = Ind_rewardloc;
end

vel_test = scores_test{nl,8}(1,:);
angvel_test = scores_test{nl,8}(4,:);

% find actual locations bins
bin_actual_test = scores_test{nl,4};
loc_actual_test_binned = loc_bins(bin_actual_test);

% in each time bin, calculate averaged Pxn in the reward zone
% (exclude the +1 bins around the rat's current bin)
t_start_test = scores_test{nl,6};
Pxn_test = scores_test{nl,3};
Pxn_test_stop_all = nan(1,size(Pxn_test,2));
ind1 = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2));
ind2 = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4));
ind = unique([ind1,ind2]);
reward_zone_range = ang_reward_zone_bin(stop_loc_id,:);
Pxn_test_stop_all(1,ind) = max(Pxn_test(reward_zone_range(1):reward_zone_range(2),ind));
%% group test lap data

% Before getting stop-1 location
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2) & ...
    loc_actual_test_binned' < stop_minus1_loc);
Pxn_test_stop = nan(1,size(Pxn_test,2));
angvel_test_stop = nan(1,size(Pxn_test,2));
Pxn_test_stop(1,ind) = Pxn_test_stop_all(1,ind);
angvel_test_stop(1,ind) = angvel_test(ind);
nwin_test_stop = length(ind);

% Before getting stop-1 location
% + running speed > 5cm/s
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2) & ...
    loc_actual_test_binned' < stop_minus1_loc &...
    angvel_test >= ang_vel_limit);
Pxn_test_stop_vel = nan(1,size(Pxn_test,2));
angvel_test_stop_vel = nan(1,size(Pxn_test,2));
Pxn_test_stop_vel(1,ind) = Pxn_test_stop_all(1,ind);
angvel_test_stop_vel(1,ind) = angvel_test(ind);
nwin_test_stop_vel = length(ind);

% Before getting stop-1 location
% + running speed > 5cm/s
% + approaching reward location (5 locations away from reward location):
if ~isnan(Sign_correct_test{1}(nl))
    ang_reward_approach_test = Ang_RewardLoc_ontrack(max(Ind_rewardloc_test{1}(nl)-ind_approach,1));
else
    ang_reward_approach_test = Ang_RewardLoc_ontrack(max(Ind_rewardloc-ind_approach,1));
end
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2) & ...
    loc_actual_test_binned' < stop_minus1_loc &...
    angvel_test >= ang_vel_limit &...
    loc_actual_test_binned' >= ang_reward_approach_test);
ind0 = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2) & ...
    loc_actual_test_binned' <= ang_reward_approach_test);  % find the last point around -5 location
if ~isempty(ind0)
    ind0 = ind0(end);
    ind = ind(ind>=ind0);
end
Pxn_test_stop_approach_vel = nan(1,size(Pxn_test,2));
angvel_test_stop_approach_vel = nan(1,size(Pxn_test,2));
Pxn_test_stop_approach_vel(1,ind) = Pxn_test_stop_all(1,ind);
angvel_test_stop_approach_vel(1,ind) = angvel_test(ind);
nwin_test_stop_approach_vel = length(ind);

%% Save data
sign_sample = 1;
sign_test = 2;


% add 2 cols at the end of these matrics
data0_prereward = [ns,nl,Ind_rewardloc,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    nanmean(angvel_sample_stop),...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    nanmean(angvel_test_stop),...
    nwin_sample_stop,nwin_test_stop];
data0_vel_prereward = [ns,nl,Ind_rewardloc,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    nanmean(angvel_sample_stop_vel),...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    nanmean(angvel_test_stop_vel),...
    nwin_sample_stop_vel,nwin_test_stop_vel];
data0_approach_vel_prereward = [ns,nl,Ind_rewardloc,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    nanmean(angvel_sample_stop_approach_vel),...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    nanmean(angvel_test_stop_approach_vel),...
    nwin_sample_stop_approach_vel,nwin_test_stop_approach_vel];
