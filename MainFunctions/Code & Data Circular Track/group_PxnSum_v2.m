%% Group Sum Pxn of continuous decodings in each sample-test trial
% Modified from group_PxnSum.m on 5/26/2019
% get sum pxn until stop location

% do 8 measurements:
% 1. data0_prereward: 
%    decoding errors before getting reward
% 2. data0_vel_prereward:
%    decoding errors with vel>5cm/s
% 3. data0_approach_vel_prereward:
%    decoding errors with vel>5cm/s, and the rat is approaching the reward location
% 4. data0_approach_vel_prereward_2bin:
%    decoding errors with vel>5cm/s, and the rat is approaching the reward location, + 2bins ahead of the rat
% 5. data0_postreward: 
%    decoding errors after getting reward
% 6. data0_vel_postreward:
%    decoding errors with vel>5cm/s
% 7. data0_approach_vel_postreward:
%    decoding errors with vel>5cm/s, and the rat is approaching the animal box
% 8. data0_approach_vel_postreward_2bin:
%    decoding errors with vel>5cm/s, and the rat is approaching the animal box, + 2bins ahead of the rat

%% 
%disp(strcat('ntrial = ',num2str(nl)))
if ~isnan(Sign_correct_sample{1}(nl))
    ang_stop_samp = Ang_RewardLoc_ontrack(Ind_rewardloc_sample{1}(nl));
else
    ang_stop_samp = Ang_RewardLoc_ontrack(Ind_rewardloc);
end
if ~isnan(Sign_correct_test{1}(nl))
    ang_stop_test = Ang_RewardLoc_ontrack(Ind_rewardloc_test{1}(nl));
else
    ang_stop_test = Ang_RewardLoc_ontrack(Ind_rewardloc);
end
ang_vel_limit = 5/(Diam_inner/2);

%% Sample Lap data
loc_bins = scores_sample{nl,5};
vel_sample = scores_sample{nl,8}(1,:);
angvel_sample = scores_sample{nl,8}(4,:);
[~,ang_stop_bin] = min(abs(loc_bins-ang_stop_samp));
[~,ang_StartZone_dep_bin] = min(abs(loc_bins-Ang_StartZone_depart_ontrack));
[~,ang_StartZone_ari_bin] = min(abs(loc_bins-Ang_StartZone_arrive_ontrack));

% find actual locations bins
bin_actual_sample = scores_sample{nl,4};
loc_actual_sample_binned = loc_bins(bin_actual_sample);

% Sum up Pxn in each time bin between [current location, stop location]
t_start_sample = scores_sample{nl,6};
Pxn_sample = scores_sample{nl,3};
Pxn_sample_ahead_sum = nan(1,length(t_start_sample));
Pxn_sample_behind_sum = nan(1,length(t_start_sample));
Pxn_sample_all_sum = nan(1,length(t_start_sample));
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2));
% 1. swipe ahead probability
for i = ind
    temp = sort(Pxn_sample(bin_actual_sample(i):ang_stop_bin,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_sample_ahead_sum(i) = sum(temp(1:nbin_sum));
    end
end
% 2. swipe behind probability
for i = ind
    temp = sort(Pxn_sample(ang_StartZone_dep_bin:bin_actual_sample(i)-1,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_sample_behind_sum(i) = sum(temp(1:nbin_sum));
    end
end
% 3. swipe ahead+behind probability
Pxn_sample_all_sum(1,ind) = sum(Pxn_sample(ang_StartZone_dep_bin:ang_stop_bin,ind));

ind = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4));
% 1. swipe ahead probability
for i = ind
    temp = sort(Pxn_sample(bin_actual_sample(i):ang_StartZone_ari_bin,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_sample_ahead_sum(i) = sum(temp(1:nbin_sum));
    end
end
% 2. swipe behind probability
for i = ind
    temp = sort(Pxn_sample(ang_stop_bin+1:bin_actual_sample(i),i)-1,'descend');
    if(length(temp) >= nbin_sum)
        Pxn_sample_behind_sum(i) = sum(temp(1:nbin_sum));
    end
end
% 3. swipe ahead+behind probability
Pxn_sample_all_sum(1,ind) = sum(Pxn_sample(ang_stop_bin+1:ang_StartZone_ari_bin,ind));

% Sum up Pxn in each time bin between [current location+2bins, reward location]
t_start_sample = scores_sample{nl,6};
Pxn_sample = scores_sample{nl,3};
Pxn_sample_ahead_sum_2bin = nan(1,length(t_start_sample));
Pxn_sample_behind_sum_2bin = nan(1,length(t_start_sample));
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2));
% 1. swipe ahead probability
for i = ind
    temp = sort(Pxn_sample(bin_actual_sample(i)+2:ang_stop_bin,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_sample_ahead_sum_2bin(i) = sum(temp(1:nbin_sum));
    end
end
% 2. swipe behind probability
for i = ind
    temp = sort(Pxn_sample(ang_StartZone_dep_bin:bin_actual_sample(i)-2,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_sample_behind_sum_2bin(i) = sum(temp(1:nbin_sum));
    end
end

ind = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4));
% 1. swipe ahead probability
for i = ind
    temp = sort(Pxn_sample(bin_actual_sample(i)+2:ang_StartZone_ari_bin,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_sample_ahead_sum_2bin(i) = sum(temp(1:nbin_sum));
    end
end
% 2. swipe behind probability
for i = ind
    temp = sort(Pxn_sample(ang_stop_bin:bin_actual_sample(i)-2,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_sample_behind_sum_2bin(i) = sum(temp(1:nbin_sum));
    end
end

%% group sample lap data
% Before getting reward
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2));
Pxn_sample_ahead_prestop_mean = nanmean(Pxn_sample_ahead_sum(ind));
Pxn_sample_behind_prestop_mean = nanmean(Pxn_sample_behind_sum(ind));
angvel_sample_prestop_mean = nanmean(angvel_sample(ind));

% Before getting reward + running speed > 5cm/s
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2) &...
    angvel_sample >= ang_vel_limit);
Pxn_sample_ahead_vel_prestop_mean = nanmean(Pxn_sample_ahead_sum(ind));
Pxn_sample_behind_vel_prestop_mean = nanmean(Pxn_sample_behind_sum(ind));
angvel_sample_vel_prestop_mean = nanmean(angvel_sample(ind));

% Before getting reward + running speed > 5cm/s
% + approaching reward location (5 locations away from reward location):
if ~isnan(Sign_correct_sample{1}(nl))
    ang_reward_approach_sample = Ang_RewardLoc_ontrack(max(Ind_rewardloc_sample{1}(nl)-ind_approach,1));
else
    ang_reward_approach_sample = Ang_RewardLoc_ontrack(max(Ind_rewardloc-ind_approach,1));
end
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2) &...
    angvel_sample >= ang_vel_limit &...
    loc_actual_sample_binned' >= ang_reward_approach_sample);
ind0 = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2) & ...
    loc_actual_sample_binned' <= ang_reward_approach_sample);  % find the last point around -5 location
if ~isempty(ind0)
    ind0 = ind0(end);
    ind = ind(ind>=ind0);
end
Pxn_sample_ahead_approach_vel_prestop_mean = nanmean(Pxn_sample_ahead_sum(ind));
Pxn_sample_behind_approach_vel_prestop_mean = nanmean(Pxn_sample_behind_sum(ind));
Pxn_sample_all_approach_vel_prestop_mean = nanmean(Pxn_sample_all_sum(ind));
angvel_sample_approach_vel_prestop_mean = nanmean(angvel_sample(ind));

% Before getting reward + running speed > 5cm/s
% + approaching reward location (5 locations away from reward location):
% + 2bins ahead of the rat
if ~isnan(Sign_correct_sample{1}(nl))
    ang_reward_approach_sample = Ang_RewardLoc_ontrack(max(Ind_rewardloc_sample{1}(nl)-ind_approach,1));
else
    ang_reward_approach_sample = Ang_RewardLoc_ontrack(max(Ind_rewardloc-ind_approach,1));
end
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2) &...
    angvel_sample >= ang_vel_limit &...
    loc_actual_sample_binned' >= ang_reward_approach_sample);
ind0 = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2) & ...
    loc_actual_sample_binned' <= ang_reward_approach_sample);  % find the last point around -5 location
if ~isempty(ind0)
    ind0 = ind0(end);
    ind = ind(ind>=ind0);
end
Pxn_sample_ahead_approach_vel_prestop_2bin_mean = nanmean(Pxn_sample_ahead_sum_2bin(ind));
Pxn_sample_behind_approach_vel_prestop_2bin_mean = nanmean(Pxn_sample_behind_sum_2bin(ind));
angvel_sample_approach_vel_prestop_2bin_mean = nanmean(angvel_sample(ind));


% After getting reward
ind = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4));
Pxn_sample_ahead_poststop_mean = nanmean(Pxn_sample_ahead_sum(ind));
Pxn_sample_behind_poststop_mean = nanmean(Pxn_sample_behind_sum(ind));
angvel_sample_poststop_mean = nanmean(angvel_sample(ind));

% After getting reward + running speed > 5cm/s
ind = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4) &...
    angvel_sample >= ang_vel_limit);
Pxn_sample_ahead_vel_poststop_mean = nanmean(Pxn_sample_ahead_sum(ind));
Pxn_sample_behind_vel_poststop_mean = nanmean(Pxn_sample_behind_sum(ind));
angvel_sample_vel_poststop_mean = nanmean(angvel_sample(ind));

% After getting reward + running speed > 5cm/s
% + approaching reward location (5 locations away from reward location):
if ~isnan(Sign_correct_sample{1}(nl))
    ang_reward_approach_sample = Ang_RewardLoc_ontrack(min(Ind_rewardloc_sample{1}(nl)+ind_approach,N_loc));
else
    ang_reward_approach_sample = Ang_RewardLoc_ontrack(min(Ind_rewardloc+ind_approach,N_loc));
end
ind = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4) &...
    angvel_sample >= ang_vel_limit &...
    loc_actual_sample_binned' <= ang_reward_approach_sample);
ind0 = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4) &...
    loc_actual_sample_binned' >= ang_reward_approach_sample);
if ~isempty(ind0)
    ind0 = ind0(1);
    ind = ind(ind<=ind0);
end
Pxn_sample_ahead_approach_vel_poststop_mean = nanmean(Pxn_sample_ahead_sum(ind));
Pxn_sample_behind_approach_vel_poststop_mean = nanmean(Pxn_sample_behind_sum(ind));
Pxn_sample_all_approach_vel_poststop_mean = nanmean(Pxn_sample_all_sum(ind));
angvel_sample_approach_vel_poststop_mean = nanmean(angvel_sample(ind));

% After getting reward + running speed > 5cm/s
% + approaching reward location (5 locations away from reward location):
% + 2bins ahead of the rat
if ~isnan(Sign_correct_sample{1}(nl))
    ang_reward_approach_sample = Ang_RewardLoc_ontrack(min(Ind_rewardloc_sample{1}(nl)+ind_approach,N_loc));
else
    ang_reward_approach_sample = Ang_RewardLoc_ontrack(min(Ind_rewardloc+ind_approach,N_loc));
end
ind = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4) &...
    angvel_sample >= ang_vel_limit &...
    loc_actual_sample_binned' <= ang_reward_approach_sample);
ind0 = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4) &...
    loc_actual_sample_binned' >= ang_reward_approach_sample);
if ~isempty(ind0)
    ind0 = ind0(1);
    ind = ind(ind<=ind0);
end
Pxn_sample_ahead_approach_vel_poststop_2bin_mean = nanmean(Pxn_sample_ahead_sum_2bin(ind));
Pxn_sample_behind_approach_vel_poststop_2bin_mean = nanmean(Pxn_sample_behind_sum_2bin(ind));
angvel_sample_approach_vel_poststop_2bin_mean = nanmean(angvel_sample(ind));

%% Test Lap data
loc_bins = scores_test{nl,5};
vel_test = scores_test{nl,8}(1,:);
angvel_test = scores_test{nl,8}(4,:);
[~,ang_stop_bin] = min(abs(loc_bins-ang_stop_test));
[~,ang_StartZone_dep_bin] = min(abs(loc_bins-Ang_StartZone_depart_ontrack));
[~,ang_StartZone_ari_bin] = min(abs(loc_bins-Ang_StartZone_arrive_ontrack));

% find actual locations bins
bin_actual_test = scores_test{nl,4};
loc_actual_test_binned = loc_bins(bin_actual_test);

% Sum up Pxn in each time bin between [current location, reward location]
t_start_test = scores_test{nl,6};
Pxn_test = scores_test{nl,3};
Pxn_test_ahead_sum = nan(1,length(t_start_test));
Pxn_test_behind_sum = nan(1,length(t_start_test));
Pxn_test_all_sum = nan(1,length(t_start_test));
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2));
% 1. swipe ahead probability
for i = ind
    temp = sort(Pxn_test(bin_actual_test(i):ang_stop_bin,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_test_ahead_sum(i) = sum(temp(1:nbin_sum));
    end
end
% 2. swipe behind probability
for i = ind
    temp = sort(Pxn_test(ang_StartZone_dep_bin:bin_actual_test(i)-1,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_test_behind_sum(i) = sum(temp(1:nbin_sum));
    end
end
% 3. swipe ahead+behind probability
Pxn_test_all_sum(1,ind) = sum(Pxn_test(ang_StartZone_dep_bin:ang_stop_bin,ind));

ind = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4));
% 1. swipe ahead probability
for i = ind
    temp = sort(Pxn_test(bin_actual_test(i):ang_StartZone_ari_bin,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_test_ahead_sum(i) = sum(temp(1:nbin_sum));
    end
end
% 2. swipe behind probability
for i = ind
    temp = sort(Pxn_test(ang_stop_bin+1:bin_actual_test(i)-1,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_test_behind_sum(i) = sum(temp(1:nbin_sum));
    end
end
% 3. swipe ahead+behind probability
Pxn_test_all_sum(1,ind) = sum(Pxn_test(ang_stop_bin+1:ang_StartZone_ari_bin,ind));

% Sum up Pxn in each time bin between [current location+2bins, reward location]
t_start_test = scores_test{nl,6};
Pxn_test = scores_test{nl,3};
Pxn_test_ahead_sum_2bin = nan(1,length(t_start_test));
Pxn_test_behind_sum_2bin = nan(1,length(t_start_test));
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2));
% 1. swipe ahead probability
for i = ind
    temp = sort(Pxn_test(bin_actual_test(i)+2:ang_stop_bin,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_test_ahead_sum_2bin(i) = sum(temp(1:nbin_sum));
    end
end
% 2. swipe behind probability
for i = ind
    temp = sort(Pxn_test(ang_StartZone_dep_bin:bin_actual_test(i)-2,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_test_behind_sum_2bin(i) = sum(temp(1:nbin_sum));
    end
end

ind = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4));
% 1. swipe ahead probability
for i = ind
    temp = sort(Pxn_test(bin_actual_test(i)+2:ang_StartZone_ari_bin,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_test_ahead_sum_2bin(i) = sum(temp(1:nbin_sum));
    end
end
% 2. swipe behind probability
for i = ind
    temp = sort(Pxn_test(ang_stop_bin:bin_actual_test(i)-2,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_test_behind_sum_2bin(i) = sum(temp(1:nbin_sum));
    end
end

%% group test lap data
% Before getting reward
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2));
Pxn_test_ahead_prestop_mean = nanmean(Pxn_test_ahead_sum(ind));
Pxn_test_behind_prestop_mean = nanmean(Pxn_test_behind_sum(ind));
angvel_test_prestop_mean = nanmean(angvel_test(ind));

% Before getting reward + running speed > 5cm/s
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2) &...
    angvel_test >= ang_vel_limit);
Pxn_test_ahead_vel_prestop_mean = nanmean(Pxn_test_ahead_sum(ind));
Pxn_test_behind_vel_prestop_mean = nanmean(Pxn_test_behind_sum(ind));
angvel_test_vel_prestop_mean = nanmean(angvel_test(ind));

% Before getting reward + running speed > 5cm/s
% + approaching reward location (5 locations away from reward location):
if ~isnan(Sign_correct_test{1}(nl))
    ang_reward_approach_test = Ang_RewardLoc_ontrack(max(Ind_rewardloc_test{1}(nl)-ind_approach,1));
else
    ang_reward_approach_test = Ang_RewardLoc_ontrack(max(Ind_rewardloc-ind_approach,1));
end
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2) &...
    angvel_test >= ang_vel_limit &...
    loc_actual_test_binned' >= ang_reward_approach_test);
ind0 = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2) & ...
    loc_actual_test_binned' <= ang_reward_approach_test);  % find the last point around -5 location
if ~isempty(ind0)
    ind0 = ind0(end);
    ind = ind(ind>=ind0);
end
Pxn_test_ahead_approach_vel_prestop_mean = nanmean(Pxn_test_ahead_sum(ind));
Pxn_test_behind_approach_vel_prestop_mean = nanmean(Pxn_test_behind_sum(ind));
Pxn_test_all_approach_vel_prestop_mean = nanmean(Pxn_test_all_sum(ind));
angvel_test_approach_vel_prestop_mean = nanmean(angvel_test(ind));

% Before getting reward + running speed > 5cm/s
% + approaching reward location (5 locations away from reward location):
% + 2bins ahead of the rat
if ~isnan(Sign_correct_test{1}(nl))
    ang_reward_approach_test = Ang_RewardLoc_ontrack(max(Ind_rewardloc_test{1}(nl)-ind_approach,1));
else
    ang_reward_approach_test = Ang_RewardLoc_ontrack(max(Ind_rewardloc-ind_approach,1));
end
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2) &...
    angvel_test >= ang_vel_limit &...
    loc_actual_test_binned' >= ang_reward_approach_test);
ind0 = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2) & ...
    loc_actual_test_binned' <= ang_reward_approach_test);  % find the last point around -5 location
if ~isempty(ind0)
    ind0 = ind0(end);
    ind = ind(ind>=ind0);
end
Pxn_test_ahead_approach_vel_prestop_2bin_mean = nanmean(Pxn_test_ahead_sum_2bin(ind));
Pxn_test_behind_approach_vel_prestop_2bin_mean = nanmean(Pxn_test_behind_sum_2bin(ind));
angvel_test_approach_vel_prestop_2bin_mean = nanmean(angvel_test(ind));


% After getting reward
ind = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4));
Pxn_test_ahead_poststop_mean = nanmean(Pxn_test_ahead_sum(ind));
Pxn_test_behind_poststop_mean = nanmean(Pxn_test_behind_sum(ind));
angvel_test_poststop_mean = nanmean(angvel_test(ind));

% After getting reward + running speed > 5cm/s
ind = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4) &...
    angvel_test >= ang_vel_limit);
Pxn_test_ahead_vel_poststop_mean = nanmean(Pxn_test_ahead_sum(ind));
Pxn_test_behind_vel_poststop_mean = nanmean(Pxn_test_behind_sum(ind));
angvel_test_vel_poststop_mean = nanmean(angvel_test(ind));

% After getting reward + running speed > 5cm/s
% + approaching reward location (5 locations away from reward location):
if ~isnan(Sign_correct_test{1}(nl))
    ang_reward_approach_test = Ang_RewardLoc_ontrack(min(Ind_rewardloc_test{1}(nl)+ind_approach,N_loc));
else
    ang_reward_approach_test = Ang_RewardLoc_ontrack(min(Ind_rewardloc+ind_approach,N_loc));
end
ind = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4) &...
    angvel_test >= ang_vel_limit &...
    loc_actual_test_binned' <= ang_reward_approach_test);
ind0 = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4) &...
    loc_actual_test_binned' >= ang_reward_approach_test);
if ~isempty(ind0)
    ind0 = ind0(1);
    ind = ind(ind<=ind0);
end
Pxn_test_ahead_approach_vel_poststop_mean = nanmean(Pxn_test_ahead_sum(ind));
Pxn_test_behind_approach_vel_poststop_mean = nanmean(Pxn_test_behind_sum(ind));
Pxn_test_all_approach_vel_poststop_mean = nanmean(Pxn_test_all_sum(ind));
angvel_test_approach_vel_poststop_mean = nanmean(angvel_test(ind));

% After getting reward + running speed > 5cm/s
% + approaching reward location (5 locations away from reward location):
% + 2bins ahead of the rat
if ~isnan(Sign_correct_test{1}(nl))
    ang_reward_approach_test = Ang_RewardLoc_ontrack(min(Ind_rewardloc_test{1}(nl)+ind_approach,N_loc));
else
    ang_reward_approach_test = Ang_RewardLoc_ontrack(min(Ind_rewardloc+ind_approach,N_loc));
end
ind = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4) &...
    angvel_test >= ang_vel_limit &...
    loc_actual_test_binned' <= ang_reward_approach_test);
ind0 = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4) &...
    loc_actual_test_binned' >= ang_reward_approach_test);
if ~isempty(ind0)
    ind0 = ind0(1);
    ind = ind(ind<=ind0);
end
Pxn_test_ahead_approach_vel_poststop_2bin_mean = nanmean(Pxn_test_ahead_sum_2bin(ind));
Pxn_test_behind_approach_vel_poststop_2bin_mean = nanmean(Pxn_test_behind_sum_2bin(ind));
angvel_test_approach_vel_poststop_2bin_mean = nanmean(angvel_test(ind));


%% Save data
sign_sample = 1;
sign_test = 2;

data0_prestop = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_prestop_mean,Pxn_sample_behind_prestop_mean,angvel_sample_prestop_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_prestop_mean,Pxn_test_behind_prestop_mean,angvel_test_prestop_mean];
data0_vel_prestop = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_vel_prestop_mean,Pxn_sample_behind_vel_prestop_mean,angvel_sample_vel_prestop_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_vel_prestop_mean,Pxn_test_behind_vel_prestop_mean,angvel_test_vel_prestop_mean];
data0_approach_vel_prestop = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_approach_vel_prestop_mean,Pxn_sample_behind_approach_vel_prestop_mean,angvel_sample_approach_vel_prestop_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_approach_vel_prestop_mean,Pxn_test_behind_approach_vel_prestop_mean,angvel_test_approach_vel_prestop_mean,...
    Pxn_sample_all_approach_vel_prestop_mean,Pxn_test_all_approach_vel_prestop_mean];
data0_approach_vel_prestop_2bin = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_approach_vel_prestop_2bin_mean,Pxn_sample_behind_approach_vel_prestop_2bin_mean,angvel_sample_approach_vel_prestop_2bin_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_approach_vel_prestop_2bin_mean,Pxn_test_behind_approach_vel_prestop_2bin_mean,angvel_test_approach_vel_prestop_2bin_mean];


data0_poststop = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_poststop_mean,Pxn_sample_behind_poststop_mean,angvel_sample_poststop_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_poststop_mean,Pxn_test_behind_poststop_mean,angvel_test_poststop_mean];
data0_vel_poststop = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_vel_poststop_mean,Pxn_sample_behind_vel_poststop_mean,angvel_sample_vel_poststop_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_vel_poststop_mean,Pxn_test_behind_vel_poststop_mean,angvel_test_vel_poststop_mean];
data0_approach_vel_poststop = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_approach_vel_poststop_mean,Pxn_sample_behind_approach_vel_poststop_mean,angvel_sample_approach_vel_poststop_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_approach_vel_poststop_mean,Pxn_test_behind_approach_vel_poststop_mean,angvel_test_approach_vel_poststop_mean,...
    Pxn_sample_all_approach_vel_poststop_mean,Pxn_test_all_approach_vel_poststop_mean];
data0_approach_vel_poststop_2bin = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_approach_vel_poststop_2bin_mean,Pxn_sample_behind_approach_vel_poststop_2bin_mean,angvel_sample_approach_vel_poststop_2bin_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_approach_vel_poststop_2bin_mean,Pxn_test_behind_approach_vel_poststop_2bin_mean,angvel_test_approach_vel_poststop_2bin_mean];