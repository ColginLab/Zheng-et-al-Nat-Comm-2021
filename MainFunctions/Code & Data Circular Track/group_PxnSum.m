%% Group Sum Pxn of continuous decodings in each sample-test trial
% edited on 07/25/2017
% Modified on 05/07/2019, add measurements 4 and 8

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
ang_reward = Ang_RewardLoc_ontrack(Ind_rewardloc);
ang_vel_limit = 5/(Diam_inner/2);

%% Sample Lap data
loc_bins = scores_sample{nl,5};
vel_sample = scores_sample{nl,8}(1,:);
angvel_sample = scores_sample{nl,8}(4,:);
[~,ang_reward_bin] = min(abs(loc_bins-ang_reward));
[~,ang_StartZone_dep_bin] = min(abs(loc_bins-Ang_StartZone_depart_ontrack));
[~,ang_StartZone_ari_bin] = min(abs(loc_bins-Ang_StartZone_arrive_ontrack));

% find actual locations bins
bin_actual_sample = scores_sample{nl,4};
loc_actual_sample_binned = loc_bins(bin_actual_sample);

% Sum up Pxn in each time bin between [current location, reward location]
t_start_sample = scores_sample{nl,6};
Pxn_sample = scores_sample{nl,3};
Pxn_sample_ahead_sum = nan(1,length(t_start_sample));
Pxn_sample_behind_sum = nan(1,length(t_start_sample));
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2));
% 1. swipe ahead probability
for i = ind
    temp = sort(Pxn_sample(bin_actual_sample(i):ang_reward_bin,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_sample_ahead_sum(i) = sum(temp(1:nbin_sum));
    end
end
% 2. swipe behind probability
for i = ind
    temp = sort(Pxn_sample(ang_StartZone_dep_bin:bin_actual_sample(i),i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_sample_behind_sum(i) = sum(temp(1:nbin_sum));
    end
end

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
    temp = sort(Pxn_sample(ang_reward_bin:bin_actual_sample(i),i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_sample_behind_sum(i) = sum(temp(1:nbin_sum));
    end
end

% Sum up Pxn in each time bin between [current location+2bins, reward location]
t_start_sample = scores_sample{nl,6};
Pxn_sample = scores_sample{nl,3};
Pxn_sample_ahead_sum_2bin = nan(1,length(t_start_sample));
Pxn_sample_behind_sum_2bin = nan(1,length(t_start_sample));
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2));
% 1. swipe ahead probability
for i = ind
    temp = sort(Pxn_sample(bin_actual_sample(i)+2:ang_reward_bin,i),'descend');
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
    temp = sort(Pxn_sample(ang_reward_bin:bin_actual_sample(i)-2,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_sample_behind_sum_2bin(i) = sum(temp(1:nbin_sum));
    end
end

%% group sample lap data
% Before getting reward
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2));
Pxn_sample_ahead_prereward_mean = nanmean(Pxn_sample_ahead_sum(ind));
Pxn_sample_behind_prereward_mean = nanmean(Pxn_sample_behind_sum(ind));
angvel_sample_prereward_mean = nanmean(angvel_sample(ind));

% Before getting reward + running speed > 5cm/s
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2) &...
    angvel_sample >= ang_vel_limit);
Pxn_sample_ahead_vel_prereward_mean = nanmean(Pxn_sample_ahead_sum(ind));
Pxn_sample_behind_vel_prereward_mean = nanmean(Pxn_sample_behind_sum(ind));
angvel_sample_vel_prereward_mean = nanmean(angvel_sample(ind));

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
Pxn_sample_ahead_approach_vel_prereward_mean = nanmean(Pxn_sample_ahead_sum(ind));
Pxn_sample_behind_approach_vel_prereward_mean = nanmean(Pxn_sample_behind_sum(ind));
angvel_sample_approach_vel_prereward_mean = nanmean(angvel_sample(ind));

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
Pxn_sample_ahead_approach_vel_prereward_2bin_mean = nanmean(Pxn_sample_ahead_sum_2bin(ind));
Pxn_sample_behind_approach_vel_prereward_2bin_mean = nanmean(Pxn_sample_behind_sum_2bin(ind));
angvel_sample_approach_vel_prereward_2bin_mean = nanmean(angvel_sample(ind));


% After getting reward
ind = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4));
Pxn_sample_ahead_postreward_mean = nanmean(Pxn_sample_ahead_sum(ind));
Pxn_sample_behind_postreward_mean = nanmean(Pxn_sample_behind_sum(ind));
angvel_sample_postreward_mean = nanmean(angvel_sample(ind));

% After getting reward + running speed > 5cm/s
ind = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4) &...
    angvel_sample >= ang_vel_limit);
Pxn_sample_ahead_vel_postreward_mean = nanmean(Pxn_sample_ahead_sum(ind));
Pxn_sample_behind_vel_postreward_mean = nanmean(Pxn_sample_behind_sum(ind));
angvel_sample_vel_postreward_mean = nanmean(angvel_sample(ind));

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
Pxn_sample_ahead_approach_vel_postreward_mean = nanmean(Pxn_sample_ahead_sum(ind));
Pxn_sample_behind_approach_vel_postreward_mean = nanmean(Pxn_sample_behind_sum(ind));
angvel_sample_approach_vel_postreward_mean = nanmean(angvel_sample(ind));

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
Pxn_sample_ahead_approach_vel_postreward_2bin_mean = nanmean(Pxn_sample_ahead_sum_2bin(ind));
Pxn_sample_behind_approach_vel_postreward_2bin_mean = nanmean(Pxn_sample_behind_sum_2bin(ind));
angvel_sample_approach_vel_postreward_2bin_mean = nanmean(angvel_sample(ind));

%% Test Lap data
loc_bins = scores_test{nl,5};
vel_test = scores_test{nl,8}(1,:);
angvel_test = scores_test{nl,8}(4,:);
[~,ang_reward_bin] = min(abs(loc_bins-ang_reward));
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
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2));
% 1. swipe ahead probability
for i = ind
    temp = sort(Pxn_test(bin_actual_test(i):ang_reward_bin,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_test_ahead_sum(i) = sum(temp(1:nbin_sum));
    end
end
% 2. swipe behind probability
for i = ind
    temp = sort(Pxn_test(ang_StartZone_dep_bin:bin_actual_test(i),i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_test_behind_sum(i) = sum(temp(1:nbin_sum));
    end
end

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
    temp = sort(Pxn_test(ang_reward_bin:bin_actual_test(i),i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_test_behind_sum(i) = sum(temp(1:nbin_sum));
    end
end

% Sum up Pxn in each time bin between [current location+2bins, reward location]
t_start_test = scores_test{nl,6};
Pxn_test = scores_test{nl,3};
Pxn_test_ahead_sum_2bin = nan(1,length(t_start_test));
Pxn_test_behind_sum_2bin = nan(1,length(t_start_test));
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2));
% 1. swipe ahead probability
for i = ind
    temp = sort(Pxn_test(bin_actual_test(i)+2:ang_reward_bin,i),'descend');
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
    temp = sort(Pxn_test(ang_reward_bin:bin_actual_test(i)-2,i),'descend');
    if(length(temp) >= nbin_sum)
        Pxn_test_behind_sum_2bin(i) = sum(temp(1:nbin_sum));
    end
end

%% group test lap data
% Before getting reward
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2));
Pxn_test_ahead_prereward_mean = nanmean(Pxn_test_ahead_sum(ind));
Pxn_test_behind_prereward_mean = nanmean(Pxn_test_behind_sum(ind));
angvel_test_prereward_mean = nanmean(angvel_test(ind));

% Before getting reward + running speed > 5cm/s
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2) &...
    angvel_test >= ang_vel_limit);
Pxn_test_ahead_vel_prereward_mean = nanmean(Pxn_test_ahead_sum(ind));
Pxn_test_behind_vel_prereward_mean = nanmean(Pxn_test_behind_sum(ind));
angvel_test_vel_prereward_mean = nanmean(angvel_test(ind));

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
Pxn_test_ahead_approach_vel_prereward_mean = nanmean(Pxn_test_ahead_sum(ind));
Pxn_test_behind_approach_vel_prereward_mean = nanmean(Pxn_test_behind_sum(ind));
angvel_test_approach_vel_prereward_mean = nanmean(angvel_test(ind));

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
Pxn_test_ahead_approach_vel_prereward_2bin_mean = nanmean(Pxn_test_ahead_sum_2bin(ind));
Pxn_test_behind_approach_vel_prereward_2bin_mean = nanmean(Pxn_test_behind_sum_2bin(ind));
angvel_test_approach_vel_prereward_2bin_mean = nanmean(angvel_test(ind));


% After getting reward
ind = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4));
Pxn_test_ahead_postreward_mean = nanmean(Pxn_test_ahead_sum(ind));
Pxn_test_behind_postreward_mean = nanmean(Pxn_test_behind_sum(ind));
angvel_test_postreward_mean = nanmean(angvel_test(ind));

% After getting reward + running speed > 5cm/s
ind = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4) &...
    angvel_test >= ang_vel_limit);
Pxn_test_ahead_vel_postreward_mean = nanmean(Pxn_test_ahead_sum(ind));
Pxn_test_behind_vel_postreward_mean = nanmean(Pxn_test_behind_sum(ind));
angvel_test_vel_postreward_mean = nanmean(angvel_test(ind));

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
Pxn_test_ahead_approach_vel_postreward_mean = nanmean(Pxn_test_ahead_sum(ind));
Pxn_test_behind_approach_vel_postreward_mean = nanmean(Pxn_test_behind_sum(ind));
angvel_test_approach_vel_postreward_mean = nanmean(angvel_test(ind));

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
Pxn_test_ahead_approach_vel_postreward_2bin_mean = nanmean(Pxn_test_ahead_sum_2bin(ind));
Pxn_test_behind_approach_vel_postreward_2bin_mean = nanmean(Pxn_test_behind_sum_2bin(ind));
angvel_test_approach_vel_postreward_2bin_mean = nanmean(angvel_test(ind));


%% Save data
sign_sample = 1;
sign_test = 2;

data0_prereward = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_prereward_mean,Pxn_sample_behind_prereward_mean,angvel_sample_prereward_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_prereward_mean,Pxn_test_behind_prereward_mean,angvel_test_prereward_mean];
data0_vel_prereward = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_vel_prereward_mean,Pxn_sample_behind_vel_prereward_mean,angvel_sample_vel_prereward_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_vel_prereward_mean,Pxn_test_behind_vel_prereward_mean,angvel_test_vel_prereward_mean];
data0_approach_vel_prereward = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_approach_vel_prereward_mean,Pxn_sample_behind_approach_vel_prereward_mean,angvel_sample_approach_vel_prereward_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_approach_vel_prereward_mean,Pxn_test_behind_approach_vel_prereward_mean,angvel_test_approach_vel_prereward_mean];
data0_approach_vel_prereward_2bin = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_approach_vel_prereward_2bin_mean,Pxn_sample_behind_approach_vel_prereward_2bin_mean,angvel_sample_approach_vel_prereward_2bin_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_approach_vel_prereward_2bin_mean,Pxn_test_behind_approach_vel_prereward_2bin_mean,angvel_test_approach_vel_prereward_2bin_mean];


data0_postreward = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_postreward_mean,Pxn_sample_behind_postreward_mean,angvel_sample_postreward_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_postreward_mean,Pxn_test_behind_postreward_mean,angvel_test_postreward_mean];
data0_vel_postreward = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_vel_postreward_mean,Pxn_sample_behind_vel_postreward_mean,angvel_sample_vel_postreward_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_vel_postreward_mean,Pxn_test_behind_vel_postreward_mean,angvel_test_vel_postreward_mean];
data0_approach_vel_postreward = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_approach_vel_postreward_mean,Pxn_sample_behind_approach_vel_postreward_mean,angvel_sample_approach_vel_postreward_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_approach_vel_postreward_mean,Pxn_test_behind_approach_vel_postreward_mean,angvel_test_approach_vel_postreward_mean];
data0_approach_vel_postreward_2bin = [ns, nl,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    Pxn_sample_ahead_approach_vel_postreward_2bin_mean,Pxn_sample_behind_approach_vel_postreward_2bin_mean,angvel_sample_approach_vel_postreward_2bin_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    Pxn_test_ahead_approach_vel_postreward_2bin_mean,Pxn_test_behind_approach_vel_postreward_2bin_mean,angvel_test_approach_vel_postreward_2bin_mean];