%% Group averaged Pxn, check if it can predict correct vs. error trials
% edited on 07/31/2017
% modified on 10/6/2017 : exclude the p(x|n) at rat's current location

%% 
%disp(strcat('ntrial = ',num2str(nl)))
ang_reward = Ang_RewardLoc_ontrack(Ind_rewardloc);
ang_vel_limit = 5/(Diam_inner/2);

stop_loc_id = Ind_rewardloc_test{1}(nl);
if ~isnan(stop_loc_id)
    stop_minus1_loc = Ang_RewardLoc_ontrack(stop_loc_id-1);
else
    stop_minus1_loc = Ang_RewardLoc_ontrack(Ind_rewardloc-1);
end


% find each reward location zone
loc_bins = scores_sample{nl,5};
df_reward = mean(diff(Ang_RewardLoc_ontrack));
ang_reward_zone_border = mean([Ang_RewardLoc_ontrack(1:end-1),Ang_RewardLoc_ontrack(2:end)],2);
ang_reward_zone_border = [Ang_RewardLoc_ontrack(1)-df_reward;ang_reward_zone_border];
ang_reward_zone_border = [ang_reward_zone_border;Ang_RewardLoc_ontrack(end)+df_reward];
for i = 1:length(ang_reward_zone_border)
    [~,ang_reward_zone_border_bin(i,1)] = min(abs(ang_reward_zone_border(i)-loc_bins));
end
ang_reward_zone = [ang_reward_zone_border(1:end-1),ang_reward_zone_border(2:end)];
ang_reward_zone_bin = [ang_reward_zone_border_bin(1:end-1),ang_reward_zone_border_bin(2:end)];

%% Sample Lap data
vel_sample = scores_sample{nl,8}(1,:);
angvel_sample = scores_sample{nl,8}(4,:);
[~,ang_reward_bin] = min(abs(loc_bins-ang_reward));
[~,ang_StartZone_depart_bin] = min(abs(loc_bins-Ang_StartZone_depart_ontrack));
[~,ang_StartZone_arrive_bin] = min(abs(loc_bins-Ang_StartZone_arrive_ontrack));

% find actual locations bins
bin_actual_sample = scores_sample{nl,4};
loc_actual_sample_binned = loc_bins(bin_actual_sample);

% v1: Pxn_sample_2
% extract max 10 Pxn in each time bin between [current location, box_arrive]
t_start_sample = scores_sample{nl,6};
Pxn_sample = scores_sample{nl,3};
Pxn_sample_2 = nan(size(Pxn_sample));
ind1 = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2));
ind2 = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4));
ind = unique([ind1,ind2]);
for i = ind
    if ang_StartZone_arrive_bin-bin_actual_sample(i)+1 >= nbin_sum
        temp = nan(size(Pxn_sample,1),1);
        temp(bin_actual_sample(i):ang_StartZone_arrive_bin,1) = Pxn_sample(bin_actual_sample(i):ang_StartZone_arrive_bin,i);
        [B,I] = sort(temp,'descend');
        ind_temp = find(~isnan(B));
        if ~isempty(ind_temp)
            I = I(ind_temp(1):ind_temp(1)+nbin_sum-1);
            I = sort(I);
            Pxn_sample_2(I,i) = Pxn_sample(I,i);
        end
    end
end

% v2: Pxn_sample_3
% extract max 10 Pxn in each time bin between [box_depart, box_arrive]
% exclude the +-1 bins around the rat's current bin
t_start_sample = scores_sample{nl,6};
Pxn_sample = scores_sample{nl,3};
Pxn_sample_3 = nan(size(Pxn_sample));
ind1 = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2));
ind2 = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4));
ind = unique([ind1,ind2]);
for i = ind
    temp = nan(size(Pxn_sample,1),1);
    % exclude the +-1 bins around the rat's current bin
    ind_temp = setdiff(ang_StartZone_depart_bin:ang_StartZone_arrive_bin,bin_actual_sample(i)-1:bin_actual_sample(i)+1);
    
    temp(ind_temp,1) = Pxn_sample(ind_temp,i);
    [B,I] = sort(temp,'descend');
    ind_temp = find(~isnan(B));
    if ~isempty(ind_temp)
        I = I(ind_temp(1):ind_temp(1)+nbin_sum-1);
        I = sort(I);
        Pxn_sample_3(I,i) = Pxn_sample(I,i);
    end
end

% v3: Pxn_sample_4
% extract max 10 Pxn in each time bin between [current location+2bins, box_arrive]
% (exclude the +1 bins around the rat's current bin)
t_start_sample = scores_sample{nl,6};
Pxn_sample = scores_sample{nl,3};
Pxn_sample_4 = nan(size(Pxn_sample));
ind1 = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2));
ind2 = find(t_start_sample >= ts_start_stop_sample_nl(3) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(4));
ind = unique([ind1,ind2]);
for i = ind
    if ang_StartZone_arrive_bin-(bin_actual_sample(i)+2)+1 >= nbin_sum
        temp = nan(size(Pxn_sample,1),1);
        temp(bin_actual_sample(i)+2:ang_StartZone_arrive_bin,1) = Pxn_sample(bin_actual_sample(i)+2:ang_StartZone_arrive_bin,i);
        [B,I] = sort(temp,'descend');
        ind_temp = find(~isnan(B));
        if ~isempty(ind_temp)
            I = I(ind_temp(1):ind_temp(1)+nbin_sum-1);
            I = sort(I);
            Pxn_sample_4(I,i) = Pxn_sample(I,i);
        end
    end
end

%% group sample lap data

% !!!!!! choose data to use  !!!!!!
if strcmp(code_version,'v1')
    Pxn_sample_touse = Pxn_sample_2;  % v1
elseif strcmp(code_version,'v2')
    Pxn_sample_touse = Pxn_sample_3;  % v2
elseif strcmp(code_version,'v3')
    Pxn_sample_touse = Pxn_sample_4;  % v3
end

% Before getting stop-1 location
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2) & ...
    loc_actual_sample_binned' < stop_minus1_loc);
Pxn_sample_prereward_mean = nanmean(Pxn_sample_touse(:,ind),2);
Pxn_sample_prereward_binsum = nan(N_loc,1);
for i = 1:N_loc
    temp = Pxn_sample_prereward_mean(ang_reward_zone_bin(i,1):ang_reward_zone_bin(i,2));
    if sum(~isnan(temp)) > 0
        Pxn_sample_prereward_binsum(i,1) = nansum(temp);
    end
end
angvel_sample_prereward_mean = nanmean(angvel_sample(ind));

% Before getting stop-1 location
% + running speed > 5cm/s
ind = find(t_start_sample >= ts_start_stop_sample_nl(1) & ...
    t_start_sample+dt <= ts_start_stop_sample_nl(2) & ...
    loc_actual_sample_binned' < stop_minus1_loc &...
    angvel_sample >= ang_vel_limit);
Pxn_sample_vel_prereward_mean = nanmean(Pxn_sample_touse(:,ind),2);
Pxn_sample_vel_prereward_binsum = nan(N_loc,1);
for i = 1:N_loc
    temp = Pxn_sample_vel_prereward_mean(ang_reward_zone_bin(i,1):ang_reward_zone_bin(i,2));
    if sum(~isnan(temp)) > 0
        Pxn_sample_vel_prereward_binsum(i,1) = nansum(temp);
    end
end
angvel_sample_vel_prereward_mean = nanmean(angvel_sample(ind));

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
Pxn_sample_approach_vel_prereward_mean = nanmean(Pxn_sample_touse(:,ind),2);
Pxn_sample_approach_vel_prereward_binsum = nan(N_loc,1);
for i = 1:N_loc
    temp = Pxn_sample_approach_vel_prereward_mean(ang_reward_zone_bin(i,1):ang_reward_zone_bin(i,2));
    if sum(~isnan(temp)) > 0
        Pxn_sample_approach_vel_prereward_binsum(i,1) = nansum(temp);
    end
end
angvel_sample_approach_vel_prereward_mean = nanmean(angvel_sample(ind));


%% Test Lap data
vel_test = scores_test{nl,8}(1,:);
angvel_test = scores_test{nl,8}(4,:);
[~,ang_reward_bin] = min(abs(loc_bins-ang_reward));
[~,ang_StartZone_arrive_bin] = min(abs(loc_bins-Ang_StartZone_arrive_ontrack));

% find actual locations bins
bin_actual_test = scores_test{nl,4};
loc_actual_test_binned = loc_bins(bin_actual_test);

% v1: Pxn_test_2
% extract max 10 Pxn in each time bin between [current location, box_arrive]
t_start_test = scores_test{nl,6};
Pxn_test = scores_test{nl,3};
Pxn_test_2 = nan(size(Pxn_test));
ind1 = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2));
ind2 = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4));
ind = unique([ind1,ind2]);
for i = ind
    if ang_StartZone_arrive_bin-bin_actual_test(i)+1 >= nbin_sum
        temp = nan(size(Pxn_test,1),1);
        temp(bin_actual_test(i):ang_StartZone_arrive_bin,1) = Pxn_test(bin_actual_test(i):ang_StartZone_arrive_bin,i);
        [B,I] = sort(temp,'descend');
        ind_temp = find(~isnan(B));
        if ~isempty(ind_temp)
            I = I(ind_temp(1):ind_temp(1)+nbin_sum-1);
            I = sort(I);
            Pxn_test_2(I,i) = Pxn_test(I,i);
        end
    end
end

% v2: Pxn_test_3
% extract max 10 Pxn in each time bin between [box_depart, box_arrive]
% exclude the +-1 bins around the rat's current bin
t_start_test = scores_test{nl,6};
Pxn_test = scores_test{nl,3};
Pxn_test_3 = nan(size(Pxn_test));
ind1 = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2));
ind2 = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4));
ind = unique([ind1,ind2]);
for i = ind
    temp = nan(size(Pxn_test,1),1);
    % exclude the +-1 bins around the rat's current bin
    ind_temp = setdiff(ang_StartZone_depart_bin:ang_StartZone_arrive_bin,bin_actual_test(i)-1:bin_actual_test(i)+1);
    
    temp(ind_temp,1) = Pxn_test(ind_temp,i);
    [B,I] = sort(temp,'descend');
    ind_temp = find(~isnan(B));
    if ~isempty(ind_temp)
        I = I(ind_temp(1):ind_temp(1)+nbin_sum-1);
        I = sort(I);
        Pxn_test_3(I,i) = Pxn_test(I,i);
    end
end

% v3: Pxn_test_4
% extract max 10 Pxn in each time bin between [current location+2bins, box_arrive]
% (exclude the +1 bins around the rat's current bin)
t_start_test = scores_test{nl,6};
Pxn_test = scores_test{nl,3};
Pxn_test_4 = nan(size(Pxn_test));
ind1 = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2));
ind2 = find(t_start_test >= ts_start_stop_test_nl(3) & ...
    t_start_test+dt <= ts_start_stop_test_nl(4));
ind = unique([ind1,ind2]);
for i = ind
    if ang_StartZone_arrive_bin-(bin_actual_test(i)+2)+1 >= nbin_sum
        temp = nan(size(Pxn_test,1),1);
        temp(bin_actual_test(i)+2:ang_StartZone_arrive_bin,1) = Pxn_test(bin_actual_test(i)+2:ang_StartZone_arrive_bin,i);
        [B,I] = sort(temp,'descend');
        ind_temp = find(~isnan(B));
        if ~isempty(ind_temp)
            I = I(ind_temp(1):ind_temp(1)+nbin_sum-1);
            I = sort(I);
            Pxn_test_4(I,i) = Pxn_test(I,i);
        end
    end
end
%% group test lap data

% !!!!!! choose data to use  !!!!!!
if strcmp(code_version,'v1')
    Pxn_test_touse = Pxn_test_2;  % v1
elseif strcmp(code_version,'v2')
    Pxn_test_touse = Pxn_test_3;  % v2
elseif strcmp(code_version,'v3')
    Pxn_test_touse = Pxn_test_4;  % v3
end

% Before getting stop-1 location
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2) & ...
    loc_actual_test_binned' < stop_minus1_loc);
Pxn_test_prereward_mean = nanmean(Pxn_test_touse(:,ind),2);
Pxn_test_prereward_binsum = nan(N_loc,1);
for i = 1:N_loc
    temp = Pxn_test_prereward_mean(ang_reward_zone_bin(i,1):ang_reward_zone_bin(i,2));
    if sum(~isnan(temp)) > 0
        Pxn_test_prereward_binsum(i,1) = nansum(temp);
    end
end
angvel_test_prereward_mean = nanmean(angvel_test(ind));

% Before getting stop-1 location
% + running speed > 5cm/s
ind = find(t_start_test >= ts_start_stop_test_nl(1) & ...
    t_start_test+dt <= ts_start_stop_test_nl(2) & ...
    loc_actual_test_binned' < stop_minus1_loc &...
    angvel_test >= ang_vel_limit);
Pxn_test_vel_prereward_mean = nanmean(Pxn_test_touse(:,ind),2);
Pxn_test_vel_prereward_binsum = nan(N_loc,1);
for i = 1:N_loc
    temp = Pxn_test_vel_prereward_mean(ang_reward_zone_bin(i,1):ang_reward_zone_bin(i,2));
    if sum(~isnan(temp)) > 0
        Pxn_test_vel_prereward_binsum(i,1) = nansum(temp);
    end
end
angvel_test_vel_prereward_mean = nanmean(angvel_test(ind));

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
Pxn_test_approach_vel_prereward_mean = nanmean(Pxn_test_touse(:,ind),2);
Pxn_test_approach_vel_prereward_binsum = nan(N_loc,1);
for i = 1:N_loc
    temp = Pxn_test_approach_vel_prereward_mean(ang_reward_zone_bin(i,1):ang_reward_zone_bin(i,2));
    if sum(~isnan(temp)) > 0
        Pxn_test_approach_vel_prereward_binsum(i,1) = nansum(temp);
    end
end
angvel_test_approach_vel_prereward_mean = nanmean(angvel_test(ind));


%% Save data
sign_sample = 1;
sign_test = 2;

data0_prereward = [ns,nl,Ind_rewardloc,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    angvel_sample_prereward_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    angvel_test_prereward_mean];
data0_vel_prereward = [ns,nl,Ind_rewardloc,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    angvel_sample_vel_prereward_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    angvel_test_vel_prereward_mean];
data0_approach_vel_prereward = [ns,nl,Ind_rewardloc,...
    sign_sample, Sign_correct_sample{1}(nl),Ind_rewardloc_sample{1}(nl)...
    angvel_sample_approach_vel_prereward_mean,...
    sign_test, Sign_correct_test{1}(nl),Ind_rewardloc_test{1}(nl)...
    angvel_test_approach_vel_prereward_mean];
