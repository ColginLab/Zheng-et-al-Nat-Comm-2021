% v1: make segments evenly from Loc1 to stop location
speed_seg = nan(1,n_seg);

ts0_start = Ts0{1,1}(nl,1)/1000000;
ts0_stop = Ts0{1,1}(nl,2)/1000000;
data0 = data_angle0{nl};

loc1 = Ang_RewardLoc_ontrack(1);  % count from location1

% method1
% loc_stop = Ang_RewardLoc_ontrack(Ind_rewardloc0{1}(nl)); % count to stop location

% method2
[~,ind] = min(abs(data0(:,1)-ts0_stop));
loc_stop = data0(ind,2); % count to stop location

loc_seg = loc1:(loc_stop-loc1)/n_seg:loc_stop;

if isnan(ts0_stop)
    ind = find(data0(:,2)<= loc_stop & data0(:,1)>=ts0_start);
    ind = ind(end);
    ts0_stop = data0(ind,1);
end

for ind_seg = 1:n_seg
    loc_seg_i1 = loc_seg(ind_seg);
    loc_seg_i2 = loc_seg(ind_seg+1);
    
    % find the time point when the rat was in each segment
    ind = find(data0(:,2)>=loc_seg_i1 & data0(:,2)<=loc_seg_i2 & data0(:,1)<=ts0_stop & data0(:,1)>=ts0_start);
    
    % average the running speed within each segments
    speed_seg(ind_seg) = mean(data0(ind,3));
end