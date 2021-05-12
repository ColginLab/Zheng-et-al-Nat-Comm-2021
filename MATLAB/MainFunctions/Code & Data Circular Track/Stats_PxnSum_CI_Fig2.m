function Stats_PxnSum_CI_Fig2(AnalysisDir, FiguresDir)

%% Do stats for group data of sum Pxn
% in sample/test trial
% Edited on 07/30/2017
% Update on 04/13/2021
file_output1 = 'group_PxnSum_5LocAway_25cells_20201217_5bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
file_output1 = strcat(AnalysisDir,file_output1);
load(file_output1)

%file_output1 = 'group_PxnSum_5LocAway_25cells_20201217_5bins_v3_ds.mat'; %4 rats, method(1-2),use downsampled place cells
%load(file_output1)

userat = [139,150,149,148];

ind1 = ismember(ns_RatID(:,2),userat);
ind1 = find(ind1>0);

ind2 = max(ns_RatID(ind1,1));
para = [4,3,1,1;...
    4,3,1,2;...
    4,3,1,3;...
    4,3,2,1;...
    4,3,2,2;...
    4,3,2,3;...
    3,3,3,1;...
    3,3,3,2;...
    3,3,3,3];
FontSize = 20;
FigLabels = {'Fig2a', 'Fig2b', 'Fig2c', 'Fig2d', 'Fig2e', 'Fig2f', 'Fig2g', 'Fig2h', 'Fig2i'};
ct = 1;
for i = 1:size(para,1)
    para_method = para(i,1);
    para_reward = para(i,2);
    para_ah_bh = para(i,3);
    para_crterr = para(i,4);
    
    if para_method == 1
        % Method 1. from starting zone to reward location
        Group_Pxn_prereward = group_Pxn_prereward(ind1,:);
        Group_Pxn_postreward = group_Pxn_postreward(ind1,:);
        Group_Pxn_prestop = group_Pxn_prestop(ind1,:);
        Group_Pxn_poststop = group_Pxn_poststop(ind1,:);
    elseif para_method == 2
        % Method 2. from starting zone to reward location, filter running speed above 5cm/s
        Group_Pxn_prereward = group_Pxn_vel_prereward(ind1,:);
        Group_Pxn_postreward = group_Pxn_vel_postreward(ind1,:);
        Group_Pxn_prestop = group_Pxn_vel_prestop(ind1,:);
        Group_Pxn_poststop = group_Pxn_vel_poststop(ind1,:);
    elseif para_method == 3
        % Method 3. from 5 locations earlier to reward location, filter running speed above 5cm/s
        Group_Pxn_prereward = group_Pxn_approach_vel_prereward(ind1,:);
        Group_Pxn_postreward = group_Pxn_approach_vel_postreward(ind1,:);
        Group_Pxn_prestop = group_Pxn_approach_vel_prestop(ind1,:);
        Group_Pxn_poststop = group_Pxn_approach_vel_poststop(ind1,:);
    elseif para_method == 4
        % Method 4. from 5 locations earlier to reward location, filter running speed above 5cm/s, 2bins ahead of the rat
        Group_Pxn_prereward = group_Pxn_approach_vel_prereward_2bin(ind1,:);
        Group_Pxn_postreward = group_Pxn_approach_vel_postreward_2bin(ind1,:);
        Group_Pxn_prestop = group_Pxn_approach_vel_prestop_2bin(ind1,:);
        Group_Pxn_poststop = group_Pxn_approach_vel_poststop_2bin(ind1,:);
    end
    %% choose to do pre-reward or post-reward?
    if para_reward == 1
        % Pre-reward
        Group_Pxn0 = Group_Pxn_prereward;
        title_reward = 'Pre-reward';
    elseif para_reward == 2
        % Post-reward
        Group_Pxn0 = Group_Pxn_postreward;
        title_reward = 'Post-reward';
    elseif para_reward == 3
        % Pre-stop
        Group_Pxn0 = Group_Pxn_prestop;
        title_reward = 'Pre-stop';
    elseif para_reward == 4
        % Post-stop
        Group_Pxn0 = Group_Pxn_poststop;
        title_reward = 'Post-stop';
    end
    StatsData_LineGraph = [];
    StatsData_BarGraph = [];
    %% choose to do probability of ahead or behind the rat?
    if para_ah_bh == 1
        i_samp = 6; i_test = 12; dirType = 'ahead'; % ahead probability
    elseif para_ah_bh == 2
        i_samp = 7; i_test = 13; dirType = 'behind'; % behind probability
    elseif para_ah_bh == 3
        i_samp = 15; i_test = 16; dirType = 'ahead+behind'; % ahead probability   % only for Method3 data
    end
    %%  Group data:
    % all samples and all tests;
    % all correct and all error laps
    
    isession = max(Group_Pxn0(:,1));
    pxn_crt_err = nan(isession,4);
    pxn_all = nan(isession,2);
    angvel_crt_err = nan(isession,4);
    angvel_all = nan(isession,2);
    for ns = 1:isession
        ind = find(Group_Pxn0(:,1) == ns);
        if ~isempty(ind)
            group_pxn = Group_Pxn0(ind,:);
            ind_nan = find(isnan(mean(group_pxn,2)));
            group_pxn(ind_nan,:) = [];
            ind_correct = find(group_pxn(:,10)==1);
            ind_error = find(group_pxn(:,10)==0);
            if ~isempty(ind_correct) & ~isempty(ind_error)
                pxn_crt_err(ns,1) = mean(group_pxn(ind_correct,i_samp)); % correct sample
                pxn_crt_err(ns,2) = mean(group_pxn(ind_correct,i_test)); % correct test
                pxn_crt_err(ns,3) = mean(group_pxn(ind_error,i_samp)); % error sample
                pxn_crt_err(ns,4) = mean(group_pxn(ind_error,i_test)); % error test
                
                angvel_crt_err(ns,1) = mean(group_pxn(ind_correct,8)); % correct sample
                angvel_crt_err(ns,2) = mean(group_pxn(ind_correct,14)); % correct test
                angvel_crt_err(ns,3) = mean(group_pxn(ind_error,8)); % error sample
                angvel_crt_err(ns,4) = mean(group_pxn(ind_error,14)); % error test
            end
            pxn_all(ns,1) = mean(group_pxn([ind_correct;ind_error],i_samp)); % all samples
            pxn_all(ns,2) = mean(group_pxn([ind_correct;ind_error],i_test)); % all tests
            
            angvel_all(ns,1) = mean(group_pxn([ind_correct;ind_error],8)); % all samples
            angvel_all(ns,2) = mean(group_pxn([ind_correct;ind_error],14)); % all tests
        end
    end
    
    %%  Group data:
    % samples 1-8, or 1-4, or 5-8
    % sample1 and the first correct test
    
    nl_limit_early = [1,4];
    nl_limit_late = [5,8];
    
    isession = max(Group_Pxn0(:,1));
    
    pxn_allsamps = nan(isession,8);
    angvel_allsamps = nan(isession,8);
    pxn_alltests = nan(isession,8);
    angvel_alltests = nan(isession,8);
    
    pxn_crtsamps = nan(isession,8);
    angvel_crtsamps = nan(isession,8);
    pxn_crttests = nan(isession,8);
    angvel_crttests = nan(isession,8);
    
    pxn_errsamps = nan(isession,8);
    angvel_errsamps = nan(isession,8);
    pxn_errtests = nan(isession,8);
    angvel_errtests = nan(isession,8);
    
    pxn_samp1_crt = nan(isession,4);
    pxn_crt1_err1 = nan(isession,4);
    angvel_samp1_crt = nan(isession,4);
    
    for ns = 1:isession
        ind = find(Group_Pxn0(:,1) == ns);
        if ~isempty(ind)
            group_pxn = Group_Pxn0(ind,:);
            
            % group all sample and tests
            pxn_allsamps(ns,:) = group_pxn(:,i_samp)';
            angvel_allsamps(ns,:) = group_pxn(:,8)';
            
            pxn_alltests(ns,:) = group_pxn(:,i_test)';
            angvel_alltests(ns,:) = group_pxn(:,14)';
            
            % group corrent sample and tests
            ind_correct = find(group_pxn(:,10)==1);
            pxn_crtsamps(ns,ind_correct) = group_pxn(ind_correct,i_samp)';
            angvel_crtsamps(ns,ind_correct) = group_pxn(ind_correct,8)';
            
            pxn_crttests(ns,ind_correct) = group_pxn(ind_correct,i_test)';
            angvel_crttests(ns,ind_correct) = group_pxn(ind_correct,14)';
            
            % group error sample and tests
            ind_error = find(group_pxn(:,10)==0);
            pxn_errsamps(ns,ind_error) = group_pxn(ind_error,i_samp)';
            angvel_errsamps(ns,ind_error) = group_pxn(ind_error,8)';
            
            pxn_errtests(ns,ind_error) = group_pxn(ind_error,i_test)';
            angvel_errtests(ns,ind_error) = group_pxn(ind_error,14)';
            
            % group sample1, test1, early crt testN, late crt testN
            ind_early = find(ind_correct >= nl_limit_early(1) & ind_correct <= nl_limit_early(2));
            ind_correct_early = ind_correct(ind_early);
            ind_late = find(ind_correct >= nl_limit_late(1) & ind_correct <= nl_limit_late(2));
            ind_correct_late = ind_correct(ind_late);
            
            ind_early = find(ind_error >= nl_limit_early(1) & ind_error <= nl_limit_early(2));
            ind_error_early = ind_error(ind_early);
            ind_late = find(ind_error >= nl_limit_late(1) & ind_error <= nl_limit_late(2));
            ind_error_late = ind_error(ind_late);
            
            if ~isempty(ind_correct_early) && ~isempty(ind_correct_late)
                pxn_samp1_crt(ns,1) = group_pxn(1,i_samp); % sample1
                pxn_samp1_crt(ns,2) = group_pxn(ind_correct_early(1),i_test); % 1st correct test
                pxn_samp1_crt(ns,3) = group_pxn(ind_correct_early(end),i_test); % last early correct test
                pxn_samp1_crt(ns,4) = group_pxn(ind_correct_late(end),i_test); % last late correct test
                
                angvel_samp1_crt(ns,1) = group_pxn(1,8); % sample1
                angvel_samp1_crt(ns,2) = group_pxn(ind_correct_early(1),14); % 1st correct test
                angvel_samp1_crt(ns,3) = group_pxn(ind_correct_early(end),14); % last correct test
                angvel_samp1_crt(ns,4) = group_pxn(ind_correct_late(end),14); % last correct test
            end
            
            if ~isempty(ind_correct_early) && ~isempty(ind_correct_late) && ~isempty(ind_error_early) && ~isempty(ind_error_late)
                pxn_crt1_err1(ns,1) = group_pxn(ind_correct_early(1),i_test); % 1st eraly correct test
                pxn_crt1_err1(ns,2) = group_pxn(ind_correct_late(end),i_test); % last late correct test
                pxn_crt1_err1(ns,3) = group_pxn(ind_error_early(1),i_test); % 1st eraly err test
                pxn_crt1_err1(ns,4) = group_pxn(ind_error_late(end),i_test); % last late err test
            end
        end
    end
    
    
    %% Plot
    % samples 1-8;
    % sample1 and the first correct test
    
    ind = find(isnan(mean(pxn_allsamps+pxn_alltests,2)));
    pxn_allsamps(ind,:) = [];
    angvel_allsamps(ind,:) = [];
    pxn_alltests(ind,:) = [];
    angvel_alltests(ind,:) = [];
    
    if para_crterr == 1
        data1 = pxn_allsamps; data2 = pxn_alltests; title0 = 'All trials';
    elseif para_crterr == 2
        data1 = pxn_crtsamps; data2 = pxn_crttests; title0 = 'Correct trials';
    elseif para_crterr == 3
        data1 = pxn_errsamps; data2 = pxn_errtests; title0 = 'Error trials';
    end
    data0 = angvel_all;
    
    
    % Scatter plot of Pxn over 8 laps for all/crt/err trials sample and tests
    ind = repmat(1:size(data1,2),size(data1,1),1);
    data1_dot = [reshape(ind,size(ind,1)*size(ind,2),1),reshape(data1,size(data1,1)*size(data1,2),1)];
    data1_dot(isnan(sum(data1_dot,2)),:) = [];
    
    ind = repmat(1:size(data2,2),size(data2,1),1);
    data2_dot = [reshape(ind,size(ind,1)*size(ind,2),1),reshape(data2,size(data2,1)*size(data2,2),1)];
    data2_dot(isnan(sum(data2_dot,2)),:) = [];
    
    [b1,bint,r,rint,stats1] = regress(data1_dot(:,2),[ones(size(data1_dot,1),1),data1_dot(:,1)]);
    [b2,bint,r,rint,stats2] = regress(data2_dot(:,2),[ones(size(data2_dot,1),1),data2_dot(:,1)]);
    x0 = 0.5:1:size(data1,2)+0.5;
    y1 = ones(length(x0),1)*b1(1)+x0'*b1(2);
    y2 = ones(length(x0),1)*b2(1)+x0'*b2(2);
    
    [rho1,pval1] = corr(data1_dot(:,1),data1_dot(:,2));
    [rho2,pval2] = corr(data2_dot(:,1),data2_dot(:,2));
    [rho0,pval0] = corr([data1_dot(:,1);data2_dot(:,1)],[data1_dot(:,2);data2_dot(:,2)]);
    rho_pval = [rho1,pval1;rho2,pval2;rho0,pval0];
    
    sample_color = [255,165,0]./255;
    test_color = [111,57,214]./255;
    
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    scatter(data1_dot(:,1)-0.1,data1_dot(:,2),25,'filled','MarkerEdgeColor',sample_color,'MarkerFaceColor',sample_color)
    scatter(data2_dot(:,1)+0.1,data2_dot(:,2),25,'filled','MarkerEdgeColor',test_color,'MarkerFaceColor',test_color);
    plot(x0,y1,'color',sample_color,'linewidth',2);
    plot(x0,y2,'color',test_color,'linewidth',2);
    text(1,0.1,strcat('cc = ',num2str(roundn(rho1,-3)),', p = ',num2str(roundn(pval1,-3))),'Color',sample_color,'FontSize',FontSize);
    text(1,0.05,strcat('cc = ',num2str(roundn(rho2,-3)),', p = ',num2str(roundn(pval2,-3))),'Color',test_color,'FontSize',FontSize);
    hold off
    set(gca,'XTick',[1:size(data1,2)]);
    xlim([0,size(data1,2)+1])
    if para_method == 4
        ytick0 = 0:0.1:0.7;
    elseif para_method == 3
        ytick0 = 0:0.1:1;
    end
    ylim(ytick0([1,end]))
    set(gca,'YTick',ytick0);
    set(gca,'fontsize',20);
    xlabel('Trials');
    if strcmp(title_reward, 'Pre-stop') || strcmp(title_reward, 'Pre-reward')
        if strcmp(dirType,'ahead') % ahead probability
            if strcmp(title_reward, 'Pre-stop')
                ylabel([{'Summation of P(x|n)'},{'from current location to stop location'}]);
            elseif strcmp(title_reward, 'Pre-reward')
                ylabel([{'Summation of P(x|n)'},{'from current location to reward location'}]);
            end
        elseif strcmp(dirType,'behind') % behind probability
            ylabel([{'Summation of P(x|n)'},{'from start location to current location'}]);
        elseif strcmp(dirType,'ahead+behind') % behind probability
            if strcmp(title_reward, 'Pre-stop')
                ylabel([{'Summation of P(x|n)'},{'from start location to stop location'}]);
            elseif strcmp(title_reward, 'Pre-reward')
                ylabel([{'Summation of P(x|n)'},{'from start location to reward location'}]);
            end
        end
    end
    
    if strcmp(title_reward, 'Post-stop') || strcmp(title_reward, 'Post-reward')
        if strcmp(dirType,'ahead') % ahead probability
            ylabel([{'Summation of P(x|n)'},{'from current location to end location'}]);
        elseif strcmp(dirType,'behind') % behind probability
            if strcmp(title_reward, 'Post-stop')
                ylabel([{'Summation of P(x|n)'},{'from stop location to current location'}]);
            elseif strcmp(title_reward, 'Post-reward')
                ylabel([{'Summation of P(x|n)'},{'from reward location to current location'}]);
            end
        elseif strcmp(dirType,'ahead+behind') % behind probability
            if strcmp(title_reward, 'Post-stop')
                ylabel([{'Summation of P(x|n)'},{'from stop location to end location'}]);
            elseif strcmp(title_reward, 'Post-reward')
                ylabel([{'Summation of P(x|n)'},{'from reward location to end location'}]);
            end
        end
    end
    title(title0)
    saveas(h1, [FiguresDir,'\',FigLabels{ct}],'fig')
    saveas(h1, [FiguresDir,'\',FigLabels{ct}],'epsc')
    %saveas([Outdir,'\SeqExample_',TrialID{tt},'_',RatID,'_',DateID],'epsc')
    ct = ct + 1;
    % Data for SPSS
    if strcmp(title0,'Correct trials')
        data_all = [];
        
        len = size(data1_dot,1);
        data_temp = [ones(len,1),ones(len,1),data1_dot];
        data_temp(:,5) = data_temp(:,1).*data_temp(:,3);
        data_all = [data_all;data_temp];
        
        len = size(data2_dot,1);
        data_temp = [ones(len,1),ones(len,1)*2,data2_dot];
        data_temp(:,5) = data_temp(:,1).*data_temp(:,3);
        data_all = [data_all;data_temp];
    elseif strcmp(title0,'Error trials')
        len = size(data1_dot,1);
        data_temp = [ones(len,1)*2,ones(len,1),data1_dot];
        data_temp(:,5) = data_temp(:,1).*data_temp(:,3);
        data_all = [data_all;data_temp];
        
        len = size(data2_dot,1);
        data_temp = [ones(len,1)*2,ones(len,1)*2,data2_dot];
        data_temp(:,5) = data_temp(:,1).*data_temp(:,3);
        data_all = [data_all;data_temp];
    end
end
end