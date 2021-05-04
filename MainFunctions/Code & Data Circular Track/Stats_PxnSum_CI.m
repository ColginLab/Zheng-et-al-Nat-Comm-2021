%% Do stats for group data of sum Pxn
% in sample/test trial
% Edited on 07/30/2017
% Update on 1/10/2018

userat = [139,150,149,148];

ind1 = ismember(ns_RatID(:,2),userat);
ind1 = find(ind1>0);

ind2 = max(ns_RatID(ind1,1));

%% Method 1. from starting zone to reward location
% Group_Pxn_prereward = group_Pxn_prereward(ind1,:);
% Group_Pxn_postreward = group_Pxn_postreward(ind1,:);
% Group_Pxn_prestop = group_Pxn_prestop(ind1,:);
% Group_Pxn_poststop = group_Pxn_poststop(ind1,:);

%% Method 2. from starting zone to reward location, filter running speed above 5cm/s
% Group_Pxn_prereward = group_Pxn_vel_prereward(ind1,:);
% Group_Pxn_postreward = group_Pxn_vel_postreward(ind1,:);
% Group_Pxn_prestop = group_Pxn_vel_prestop(ind1,:);
% Group_Pxn_poststop = group_Pxn_vel_poststop(ind1,:);

%% Method 3. from 5 locations earlier to reward location, filter running speed above 5cm/s
% Group_Pxn_prereward = group_Pxn_approach_vel_prereward(ind1,:);
% Group_Pxn_postreward = group_Pxn_approach_vel_postreward(ind1,:);
% Group_Pxn_prestop = group_Pxn_approach_vel_prestop(ind1,:);
% Group_Pxn_poststop = group_Pxn_approach_vel_poststop(ind1,:);

%% Method 4. from 5 locations earlier to reward location, filter running speed above 5cm/s, 2bins ahead of the rat
Group_Pxn_prereward = group_Pxn_approach_vel_prereward_2bin(ind1,:);
Group_Pxn_postreward = group_Pxn_approach_vel_postreward_2bin(ind1,:);
Group_Pxn_prestop = group_Pxn_approach_vel_prestop_2bin(ind1,:);
Group_Pxn_poststop = group_Pxn_approach_vel_poststop_2bin(ind1,:);

%% choose to do pre-reward or post-reward?
% Pre-reward
% Group_Pxn0 = Group_Pxn_prereward;
% title_reward = 'Pre-reward';

% Post-reward
% Group_Pxn0 = Group_Pxn_postreward;
% title_reward = 'Post-reward';

% Pre-stop
Group_Pxn0 = Group_Pxn_prestop;
title_reward = 'Pre-stop';

% Post-stop
% Group_Pxn0 = Group_Pxn_poststop;
% title_stop = 'Post-stop';

StatsData_LineGraph = [];
StatsData_BarGraph = [];
%% choose to do probability of ahead or behind the rat?
i_samp = 6; i_test = 12; dirType = 'ahead'; % ahead probability
% i_samp = 7; i_test = 13; dirType = 'behind'; % behind probability

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

%% Stats: Compare all sample vs. all tests
ind = find(isnan(mean(pxn_all,2)));
pxn_all(ind,:) = [];
angvel_all(ind,:) = [];

[h,p_pxn_all,ci,stats] = ttest(pxn_all(:,1),pxn_all(:,2));p_pxn_all
[h,p_angvel_all,ci,stats] = ttest(angvel_all(:,1),angvel_all(:,2));p_angvel_all

p_pxn_all = signrank(pxn_all(:,1),pxn_all(:,2));p_pxn_all
p_angvel_all = signrank(angvel_all(:,1),angvel_all(:,2));p_angvel_all

data0 = pxn_all;
ylabel0 = 'Pxn';

% data0 = angvel_all;
% ylabel0 = 'Angular speed (rad/s)';

sample_color = [255,165,0]./255;
test_color = [111,57,214]./255;
ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
data0_mean = mean(data0);
data0_sem = std(data0)./sqrt(size(data0,1));
hold on
h1 = errorbar(1,data0_mean(1),data0_sem(1),'s-','color',sample_color,'LineWidth',2);
h2 = errorbar(4,data0_mean(2),data0_sem(2),'s-','color',test_color,'LineWidth',2);
errorbar_tick(h1, 0);
errorbar_tick(h2, 0);
plot([2,3],data0,'k');
plot(2,data0(:,1),'o','color',sample_color);
plot(3,data0(:,2),'o','color',test_color);
hold off
set(gca,'XTick',[1,4]);
xlim([0,5])
set(gca,'XTickLabel',{'Sample trials','Test trials'});
set(gca,'fontsize',20);


ylabel(ylabel0);
% ylim([-0.3,0.2])
% set(gca,'YTick',-0.3:0.1:0.2);
title('Pre-reward')
title('Pre-reward,Vel >5cm/s')
title('Pre-reward,Vel >5cm/s and approaching stop location')


%% Stats: Compare correct vs. error trials
ind = find(isnan(mean(pxn_crt_err,2)));
pxn_crt_err(ind,:) = [];
angvel_crt_err(ind,:) = [];


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

% data1 = pxn_allsamps; data2 = pxn_alltests;
data1 = pxn_crtsamps; data2 = pxn_crttests;
% data1 = pxn_errsamps; data2 = pxn_errtests;
data0 = angvel_all;


% Scatter plot of Pxn over 8 laps for sample and tests
ind = repmat(1:size(data1,2),size(data1,1),1);
data1_dot = [reshape(ind,size(ind,1)*size(ind,2),1),reshape(data1,size(data1,1)*size(data1,2),1)];
data1_dot(isnan(sum(data1_dot,2)),:) = [];

ind = repmat(1:size(data2,2),size(data2,1),1);
data2_dot = [reshape(ind,size(ind,1)*size(ind,2),1),reshape(data2,size(data2,1)*size(data2,2),1)];
data2_dot(isnan(sum(data2_dot,2)),:) = [];

[b1,bint,r,rint,stats1] = regress(data1_dot(:,2),[ones(size(data1_dot,1),1),data1_dot(:,1)]);
[b2,bint,r,rint,stats2] = regress(data2_dot(:,2),[ones(size(data2_dot,1),1),data2_dot(:,1)]);
x0 = 1:size(data1,2);
y1 = ones(size(data1,2),1)*b1(1)+x0'*b1(2);
y2 = ones(size(data2,2),1)*b2(1)+x0'*b2(2);

[rho1,pval1] = corr(data1_dot(:,1),data1_dot(:,2));
[rho2,pval2] = corr(data2_dot(:,1),data2_dot(:,2));

sample_color = [255,165,0]./255;
test_color = [111,57,214]./255;
figure('Units','normalized','Position',[0 0 0.3 0.9]);
hold on
scatter(data1_dot(:,1),data1_dot(:,2),'MarkerEdgeColor',sample_color);
scatter(data2_dot(:,1),data2_dot(:,2),'filled','MarkerEdgeColor',test_color,'MarkerFaceColor',test_color);
% lh = legend('Samples','Tests','Location','SouthEast');
% set(lh,'Box','off')
plot(x0,y1,'color',sample_color);
plot(x0,y2,'color',test_color);
text(1,0.1,strcat('Samples r = ',num2str(roundn(rho1,-3))),'Color',sample_color,'FontSize',20);
text(1,0.05,strcat('Tests r = ',num2str(roundn(rho2,-3))),'Color',test_color,'FontSize',20);
hold off
set(gca,'XTick',[1:size(data1,2)]);
xlim([0,size(data1,2)+1])
ytick0 = 0:0.1:0.7;
ylim(ytick0([1,end]))
set(gca,'YTick',ytick0);
set(gca,'fontsize',20);
xlabel('Trials');
if strcmp(dirType,'ahead') % ahead probability
    ylabel([{'Summation of P(x|n)'},{'from current location to reward location'}]);
elseif strcmp(dirType,'behind') % behind probability
    ylabel([{'Summation of P(x|n)'},{'from start location to current location'}]);
end


% Plot Pxn over 8 laps for sample and tests
% Bootstrapped confidence interval parameters
s2=statset('bootci');
s2.UseSubstreams=true;
s2.Streams=RandStream('mlfg6331_64','Seed',1);
nBoot = 5000;

sample_color = [255,165,0]./255;
test_color = [111,57,214]./255;
ffa2 = figure('Units','normalized','Position',[0 0 0.4 0.6]);
data1_mean = nanmean(data1);
data1_CI = bootci(nBoot,{@nanmean,data1},'Options',s2);
data1_CI(1,:) = data1_mean-data1_CI(1,:);
data1_CI(2,:) = data1_CI(2,:)-data1_mean;
data1_sem = std(data1)./sqrt(size(data1,1));
data2_mean = nanmean(data2);
data2_CI = bootci(nBoot,{@nanmean,data2},'Options',s2);
data2_CI(1,:) = data2_mean-data2_CI(1,:);
data2_CI(2,:) = data2_CI(2,:)-data2_mean;
data2_sem = std(data2)./sqrt(size(data2,1));
reset(s2.Streams)
hold on
h1 = errorbar(1:length(data1_mean),data1_mean,data1_CI(1,:),data1_CI(2,:),'s-','color',sample_color,'LineWidth',2);
h2 = errorbar(1:length(data2_mean),data2_mean,data2_CI(1,:),data2_CI(2,:),'s-','color',test_color,'LineWidth',2);
% errorbar_tick(h1, 0);
% errorbar_tick(h2, 0);
lh = legend('Samples','Tests','Location','SouthEast');
set(lh,'Box','off')
hold off
set(gca,'XTick',[1:length(data1_mean)]);
xlim([0,length(data1_mean)+1])
ytick0 = 0:0.1:0.5;
ylim(ytick0([1,end]))
set(gca,'YTick',ytick0);
set(gca,'fontsize',20);
xlabel('Trials');
if strcmp(dirType,'ahead') % ahead probability
    ylabel([{'Summation of P(x|n)'},{'from current location to reward location'}]);
elseif strcmp(dirType,'behind') % behind probability
    ylabel([{'Summation of P(x|n)'},{'from start location to current location'}]);
end

% title(title_reward)
% title(strcat(title_reward,',Vel >5cm/s'))
% title(strcat(title_reward,',Vel >5cm/s and approaching stop location'))

% Stats Data
pxnSum = [data1, data2];
DirType = repmat({dirType},size(pxnSum,1),1);
statsData = [cell2table(DirType), array2table(pxnSum,'VariableNames',{'T1_samp', 'T2_samp', 'T3_samp', 'T4_samp', 'T5_samp','T6_samp', 'T7_samp', 'T8_samp',...
                                                                      'T1_test', 'T2_test', 'T3_test', 'T4_test', 'T5_test','T6_test', 'T7_test', 'T8_test'})];
StatsData_LineGraph = cat(1,StatsData_LineGraph,statsData);

pxn_samp1_crt_plot = pxn_samp1_crt;

angvel_samp1_crt_plot = angvel_samp1_crt;

ind = find(isnan(mean(pxn_samp1_crt_plot,2)));
pxn_samp1_crt_plot(ind,:) = [];
angvel_samp1_crt_plot(ind,:) = [];

[h,p_pxn_samp1_crt,ci,stats] = ttest(pxn_samp1_crt_plot(:,1),pxn_samp1_crt_plot(:,3));p_pxn_samp1_crt
[h,p_angvel_samp1_crt,ci,stats] = ttest(angvel_samp1_crt_plot(:,1),angvel_samp1_crt_plot(:,3));p_angvel_samp1_crt

p_pxn_samp1_crt = signrank(pxn_samp1_crt_plot(:,1),pxn_samp1_crt_plot(:,2));p_pxn_samp1_crt
p_angvel_samp1_crt = signrank(angvel_samp1_crt_plot(:,1),angvel_samp1_crt_plot(:,2));p_angvel_samp1_crt

data0 = pxn_samp1_crt_plot;
% data0 = angvel_samp1_crt_plot;


% plot Sample1/ E-crt test1/ L-crt testN
sample_color = [255,165,0]./255;
test_color = [111,57,214]./255;
ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
data0_mean = mean(data0);
data0_sem = std(data0)./sqrt(size(data0,1));
hold on
bar(1,data0_mean(1),'FaceColor',sample_color,'EdgeColor', sample_color);
bar(2,data0_mean(2),'FaceColor',test_color,'EdgeColor', test_color);
bar(3,data0_mean(4),'FaceColor',test_color,'EdgeColor', test_color);
h1 = errorbar(1,data0_mean(1),data0_sem(1),'color','k','LineWidth',2);
h2 = errorbar(2,data0_mean(2),data0_sem(2),'color','k','LineWidth',2);
h3 = errorbar(3,data0_mean(4),data0_sem(4),'color','k','LineWidth',2);
errorbar_tick(h1, 0);
errorbar_tick(h2, 0);
errorbar_tick(h3, 0);
% plot([1,2,3],data0,'ko');
hold off
set(gca,'XTick',[1,2,3]);
xlim([0.5,3.5])
ytick0 = 0:0.1:0.5;
ylim(ytick0([1,end]))
set(gca,'YTick',ytick0);
set(gca,'XTickLabel',{'Sample1','1st crt test','Last crt test'});
set(gca,'fontsize',20);

ylabel([{'Summation of P(x|n)'},{'from current location to reward location'}]);
% ylabel('Angular speed (rad/s)');

title(title_reward)
title(strcat(title_reward,',Vel >5cm/s'))
title(strcat(title_reward,',Vel >5cm/s and approaching stop location'))


% plot Sample1/ E-crt test1/ E-crt testN/ L-crt testN
sample_color = [255,165,0]./255;
test_color = [111,57,214]./255;
ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
data0_mean = mean(data0);
data0_sem = std(data0)./sqrt(size(data0,1));
hold on
bar(1,data0_mean(1),'FaceColor',sample_color,'EdgeColor', sample_color);
bar(2,data0_mean(2),'FaceColor',test_color,'EdgeColor', test_color);
bar(3,data0_mean(3),'FaceColor',test_color,'EdgeColor', test_color);
bar(4,data0_mean(4),'FaceColor',test_color,'EdgeColor', test_color);
h1 = errorbar(1,data0_mean(1),data0_sem(1),'color','k','LineWidth',2);
h2 = errorbar(2,data0_mean(2),data0_sem(2),'color','k','LineWidth',2);
h3 = errorbar(3,data0_mean(3),data0_sem(3),'color','k','LineWidth',2);
h4 = errorbar(4,data0_mean(4),data0_sem(4),'color','k','LineWidth',2);
errorbar_tick(h1, 0);
errorbar_tick(h2, 0);
errorbar_tick(h3, 0);
errorbar_tick(h4, 0);
% plot([1,2,3],data0,'ko');
hold off
set(gca,'XTick',[1,2,3,4]);
xlim([0.5,4.5])
ytick0 = 0:0.1:0.5;
ylim(ytick0([1,end]))
set(gca,'YTick',ytick0);
set(gca,'XTickLabel',{'Sample1','E-crt test1','E-crt testN','L-crt testN'});
set(gca,'fontsize',20);

ylabel([{'Summation of P(x|n)'},{'from current location to reward location'}]);
% ylabel('Angular speed (rad/s)');

title(title_reward)
title(strcat(title_reward,',Vel >5cm/s'))
title(strcat(title_reward,',Vel >5cm/s and approaching stop location'))


% plot E-crt test1/ L-crt testN/ E-err test1/ L-err testN
pxn_crt1_err1_plot = pxn_crt1_err1;
ind = find(isnan(mean(pxn_crt1_err1_plot,2)));
pxn_crt1_err1_plot(ind,:) = [];
[h,p_pxn_crt1_err1,ci,stats] = ttest(pxn_crt1_err1_plot(:,1),pxn_crt1_err1_plot(:,3));p_pxn_crt1_err1
p_pxn_crt1_err1 = signrank(pxn_crt1_err1_plot(:,1),pxn_crt1_err1_plot(:,2));p_pxn_crt1_err1

data0 = pxn_crt1_err1_plot;

test_color = [111,57,214]./255;
ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
data0_mean = mean(data0);
data0_sem = std(data0)./sqrt(size(data0,1));
hold on
bar(1,data0_mean(1),'FaceColor',test_color,'EdgeColor', test_color);
bar(2,data0_mean(2),'FaceColor',test_color,'EdgeColor', test_color);
bar(3,data0_mean(3),'FaceColor',test_color,'EdgeColor', test_color);
bar(4,data0_mean(4),'FaceColor',test_color,'EdgeColor', test_color);
h1 = errorbar(1,data0_mean(1),data0_sem(1),'color','k','LineWidth',2);
h2 = errorbar(2,data0_mean(2),data0_sem(2),'color','k','LineWidth',2);
h3 = errorbar(3,data0_mean(3),data0_sem(3),'color','k','LineWidth',2);
h4 = errorbar(4,data0_mean(4),data0_sem(4),'color','k','LineWidth',2);
errorbar_tick(h1, 0);
errorbar_tick(h2, 0);
errorbar_tick(h3, 0);
errorbar_tick(h4, 0);
% plot([1,2,3],data0,'ko');
hold off
set(gca,'XTick',[1,2,3,4]);
xlim([0.5,4.5])
ytick0 = 0:0.1:0.5;
ylim(ytick0([1,end]))
set(gca,'YTick',ytick0);
set(gca,'XTickLabel',{'E-crt test1','L-crt testN','E-err test1','L-err testN'});
set(gca,'fontsize',20);

ylabel([{'Summation of P(x|n)'},{'from current location to reward location'}]);
% ylabel('Angular speed (rad/s)');

title(title_reward)
title(strcat(title_reward,',Vel >5cm/s'))
title(strcat(title_reward,',Vel >5cm/s and approaching stop location'))



% For repeated measures ANOVA, switch col3 and col4 in data0
% so that we can do multiple comparison by using "Simple" compares to the "last value"
data0_new = data0(:,[1,2,4,3]);

% Stats Data
pxnSum = data0(:,[1 2 4]);
DirType = repmat({dirType},size(pxnSum,1),1);
statsData = [cell2table(DirType), array2table(pxnSum,'VariableNames',{'Sample1','Test_crt1','Test_crtN'})];
StatsData_BarGraph = cat(1,StatsData_BarGraph,statsData);


%% Group data: 
% tests error-(including err-1,err-2,...
% tests error+(including err+1,err+2,...

isession = max(Group_Pxn0(:,1));

pxn_allsamps = nan(isession,8);
angvel_allsamps = nan(isession,8);

pxn_alltests = nan(isession,8);
angvel_alltests = nan(isession,8);

pxn_err_mins_plus = nan(isession,2);
angvel_err_mins_plus = nan(isession,2);

for ns = 1:isession
    ind = find(Group_Pxn0(:,1) == ns);
    if ~isempty(ind)
        group_pxn = Group_Pxn0(ind,:);
        
        pxn_alltests(ns,:) = group_pxn(:,i_test)';
        angvel_alltests(ns,:) = group_pxn(:,14)';
        
        ind_err_mins = find(group_pxn(:,10)==0 & group_pxn(:,11)<group_pxn(:,5));
        ind_err_plus = find(group_pxn(:,10)==0 & group_pxn(:,11)>group_pxn(:,5));
        
        if ~isempty(ind_err_mins)
%             pxn_err_mins_plus(ns,1) = mean(group_pxn(ind_err_mins,i_test)); % all error- tests
%             angvel_err_mins_plus(ns,1) = mean(group_pxn(ind_err_mins,14)); % all error- tests
            
            pxn_err_mins_plus(ns,1) = group_pxn(ind_err_mins(end),i_test); % all error- tests
            angvel_err_mins_plus(ns,1) = group_pxn(ind_err_mins(end),14); % all error- tests
        end
        if ~isempty(ind_err_plus)
%             pxn_err_mins_plus(ns,2) = mean(group_pxn(ind_err_plus,i_test)); % all error+ tests
%             angvel_err_mins_plus(ns,2) = mean(group_pxn(ind_err_plus,14)); % all error+ tests
            
            pxn_err_mins_plus(ns,2) = group_pxn(ind_err_plus(end),i_test); % all error+ tests
            angvel_err_mins_plus(ns,2) = group_pxn(ind_err_plus(end),14); % all error+ tests
        end
        
    end
end

%% Plot
% Plot: sample1, 1st amd last correct test, err- and err+ tests
pxn_all = [pxn_samp1_crt(:,[1,2,4]),pxn_err_mins_plus];
angvel_all = [angvel_samp1_crt(:,[1,2,4]),angvel_err_mins_plus];

pxn_all_plot = pxn_all;
angvel_all_plot = angvel_all;

ind = find(isnan(mean(pxn_all_plot,2)));
pxn_all_plot(ind,:) = [];
angvel_all_plot(ind,:) = [];

data0 = pxn_all_plot;
% data0 = angvel_all_plot;

sample_color = [255,165,0]./255;
test_color = [111,57,214]./255;
test_err_color = [120,120,120]./255;
ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
data0_mean = mean(data0);
data0_sem = std(data0)./sqrt(size(data0,1));
hold on
bar(1,data0_mean(1),'FaceColor',sample_color,'EdgeColor', sample_color);
bar(2,data0_mean(2),'FaceColor',test_color,'EdgeColor', test_color);
bar(3,data0_mean(3),'FaceColor',test_color,'EdgeColor', test_color);
bar(4,data0_mean(4),'FaceColor',test_err_color,'EdgeColor', test_err_color);
bar(5,data0_mean(5),'FaceColor',test_err_color,'EdgeColor', test_err_color);
h1 = errorbar(1,data0_mean(1),data0_sem(1),'color','k','LineWidth',2);
h2 = errorbar(2,data0_mean(2),data0_sem(2),'color','k','LineWidth',2);
h3 = errorbar(3,data0_mean(3),data0_sem(3),'color','k','LineWidth',2);
h4 = errorbar(4,data0_mean(4),data0_sem(4),'color','k','LineWidth',2);
h5 = errorbar(5,data0_mean(5),data0_sem(5),'color','k','LineWidth',2);
errorbar_tick(h1, 0);
errorbar_tick(h2, 0);
errorbar_tick(h3, 0);
errorbar_tick(h4, 0);
errorbar_tick(h5, 0);
% plot([1,2,3],data0,'ko');
hold off
set(gca,'XTick',[1,2,3,4,5]);
xlim([0.5,5.5])
ylim([0.2,0.5])
set(gca,'XTickLabel',{'Sample1','1st crt test','Last crt test','Err- tests','Err+ tests'});
set(gca,'fontsize',20);

ylabel([{'Summation of P(x|n)'},{'from current location to reward location'}]);
% ylabel('Angular speed (rad/s)');

% title('Pre-reward')
% title('Pre-reward,Vel >5cm/s')
% title('Pre-reward,Vel >5cm/s and approaching stop location')

% Plot: err- and err+ tests
pxn_err_mins_plus_plot = pxn_err_mins_plus;
angvel_err_mins_plus_plot = angvel_err_mins_plus;

ind = find(~isnan(pxn_err_mins_plus_plot(:,1)));
pxn_err_mins_mean = mean(pxn_err_mins_plus_plot(ind,1));
pxn_err_mins_sem = std(pxn_err_mins_plus_plot(ind,1))./sqrt(length(ind));
angvel_err_mins_mean = mean(angvel_err_mins_plus_plot(ind,1));
angvel_err_mins_sem = std(angvel_err_mins_plus_plot(ind,1))./sqrt(length(ind));
ind = find(~isnan(pxn_err_mins_plus_plot(:,2)));
pxn_err_plus_mean = mean(pxn_err_mins_plus_plot(ind,2));
pxn_err_plus_sem = std(pxn_err_mins_plus_plot(ind,2))./sqrt(length(ind));
angvel_err_plus_mean = mean(angvel_err_mins_plus_plot(ind,2));
angvel_err_plus_sem = std(angvel_err_mins_plus_plot(ind,2))./sqrt(length(ind));

data0_mean = [pxn_err_mins_mean,pxn_err_plus_mean];
data0_sem = [pxn_err_mins_sem,pxn_err_plus_sem];
% data0_mean = [angvel_err_mins_mean,angvel_err_plus_mean];
% data0_sem = [angvel_err_mins_sem,angvel_err_plus_sem];

sample_color = [255,165,0]./255;
test_color = [111,57,214]./255;
test_err_color = [120,120,120]./255;
ffa = figure('Units','normalized','Position',[0 0 0.4 0.6]);
hold on
bar(1,data0_mean(1),'FaceColor',test_err_color,'EdgeColor', test_err_color);
bar(2,data0_mean(2),'FaceColor',test_err_color,'EdgeColor', test_err_color);
h1 = errorbar(1,data0_mean(1),data0_sem(1),'color','k','LineWidth',2);
h2 = errorbar(2,data0_mean(2),data0_sem(2),'color','k','LineWidth',2);
errorbar_tick(h1, 0);
errorbar_tick(h2, 0);
% plot([1,2,3],data0,'ko');
hold off
set(gca,'XTick',[1,2]);
xlim([0.5,2.5])
ylim([0.2,0.6])
set(gca,'XTickLabel',{'Err- tests','Err+ tests'});
set(gca,'fontsize',20);

ylabel([{'Summation of P(x|n)'},{'from current location to reward location'}]);
% ylabel('Angular speed (rad/s)');

%% Repeated-measured ANOVA
% For line graph
within_fact = fullfact([8,2]);
within_tbl2 = array2table(within_fact,'VariableNames',{'TrialN','TrialType'});
within_tbl2.TrialType = categorical(within_tbl2.TrialType);
within_tbl2.TrialN = categorical(within_tbl2.TrialN);

%rm1 = fitrm(StatsData_LineGraph,'T1_samp-T8_test ~ DirType','WithinDesign',within_tbl2);
%StatsOut.LineGraph.ANOVA = ranova(rm1,'WithinModel','TrialN*TrialType');
%StatsOut.LineGraph.posthoc = multcompare(rm1,'TrialN','by','DirType');

% Separate ahead and behind
in = strcmp(StatsData_LineGraph.DirType,'ahead');
rm1 = fitrm(StatsData_LineGraph(in,:),'T1_samp-T8_test ~ 1','WithinDesign',within_tbl2);
StatsOut.Ahead.ANOVA = ranova(rm1,'WithinModel','TrialN*TrialType');

%in = strcmp(StatsData_LineGraph.DirType,'behind');
%rm1 = fitrm(StatsData_LineGraph(in,:),'T1_samp-T8_test ~ 1','WithinDesign',within_tbl2);
%StatsOut.Behind.ANOVA = ranova(rm1,'WithinModel','TrialN*TrialType');

% For bar graph
% TrialType = table((1:3)','VariableNames',{'TrialType'});
% rm2 = fitrm(StatsData_BarGraph,'Sample1-Test_crtN ~ DirType','WithinDesign',TrialType);
% StatsOut.BarGraph.ANOVA = ranova(rm2);
% StatsOut.BarGraph.posthoc = multcompare(rm2,'TrialType','by','DirType');

%save([fig_out,'\Stats_PxnSum.mat'],'StatsOut','StatsData_LineGraph','within_tbl2')