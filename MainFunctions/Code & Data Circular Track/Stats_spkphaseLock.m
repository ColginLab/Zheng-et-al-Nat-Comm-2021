parentfd = fileparts(mfilename('fullpath'));
outfd = [parentfd,'\GroupData Figures'];
inputfile = [parentfd,'\GroupData\data_gammaphase_20190628.mat'];
load(inputfile)

%% Extract data
TrialType  = [];
ThetaPhase = [];
FGPhase    = [];
FGCycleN   = [];
SGPhase    = [];
SGCycleN   = [];
for tt = 1:length(Data)
    TrialType = cat(1,TrialType,ones(size(Data{tt},2),1)*tt);
    ThetaPhase = cat(1,ThetaPhase,Data{tt}(1,:)');
    FGPhase = cat(1,FGPhase,Data{tt}(2,:)');
    FGCycleN = cat(1,FGCycleN,Data{tt}(3,:)');
    SGPhase = cat(1,SGPhase,Data{tt}(4,:)');
    SGCycleN = cat(1,SGCycleN,Data{tt}(5,:)');
end

StatsData = [array2table(TrialType), array2table(ThetaPhase), ...
             array2table(FGPhase), array2table(FGCycleN),...
             array2table(SGPhase), array2table(SGCycleN)];

h1 = figure('Units','normalized','Position',[0 0 .4 1]);
s = RandStream('mt19937ar','Seed',1); 
nShuffle = 5000;
%% Test Slow gamma phase coding across successive cycles         
in = StatsData.SGCycleN > -2 & StatsData.SGCycleN < 2;
alpha = StatsData.SGPhase(in)/180*pi;
x = StatsData.SGCycleN(in)+2;

% Circular-linear correlation between cycle number and slow gamma phase
TrialType = StatsData.TrialType(in);
TrialType_uq = unique(TrialType);
rho_all = NaN(1,length(TrialType_uq));
for ii = 1:length(TrialType_uq)
    in2 = TrialType == TrialType_uq(ii);
    rho_all(ii) = circ_corrcl(alpha(in2), x(in2));
end

% Obtain null distribution
rho_all_shu = NaN(nShuffle,length(TrialType_uq));
for ss = 1:nShuffle
    ind_shu = randperm(s,length(TrialType),length(TrialType));
    TrialType_shu = TrialType(ind_shu);
    for ii = 1:length(TrialType_uq)
        in2 = TrialType_shu == TrialType_uq(ii);
        rho_all_shu(ss,ii) = circ_corrcl(alpha(in2), x(in2));
    end
end
reset(s)
CI_null = [prctile(rho_all_shu(:),2.5), prctile(rho_all_shu(:),97.5)];

% plot results
subplot(2,1,1); hold on
plot([0 length(TrialType_uq)+1],[CI_null(1) CI_null(1)],'k--')
plot([0 length(TrialType_uq)+1],[CI_null(2) CI_null(2)],'k--')
for ii = 1:length(TrialType_uq)
    plot(ii,rho_all(ii),'ro')
end
axis square
ylabel('Cirular-linear regression')
set(gca,'XTick',[],'YLim',[0 .09])

%% Test fast gamma phase strength across trial types
in = StatsData.FGCycleN > -3 & StatsData.FGCycleN < 3;
alpha = StatsData.FGPhase(in)/180*pi;
x = StatsData.FGCycleN(in)+3;

% Mean vector length
TrialType = StatsData.TrialType(in);
TrialType_uq = unique(TrialType);
r_all = NaN(1,length(TrialType_uq));
for ii = 1:length(TrialType_uq)
    in2 = TrialType == TrialType_uq(ii);
    r_all(ii) = circ_r(alpha(in2));
end

% Obtain null distribution
r_all_shu = NaN(nShuffle,length(TrialType_uq));
for ss = 1:nShuffle
    ind_shu = randperm(s,length(TrialType),length(TrialType));
    TrialType_shu = TrialType(ind_shu);
    for ii = 1:length(TrialType_uq)
        in2 = TrialType_shu == TrialType_uq(ii);
        r_all_shu(ss,ii) = circ_r(alpha(in2));
    end
end
reset(s)
CI_null = [prctile(r_all_shu(:),2.5), prctile(r_all_shu(:),97.5)];

%% plot results
subplot(2,1,2); hold on
plot([0 length(TrialType_uq)+1],[CI_null(1) CI_null(1)],'k--')
plot([0 length(TrialType_uq)+1],[CI_null(2) CI_null(2)],'k--')
for ii = 1:length(TrialType_uq)
    plot(ii,r_all(ii),'ro')
end
axis square
ylabel('Mean vector length')
set(gca,'XTick',1:length(TrialType_uq),'XTickLabel',title12,'YLim',[0 .05])

%% save output
saveas(h1,[outfd,'\Stats_spkphaseLock'],'fig')
saveas(h1,[outfd,'\Stats_spkphaseLock'],'epsc')
close(h1)
save([outfd,'\Stats_spkphaseLock.mat'],'StatsData')