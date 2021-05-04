function RunSpeedDist(Runfd,Outdir)
Condition = {'Prerunning';'Sample';'Test';'Postrunning'};
outdir = strcat(Outdir,'\RunSpeedDist');
if ~isdir(outdir)
    mkdir(outdir)
end

wSpeedDist_all = [];
RatID_Label = [];
Condition_Label = [];
for rr = 1:length(Runfd)
    ind = strfind(Runfd{rr},'Rat');
    RatID = Runfd{rr}(ind:ind+5);     
    posIndfile = strcat(Runfd{rr},'\Data_angle_ontrack.mat');
    load(posIndfile)
    speedaxis = linspace(.01,80);
    for tt = 1:length(Condition)
        for laps = 1:length(data_angle{tt})
            speed = data_angle{tt}{laps}(:,3);
            wSpeedDist = kernelDensity(speed,speedaxis,2);
            wSpeedDist_all = cat(1,wSpeedDist_all,wSpeedDist);
            RatID_Label = cat(1,RatID_Label,{RatID});
            Condition_Label = cat(1,Condition_Label,Condition(tt));
        end
    end    
end

%% Plot results
RatID = unique(RatID_Label);
condition = {'Sample';'Test'};
colorLabel = {'b';'r'};
h = figure;
set(h,'OuterPosition',[255.4,376.2,1030.4,386.4])
for rr = 1:length(RatID)
    ha = zeros(length(condition),1);
    for cc = 1:length(condition)
        in = strcmp(RatID_Label,RatID{rr}) & strcmp(Condition_Label,condition{cc});
        sample = bootstrp(5000,@mean,wSpeedDist_all(in,:));
        [CI_u, CI_l] = CorrectionForMultipleComparsion(sample);
        upperLimit = CI_u-mean(wSpeedDist_all(in,:));
        lowerLimit = mean(wSpeedDist_all(in,:))-CI_l;
        
        subplot(1,length(RatID),rr); hold on
        temp = shadedErrorBar(speedaxis,mean(wSpeedDist_all(in,:)),[upperLimit;lowerLimit],strcat(colorLabel{cc},'-'));
        ha(cc) = temp.mainLine;
        axis tight square
    end
    xlabel('Speed(cm/s)')
    title(RatID{rr})
    if rr == length(RatID)
        legend(ha,condition)
    elseif rr == 1
        ylabel('Proportion')
    end
end
saveas(h,strcat(outdir,'\SpeedDist_allrat'),'epsc')
saveas(h,strcat(outdir,'\SpeedDist_allrat'),'png')
close(h)

h = figure; hold on
set(h,'OuterPosition',[255.4,376.2,419.2,386.4])
ha = zeros(length(condition),1);
for cc = 1:length(condition)
    in = strcmp(Condition_Label,condition{cc});
    sample = bootstrp(5000,@mean,wSpeedDist_all(in,:));
    [CI_u, CI_l] = CorrectionForMultipleComparsion(sample);
    upperLimit = CI_u-mean(wSpeedDist_all(in,:));
    lowerLimit = mean(wSpeedDist_all(in,:))-CI_l;

    temp = shadedErrorBar(speedaxis,mean(wSpeedDist_all(in,:)),[upperLimit;lowerLimit],strcat(colorLabel{cc},'-'));
    ha(cc) = temp.mainLine;
    axis tight square
end
set(gca,'Box','off')
xlabel('Speed(cm/s)')
ylabel('Proportion')
legend(ha,condition)
saveas(h,strcat(outdir,'\SpeedDist_combined'),'epsc')
saveas(h,strcat(outdir,'\SpeedDist_combined'),'png')
close(h)

%% ANOVA
in = ismember(Condition_Label,condition);
MeanSpeed = mean(wSpeedDist_all(in,:),2);
ratID = RatID_Label(in);
conditions = Condition_Label(in);
t=[cell2table(ratID) cell2table(conditions) array2table(MeanSpeed)];
lm = fitlm(t,'MeanSpeed ~ ratID*conditions');
Stats.ANOVA = anova(lm);

%% Export analysis results and all the necessary scripts
FunctionPath = mfilename('fullpath');
FunctionName = mfilename;
OutScript_dir = strcat(outdir,'\Scripts');
if isdir(OutScript_dir)
    rmdir(OutScript_dir,'s')
end
mkdir(OutScript_dir)
[fList] = matlab.codetools.requiredFilesAndProducts(FunctionPath);
for ff = 1:length(fList)
    [~,fname,ext] = fileparts(fList{ff});
    copyfile(fList{ff},strcat(OutScript_dir,'\',fname,ext));
end
save(strcat(outdir,'\wSpeedDist.mat'),'wSpeedDist_all','RatID_Label','Condition_Label','FunctionName')
save(strcat(outdir,'\Stats.mat'),'Stats')

function wDist = kernelDensity(x,xaxis,sigma)

delta = repmat(x,1,length(xaxis))-repmat(xaxis,length(x),1);
W = exp(-0.5*delta.*delta/sigma^2);
W = W./repmat(sum(W,2),1,size(W,2)); %normalize each kernel so the ones on the edge have same contribution
wDist = sum(W)./sum(sum(W));
wDist = wDist/sum(wDist);
        
function [CI_u, CI_l] = CorrectionForMultipleComparsion(sample)
%Simultaneous bounds are wider than separate bounds, because it is more stringent to require that the entire 
%curve be within the bounds than to require that the curve at a single predictor value be within the bounds.
    [sample_rank] = tiedrank(sample);
    sample_rank=sample_rank/size(sample_rank,1);

    sample_min=prctile(min(sample_rank,[],2),2.5);
    sample_max=prctile(max(sample_rank,[],2),97.5);

    CI_u=zeros(1,size(sample,2));
    for ii = 1:size(sample,2);
        CI_u(ii) = min(sample(sample_rank(:,ii)>=sample_max,ii));
    end

    CI_l=zeros(1,size(sample,2));
    for ii = 1:size(sample,2);
        CI_l(ii) = max(sample(sample_rank(:,ii)<=sample_min,ii));
    end