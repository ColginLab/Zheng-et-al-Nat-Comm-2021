function PlaceFieldsDist_alignedToRestBox(Outdir)

outdir = strcat(Outdir,'\PlaceFieldsDist_alignedToRestBox');
if ~isdir(outdir)
    mkdir(outdir)
end

inputdir = strcat(Outdir,'\BayesianDecoding_ripple');
if ~isdir(inputdir)
    error('Please run "RunBayesianDecoding_ripple" first')
end
inputdir_preRun = strcat(Outdir,'\BayesianDecoding_ripple_PreRunTuning');
if ~isdir(inputdir_preRun)
    error('Please run "RunBayesianDecoding_ripple_PreRunTuning" first')
end

temp = load(strcat(inputdir,'\BayesianDecodingResult.mat'));
[PosTuning,CellID] = ExtractPosTuning(temp);
temp = load(strcat(inputdir_preRun,'\BayesianDecodingResult.mat'));
[PosTuning_preRun,CellID_preRun,posaxis] = ExtractPosTuning(temp);

[CellID,ia,ib] = intersect(CellID,CellID_preRun);
PosTuning = PosTuning(ia,:);
PosTuning_preRun = PosTuning_preRun(ib,:);

% Sort cell by positions of maximal firing rate
[~,maxInd] = max(PosTuning,[],2);
[~,sortInd] = sort(maxInd);
PosTuning_after = PosTuning(sortInd,:);
PosTuning_preRun_after  = PosTuning_preRun(sortInd,:);

h = Makeplots(PosTuning_after,PosTuning_preRun_after,posaxis);
saveas(h,strcat(outdir,'\PlaceFieldDis_after'),'epsc')
saveas(h,strcat(outdir,'\PlaceFieldDis_after'),'png')
close(h)

% Sort cell by positions of maximal firing rate
[~,maxInd] = max(PosTuning_preRun,[],2);
[~,sortInd] = sort(maxInd);
PosTuning_before = PosTuning(sortInd,:);
PosTuning_preRun_beofre  = PosTuning_preRun(sortInd,:);

h2 = Makeplots(PosTuning_before,PosTuning_preRun_beofre,posaxis);
saveas(h2,strcat(outdir,'\PlaceFieldDis_before'),'epsc')
saveas(h2,strcat(outdir,'\PlaceFieldDis_before'),'png')
close(h2)

% Make correlation plot
[h3, rho] = MakeCorrelationPlot(PosTuning,PosTuning_preRun,posaxis);
saveas(h3,strcat(outdir,'\PlaceFieldDisCorr'),'epsc')
saveas(h3,strcat(outdir,'\PlaceFieldDisCorr'),'png')
close(h3)

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
save(strcat(outdir,'\PlotOutput.mat'),'PosTuning_preRun','PosTuning','posaxis','rho','CellID')

%% Helper function
function [PosTuning,CellID,posaxis] = ExtractPosTuning(input)
PosTuning = [];
CellID = [];
ratID = fieldnames(input.BayesianDecodingResult);
posbinsize = input.BayesianDecodingResult.Rat139.D20170220CT2.param.posbinsize/2/pi*360;
posaxis = 0:posbinsize:360-posbinsize;
RestBoxPos = 0;
for rr = 1:length(ratID)
    dateID = fieldnames(input.BayesianDecodingResult.(ratID{rr}));
    for dd = 1:length(dateID)
        posTuning = input.BayesianDecodingResult.(ratID{rr}).(dateID{dd}).PosTuning;
        cellID = input.BayesianDecodingResult.(ratID{rr}).(dateID{dd}).cellID;
        shift = match(RestBoxPos,posaxis)-length(posaxis)/2;
        PosTuning = cat(1,PosTuning,circshift(posTuning,[0,-shift]));
        CellID = cat(1,CellID,strcat(ratID{rr},dateID{dd},cellID));
    end
end
posaxis = posaxis-180;

function h = Makeplots(PosTuning,PosTuning_preRun,posaxis)
% Normalize position tuning
PosTuning_norm = PosTuning./repmat(max(PosTuning,[],2),1,size(PosTuning,2));
PosTuning_preRun_norm = PosTuning_preRun./repmat(max(PosTuning_preRun,[],2),1,size(PosTuning_preRun,2));

% Make plots
h = figure;
set(h,'OuterPosition',[598.6,-8.6,792.8,873.6])
h1 = subplot(4,2,[1 3 5]);
imagesc(posaxis,1:length(PosTuning_preRun_norm),PosTuning_preRun_norm);  hold on
colormap hot
plot([0,0],[1 length(PosTuning_preRun_norm)],'LineStyle','--','LineWidth',1,'Color','w')
title('Before reward presentation')
ylabel('CA1 cell indices')
axis tight xy

h2 = subplot(4,2,[2 4 6]);
imagesc(posaxis,1:length(PosTuning_norm),PosTuning_norm); hold on
plot([0,0],[1 length(PosTuning_preRun_norm)],'LineStyle','--','LineWidth',1,'Color','w')
title('After reward presentation')
axis tight xy
pos = get(h2,'Position');
hb = colorbar;
ylabel(hb,'Normalized firing rate')
set(h2,'Position',pos)

h3 = subplot(4,2,7);
ha = plot(posaxis,median(PosTuning_preRun_norm));
axis tight
set(ha,'LineWidth',1);
xlabel('Aligned angular position (deg)')
ylabel('Median normalized FR')
linkaxes([h1 h3],'x')

h4 = subplot(4,2,8);
ha = plot(posaxis,median(PosTuning_norm));
axis tight
set(ha,'LineWidth',1);
xlabel('Aligned angular position (deg)')
linkaxes([h2 h4],'x')

function [h, rho] = MakeCorrelationPlot(PosTuning,PosTuning_preRun,posaxis)
% Normalize position tuning
PosTuning_norm = PosTuning./repmat(max(PosTuning,[],2),1,size(PosTuning,2));
PosTuning_preRun_norm = PosTuning_preRun./repmat(max(PosTuning_preRun,[],2),1,size(PosTuning_preRun,2));

rho = corr(PosTuning_preRun_norm,PosTuning_norm);

h = figure; hold on
imagesc(posaxis,posaxis,rho)
axis square tight xy
plot([0 0],ylim,'w--')
plot(xlim,[0 0],'w--')
set(gca,'CLim',[0 1])
ylabel('Before reward position')
xlabel('After reward position')
hb = colorbar;
ylabel(hb,'Pearson r')
