function PlotCrossValidation(InputFile,Outdir)
% InputFile = 'I:\MATLAB\MainFunctions\Code & Data Circular Track\Analysis_EH_ds\ConfusionMatrix';
load(InputFile)
ConMat = [];
Error_sum = [];
errorAxis = 0:.05:pi;
RatID = fieldnames(ConfusionMat);
for rr = 1:length(RatID)
    DateID = fieldnames(ConfusionMat.(RatID{rr}));
    for dd = 1:length(DateID)
        mat = ConfusionMat.(RatID{rr}).(DateID{dd}).confusionMatrix;
        error = ConfusionMat.(RatID{rr}).(DateID{dd}).error;
        posAxis = ConfusionMat.(RatID{rr}).(DateID{dd}).posaxis;
        
        c = histc(error,errorAxis);
        error_sum = cumsum(c)/sum(c);
        
        ConMat = cat(3,ConMat,mat);
        Error_sum = cat(1,Error_sum,error_sum);
    end
end

%% Plot results
ffa1 = figure('Units','normalized','Position',[0 0 0.6 0.5]);
subplot(1,2,1)
imagesc(posAxis,posAxis,nanmean(ConMat,3))
axis square xy
colormap rvgray
xlabel('Real position (rad)')
ylabel('Decoded position (rad)')
cb = colorbar;
ylabel(cb,'Mean probability')
title(['n = ',num2str(size(Error_sum,1))])

subplot(1,2,2)
plot(errorAxis,Error_sum,'Color',[.8 .8 .8])
axis square
colorbar
xlim([0 pi])
xlabel('Decoded error (rad)')
ylabel('Proportion')
set(gca,'Box','off')

saveas(ffa1,[Outdir,'\CrossValidation'],'fig')
saveas(ffa1,[Outdir,'\CrossValidation'],'epsc')
close(ffa1)