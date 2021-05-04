%-------------------------------------------------------------------------------------
%plot the figures
load('resultsindividual');
t=1:size(p,2)-1;

% figure(1);  clf;

%plot learning curve
subplot(2,1,1);
hold on;
plot(t, pmode(2:end),'r-');
plot(t, p05(2:end),'k', t, p95(2:end), 'k');
if(MaxResponse == 1)
    hold on; [y, x] = find(Responses > 0);
    h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor','k');
    set(h, 'MarkerEdgeColor', 'k');
    hold on; [y, x] = find(Responses == 0);
    h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor', [0.75 0.75 0.75]);
    set(h, 'MarkerEdgeColor', 'k');
    axis([1 t(end)  0 1.05]);
else
    hold on; plot(t, Responses./MaxResponse,'ko');
    axis([1 t(end)  0 1]);
end
line([1 t(end)], [BackgroundProb  BackgroundProb ]);
title(['IO(0.95) Learning trial = ' num2str(cback) '   Learning state process variance = ' num2str(SigE^2) ]);
xlabel('Trial Number')
ylabel('Probability of a Correct Response')

%plot IO certainty
subplot(2,1,2)
plot(t,1 - pmatrix(2:end),'k')
line([ 1 t(end)],[0.90 0.90]);
line([ 1 t(end)],[0.99 0.99]);
line([ 1 t(end)],[0.95 0.95]);
axis([1 t(end)  0 1]);
grid on;
xlabel('Trial Number')
ylabel('Certainty')

pdata = [t; p05(2:end); pmid(2:end); p95(2:end)]';
