function [x, y, v] = trialblocks(I, clump, n1, n2, plotno);
%trialblocks is a helper that plots the error bars of the learning curve
%plot

J   = size(I,1);
len = size(I, 2);
x =[]; y = []; v=[];

for j = 1: fix(len/clump)

  x = [x (j-1)*clump+clump/2];
  y1 = I(:, (j-1)*clump+1:j*clump);
  y2 = reshape(y1, 1, J*clump);
  y = [y mean(y2)];
  v = [v std(y2)/sqrt(clump*J)];

end

if(plotno == 0)
	return
end

figure(1); 
subplot(n1, n2, plotno); hold on;
hold on;
errorbar(x,y,v, 'k');
box on;
