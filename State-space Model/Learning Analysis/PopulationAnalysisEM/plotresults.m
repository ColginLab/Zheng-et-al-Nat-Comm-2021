%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program plots the output from runanalysis.m
% Anne Smith, March 16th, 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        p      (vectors)         mode of prob correct estimate from forward filter
%        p05,p95   (vectors)      conf limits of prob correct estimate from forward filter
%        b      (vectors)         mode of prob correct estimate from backward filter
%        b05,b95   (vectors)      conf limits of prob correct estimate from backward filter
%
%
%       pdistn          calculates the confidence limits (b05, b95, bmid, 
%                       bmode, pmatrix) of the EM estimate of the learning
%                       state process 

load('resultspopulation');
%raster of raw data
subplot(311)

data    = Responses;
lx      = size(data,1) + 1;
datanew = data;
datanew = [datanew; zeros(1,size(data,2))];
datanew = [datanew zeros(lx,1)];

pcolor(1-datanew)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Group effects

 len   = size(Bstarhat,2);
 b_est = Bnew(end);
 cnt   = 1;
 t = 0:len-1;
 figno = 1;

 figure(1); subplot(3,1,2); box on;
  
 bs = [b_est*ones(1, len)];
 
%plot non-overlapping trials blocks (length = clump)


%compute pdf of the population
 [p05a, p95a, pmida, pmodea, pmatrix] ...
      = pdistn(Bstarnew(1,:).*bs, (bs.^2).*squeeze(Wnew(1,1,:))',  muone, .5);

%plot 90% confidence limits in gray
 xx1   = [ t fliplr(t)];
 yy1   = [p05a(1:end) fliplr(p95a(1:end))];
 h     = fill(xx1, yy1, [0.7 0.8 0.8]);

 set(h, 'EdgeColor',[ 0.6 0.7 0.7]);

%plot mode in red
 hold on; plot(t, pmodea(1:end), 'r.-')

%find learning trial for group
  cback = find(p05a < BackgroundProb);
  if(~isempty(cback))
   if(cback(end) < size(Responses,2) )
    cback = cback(end);
   else
    cback = NaN;
   end
  else
   cback = NaN;
  end


fprintf(2,'IO(0.95) The learning trial for population is %d \n', cback)

%plot IO certainty

 subplot(3,1,3)
 hold on; plot(t,1 - pmatrix(1:end),'r.-')
 line([ 0 t(end)],[0.95 0.95]);
 axis([0 t(end)  0 1]);
 grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot individuals 

for j = 1:J

 figure(1); hold on;

 ratno = j;

for time  = 1:len

 rat   = ratno + 1;

 mean1 = [ Bstarnew(1, time) Bstarnew(rat, time)];
 cov1  = [ Wnew(1, 1, time) Wnew(1, rat, time);  ...
          Wnew(rat, 1, time) Wnew(rat, rat, time) ];

%compute random numbers from beta^j*x_time (they covary)

 numbermcs = 100000;
 r         = mvnrnd(mean1, cov1, numbermcs);
 r1        = r(:,1); 
 r2        = r(:,2);

 vv        = r1.*r2;  %products

 pp        = exp(muone + vv)./(1 +  exp(muone + vv)) ;

%   [n, bin ] = hist(pp, 10000);
%   [y, i ] = max(n);

   sortedpp = sort(pp);
   b05(time)  = sortedpp(.05*numbermcs); 
   b95(time)  = sortedpp(.95*numbermcs); 

   bmed(time)  = median(pp);

   %bmode(time) = bin(i);

   cback = find(b05 < BackgroundProb);

  if(~isempty(cback))
   if(cback(end) < size(Responses,2) )
    cback = cback(end);
   else
    cback = NaN;
   end
  else
   cback = NaN;
  end
  ltp05(ratno) = cback;

end


  figure(2); hold on;
  subplot(J,1,j); plot(t, b05,'k');
  %title(['Animal ' num2str(ratno) '   IO(0.95) learning trial ' num2str(ltp05(ratno))]);
  text(40, .2, ['LT ' num2str(ltp05(ratno))]);
  hold on; plot(t, b95,'k');
  hold on; plot(t, bmed,'r');
  hold on; line([0 80],[0.5 0.5]);
  axis([0 len-1 0 1.05]);

  mm = Responses(ratno,:);
  hold on; [y, x] = find(mm > 0);
  h = plot(x,y+0.03,'.'); set(h, 'MarkerFaceColor','b');
  set(h, 'MarkerEdgeColor', 'b');

  hold on; [y, x] = find(mm == 0);
  h = plot(x,1.03,'.'); set(h, 'MarkerFaceColor', [1.0 0.0 0.0]);
  set(h, 'MarkerEdgeColor', [1.0 .00 .00]);
  axis([0 len-1 0 1.05]);



  figure(1);
  subplot(312); hold on; box on;
  plot(t, bmed(1:end),'g-')
  axis([0 len 0 1]);

end


 subplot(312)
 hold on; line([0 80],[0.5 0.5]);
 
%plotmode of group over individuals
 hold on; plot(t, pmodea(1:end), 'r.-')

 clump  = 8; 
 plotno = 2;
 [x, y, v] = trialblocks(BackgroundProb, clump, 3, 1, plotno);
 
 axis([0 len-1 0 1]);

 subplot(313)
 hold on; line([0 80],[0.5 0.5]);
 axis([0 len-1 0 1]);

 fprintf(2, 'Individual learning trials \n')
 fprintf(2, '%d \n', ltp05')
 
 figure(1); subplot(311); ylabel('Subject')
 figure(1); subplot(312); ylabel('Pr(Correct Response)')
 figure(1); subplot(313); ylabel('IO Certainty'); xlabel('Trial Number')
 