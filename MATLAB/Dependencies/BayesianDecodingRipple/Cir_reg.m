function [r2,calphase,xaxis,p,slope] = Cir_reg(Pxn,xbins,tbins,bins2use)
% CZ: shuffle the position bins instead of time bins, which was uesd in
% Davison et al Neuron 2009 paper

%changed to use actual time and space dimensions rather than
%just bin numbers
%also added in com measeure and its difference from the predicted line as
%an additional measure of accuracy of prediction


%how many surrogate r2 values to get:
drawnum = 1000;

all = [];
s = RandStream('mt19937ar','Seed',1); 
for i = 1:size(bins2use,2) %passing through time bins
   bin = bins2use(i);
   pdtemp = Pxn(:,bin)';
   if isempty(find(isnan(pdtemp), 1))
       randlocs = xbins(  randsample(s,length(pdtemp),drawnum,true,pdtemp)  );
       randlocs(:,2) = tbins(bin);
       all = [all;randlocs];
   end
end

if ~isnan(all) & ~isempty(all)
    %% calculate the original sequence
    x = all(:,2);
    y = all(:,1);
    [para,~,p] = circ_lin_regress(x, y, 8);
    xaxis = sort(unique(x));
    calphase = 2*pi*para(1,1)*xaxis+para(1,2);
    slope = para(1,1)*2*pi;
%     calphase = calphase*360/2/pi;
    over = calphase>=2*pi;
    under = calphase<0;
    calphase(over) = calphase(over)-2*pi;
    calphase(under) = calphase(under)+2*pi;
    
    % compute the residual
    yresidual = zeros(drawnum,size(xaxis,1));
    ytotal = zeros(drawnum,size(xaxis,1));
    for ii = 1:size(xaxis,1)
        yfit = repmat(calphase(ii),drawnum,1);
         yresidual(:,ii) = circdistance(y(x==xaxis(ii)), yfit, 1);
         ytotal(:,ii) = circdistance(y(x==xaxis(ii)), circ_mean(y), 1);
    end
    SSresid = sum(sum(yresidual));
    SStotal = sum(sum(ytotal));
    r2 = 1 - SSresid/SStotal;
    if r2<0
        r2 = NaN;
    end
        
end    
end
