function [p_x_n] = BayesianDecoder(spkraster,pos_tun,timewindow,step,sampFreq)
%This version of Baysian Decoder assumes a uniform prior probability over
%position. The position estimate depends only on spikes observed in the 
%given time window, and not on any previous output of the estimator.
% Reference: Davidson, TJ et. al., Neuron 2009 and Zhang, et al., Journal
% of Neurophysiology 1998
%-------------------------------------------------------------------------%
% Input
%   spkraster: the spike raster intended for decoding. Should be in number 
%   of cells X time bin
%   pos_tun: position tuning curve for each cell. Should be in number of
%   cells X posbin
%   timewindow: the length of the time for decoding. Should be in second  
%   step: time step for time window advancement. Should be in sec
%   sampFreq: the sampling frequency of spike raster. Should be in Hz
% Output
%   p_x_n: the probability distribution of position given the observed
%   spikes. position bin X timebin

%-------------------------------------------------------------------------%
if size(spkraster,1)~= size(pos_tun,1)
    error('the number of cells should be the same')
end

%make nsameple a odd number
if mod(sampFreq,2) == 0
    nsample = timewindow*sampFreq+1;
else
    nsample = timewindow*sampFreq;
end

%main part of the decoding
p_x_n = NaN(size(pos_tun,2),ceil(size(spkraster,2)/step/sampFreq));
ind2pxn = 0;
for pp = floor(nsample/2)+1:step*sampFreq:size(spkraster,2)-floor(nsample/2)
    ind2pxn = ind2pxn+1;
    
    range = pp-floor(nsample/2):pp+floor(nsample/2)-1;
    if sum(sum(spkraster(:,range),2)~=0) >= 1
        n = repmat(sum(spkraster(:,range),2),[1 size(pos_tun,2)]);
        p_x_n(:,ind2pxn) = (prod(pos_tun.^n,1))'.*exp(-timewindow*sum(pos_tun,1))'; %the poisson equation to infer posterior probability
    end
    
    p_x_n(:,ind2pxn) = p_x_n(:,ind2pxn)/nansum(p_x_n(:,ind2pxn));
        
end
% end condition
if size(spkraster,2) > max(range)
    if sum(sum(spkraster(:,max(range)+1:end))) > 0 
        ind2pxn = ind2pxn+1;
        n = repmat(sum(spkraster(:,max(range)+1:end),2),[1 size(pos_tun,2)]);
        p_x_n(:,ind2pxn) = (prod(pos_tun.^n,1))'.*exp(-timewindow*sum(pos_tun,1))';        
    end
    
    p_x_n(:,ind2pxn) = p_x_n(:,ind2pxn)/nansum(p_x_n(:,ind2pxn));
    if ind2pxn<size(p_x_n,2)
        p_x_n(:,ind2pxn+1:end) = repmat(p_x_n(:,ind2pxn),1,size(p_x_n,2)-ind2pxn);
    end
    
end

