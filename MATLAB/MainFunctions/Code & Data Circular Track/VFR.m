function [VFR] = VFR(zTFR,V,vbins,freqVec,varargin)
sigma = 2;
if nargin>4
   sigma = varargin{1};
end
delta = repmat(V,1,length(vbins))-repmat(vbins,length(V),1);
W = exp(-0.5*delta.*delta/sigma^2);
D = ones(length(freqVec),length(V))*exp(-0.5*delta.*delta/sigma^2);
VFR = zTFR*W./D;