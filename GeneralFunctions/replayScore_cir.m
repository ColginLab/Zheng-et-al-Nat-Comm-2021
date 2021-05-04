function score = replayScore_cir(pxn,xfitBin)

binrange = 5;
mid = floor(size(pxn,1)/2);
shift = xfitBin-mid;
pxn_shifted = pxn;
for ii = 1:length(xfitBin)
    pxn_shifted(:,ii) = circshift(pxn(:,ii),-shift(ii));
end

score = sum(sum(pxn_shifted(mid-binrange:mid+binrange,:)))/sum(sum(pxn_shifted));