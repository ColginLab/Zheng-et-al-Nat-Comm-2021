function [p05, p95, pmid, pmode, pmatrix] = pdistn(x, s, mu, background_prob);
%pdist is a helper function that calculates the confidence limits of the EM
%estimate of the learning state process.  For each trial, the function
%constructs the probability density for a correct response.  It then builds
%the cumulative density function from this and computes the p values of
%the confidence limits
%
%variables:
%   xx(ov)   EM estimate of learning state process
%   ss(ov)   EM estimate of learning state process variance
%   pmatrix  vector of the level of certainty the ideal observer has that performance is better than chance at each trial  
%   dels     bin size of the probability density p values
%   pr       bins of the probability density distribution
%   fp       p{k|j}, probability density of the probability of a correct response at trial k     (equation B.3)*
%   pdf      probability density function      
%   sumpdf   cumulative density function of the pdf
%   lowlimit index of the p value that gives the lower 95% confidence
%            bound
%   highlimit   index of the p value that gives the upper 95% confidence
%               bound
%   middlimit   index of the p value that gives the 
%   p05      the p value that gives the lower 95% confidence bound
%   p95      the p value that gives the upper 95% confidence bound
%   pmid     the p value that gives the 50% confidence bound
%   pmode    the p value that gives the highest probability density

pmatrix = [];
for ov = 1:size(x,2)

 xx = x(ov);
 ss = s(ov);
 
 dels=1e-4;

 pr  = dels:dels:1-dels;
 term1 = 1./(sqrt(2*pi*ss) * (pr.*(1-pr)));
 term2 = exp(-1/(2*ss) * (log (pr./((1-pr)*exp(mu))) - xx).^2);
 pdf = term1 .* term2;
 pdf = dels * pdf;
 
 
% Integrate the pdf
 sumpdf = cumtrapz(pdf);
% sumpdf = cumsum(pdf);

lowlimit  = find(sumpdf>0.05);
if(~isempty(lowlimit) )
lowlimit  = lowlimit(1);
else
lowlimit  = 1;
end

highlimit = find(sumpdf>0.95);
% highlimit = find(sumpdf>0.995);
if(~isempty(highlimit) )
if(length(highlimit)>1)
highlimit = highlimit(1)-1;
else
highlimit =  highlimit(1);
end
else
highlimit = length(pr);
end

middlimit = find(sumpdf>0.5);
if(~isempty(middlimit))
middlimit = middlimit(1);
else
middlimit = length(pr);
end


 p05(ov)   = pr(lowlimit(1));
 p95(ov)   = pr(highlimit(1));
 pmid(ov)  = pr(middlimit(1));
 [y,i]     = max(pdf);
 pmode(ov) = pr(i);
 

 pmatrix =[pmatrix; sumpdf];
 
end

inte = fix(background_prob/dels);

pmatrix = pmatrix(:, inte);