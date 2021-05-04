function [onset, offset,para] = DetectSequenceEvents_cz(pxn,varargin)
%% Extract sequential activity from continuous Bayesian decoding
% Input
% pxn: posterior probability. p x t, where p is position bin and t is time bin
% Modified from DetectSequenceEvents.m
% do not accept swipe backward in any bin


% Default parameters
maxJump_thr = 10;  % in position bins
timeWin = 10;     % in time bin
timeStep = 1;     % in time bin
Distance_thr = 15; % in position bins
jump_prop_thr = 0.5;
BothDir = false; % include both direction?

if nargin > 1
    maxJump_thr = varargin{1};
end
if nargin > 2
    timeWin = varargin{2};
end
if nargin > 3
    timeStep = varargin{3};
end
if nargin > 4
    Distance_thr = varargin{4};
end
if nargin > 5
    jump_prop_thr = varargin{5};
end
if nargin > 6
    BothDir = varargin{6};
end

onset = [];
offset = [];
for ii = 1:timeStep:size(pxn,2)
    range = ii:ii+timeWin-1;
    if max(range)>size(pxn,2)
        break
    end
    emptyBin = isnan(pxn(1,range));
    if any(emptyBin)        
        continue
    end
    [~, decodedPos] = max(pxn(:,range));
    jump = diff(decodedPos);
    jump_prop = sum(jump>0)/length(jump);
    down_prop = sum(jump<0)/length(jump);
    maxJump = nanmax(jump);
    if BothDir
        Distance = abs(decodedPos(end)-decodedPos(1));% include both directions
        if maxJump <= maxJump_thr && Distance >= Distance_thr && jump_prop >= jump_prop_thr
           onset = cat(1,onset,ii);
           offset = cat(1,offset,range(end));
        end
    else
        Distance = decodedPos(end)-decodedPos(1);% include only positive directions        
        if maxJump <= maxJump_thr && Distance >= Distance_thr && jump_prop >= jump_prop_thr &&...
            down_prop == 0
           onset = cat(1,onset,ii);
           offset = cat(1,offset,range(end));
        end
    end
    

end
% Combine overlapping events
if length(onset) > 1
    isi = onset(2:end)-offset(1:end-1);
    merge = find(isi <= 0);
    onset(merge+1) = [];
    offset(merge) = [];
end

para = [];
for ii = 1:length(onset)
    timeWin = offset(ii)-onset(ii)+1;
    [~, decodedPos] = max(pxn(:,onset(ii):offset(ii)));
    jump = abs(diff(decodedPos));
    jump_prop = sum(jump>0)/length(jump);
    maxJump = nanmax(jump);
    if BothDir
        Distance = abs(decodedPos(end)-decodedPos(1));% include both directions
    else
        Distance = decodedPos(end)-decodedPos(1);% include only positive directions
    end
    para(ii,:) = [maxJump,timeWin,timeStep,Distance,jump_prop];
end
