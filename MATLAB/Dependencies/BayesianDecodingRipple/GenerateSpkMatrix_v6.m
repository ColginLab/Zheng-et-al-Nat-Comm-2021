function [SpikeMatrix,spkFR,posaxis,peakFR_pos,posbin,peakFR_sortInd,CellID] = ...
    GenerateSpkMatrix_v6(eeg2posIndfile,posIndfile,varargin)
% input
%
% eeg2posIndfile: eeg index of position
% posIndfile: the position index file named "Data_angle_ontrack"
% varargin(1): bin size of the position; Default is 5
% varargin(2): whether or not to sort; Default is to sort
% varargin(3): spike index that matches eeg index
% varargin(4): Trials to use (1 is Prerunning; 2 is sample; 3 is test; 4 is TrialsPost)
%-----------------------------------------------
% output
%
% SpikeMatrix: the matrix that contains the time of each spike from
%       each cell; in cell ID X sampling time point for each cell content
% spkFR: position tuning curve of each cell; in cell ID X position bin
% posaxis: position value of each position bin
% peakFR_pos: position value of peak firing rate position from each cell
% posbin_ind: current position bin of each sampling time point
% peakFR_sortInd: sorted index of cell ID


%Default
binsize = 5; % in linear track is cm; in circular track is degree
run_thr = 5; % cm/sec
sortornot = true;
trialsToUse = [1 2 3 4];
kernel = gausskernel(5,2);
if nargin >2
    binsize = varargin{1};
end

if nargin >3
    sortornot = varargin{2};
end

spikeIndfile = varargin{3};

if nargin >5
    trialsToUse = varargin{4};
end

load(spikeIndfile);
load(eeg2posIndfile);
load(posIndfile);               % where does posIndfile come from?

% find the data within trials
if nargin >5
    include = false(size(eeg2pos));
    for ii = 1:length(trialsToUse)
        trial_data = data_angle{trialsToUse(ii)};
        lapN = size(trial_data,1);
        for ll = 1:lapN
            lap_ts = trial_data{ll,1}(:,1);
            % position index
            ind_range = [find((data_angle_all(:,1)-lap_ts(1))== 0),...
                         find((data_angle_all(:,1)-lap_ts(end))== 0)];
            % translate into eeg index
            eegind_range = [find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1),...
                            find(eeg2pos>=min(ind_range) & eeg2pos<=max(ind_range),1,'last')];
            include(min(eegind_range):max(eegind_range)) = true;
        end
    end
else
    include = true(size(eeg2pos));
end 

% obtain position occupancy
ang_pos = data_angle_all(eeg2pos(include),2)/2/pi*360;
speed = data_angle_all(eeg2pos(include),3);
posaxis = 0:binsize:360;
run_ind = speed>=run_thr;
[~,posbin] = histc(ang_pos,posaxis);
occupancy = hist(posbin(run_ind),1:length(posaxis)-1);
occupancy = conv_cir(occupancy,kernel)';
occupancy = occupancy+.0001; % Add offset to prevent zeros

% obtain total number of spike
CellID = fieldnames(spike_ind);
nSpike = 0;
for cc = 1:length(CellID)
    % remove any repeating spike indices. This happens when there are missing eeg
    % samples, so multiple spike timestamp map onto the same eeg index
    spkISI = diff(spike_ind.(CellID{cc}));
    if any(spkISI==0)
        spike_ind.(CellID{cc})(spkISI==0)=[];
        disp(['Remove repeating spike index from ' spikeIndfile])
    end
    nSpike = nSpike+length(spike_ind.(CellID{cc}));
end

% create sparse matrix
cellID_ind = zeros(nSpike,1);
sampleTS_ind = zeros(nSpike,1);
nSpike = 0;
for cc = 1:length(CellID)
    nSpike_range = nSpike+1:nSpike+length(spike_ind.(CellID{cc}));
    cellID_ind(nSpike_range) = cc;
    sampleTS_ind(nSpike_range) = spike_ind.(CellID{cc});    

    nSpike = nSpike+length(spike_ind.(CellID{cc}));
end
SpikeMatrix = sparse(cellID_ind,sampleTS_ind,true(size(cellID_ind)),length(CellID),length(eeg2pos));
SpikeMatrix = SpikeMatrix(:,include);

Spkposbin_count = zeros(length(CellID),length(posaxis)-1);
for cc = 1:length(CellID)
    Spkposbin_count(cc,:) = hist(posbin(SpikeMatrix(cc,:) & run_ind'),1:length(posaxis)-1);
end
spkFR = Spkposbin_count./repmat(occupancy,size(Spkposbin_count,1),1);
[~, peakFR_pos] = max(spkFR,[],2);
spkFR = spkFR*2000; %convert to Hz
peakFR_pos = peakFR_pos*binsize;

peakFR_sortInd=[];
if sortornot
    [~, peakFR_sortInd] = sort(peakFR_pos);
    peakFR_pos = peakFR_pos(peakFR_sortInd);    
    spkFR = spkFR(peakFR_sortInd,:);
    SpikeMatrix = SpikeMatrix(peakFR_sortInd,:);
    CellID = CellID(peakFR_sortInd);
end
posaxis(end)=[]; % remove last position bin, which overlaps with the first bin on the circular track    

%% plot raster
    % linear_ind = find(SpikeMatrix_sort);
    % [cellID,sampletime] = ind2sub(size(SpikeMatrix_sort),linear_ind);
    % figure;
    % hold on;
    % cmap = lines;
    % for ii=1:size(sampletime,1)
    %     line([sampletime(ii) sampletime(ii)],[cellID(ii)-.5 cellID(ii)+.5],'color',cmap(cellID(ii),:));
    % end