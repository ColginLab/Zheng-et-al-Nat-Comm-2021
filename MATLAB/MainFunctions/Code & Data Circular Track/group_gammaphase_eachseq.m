% Edited on 06/28/2019
% load sequences that have been detected before, in each trial type
% for each sequence, get the sg/fg phase and cycle number for each spikes
function group_gammaphase_eachseq(DataDir, file_analysis_name_ext, file_analysis_out_ext, FiguresDir)
%parentfd = fileparts(mfilename('fullpath'));
%data_folder =[parentfd,'\GroupData\'];
%FiguresDir = [parentfd,'\GroupData Figures'];
%file_input1 = 'group_gammaTFR_eachseq_20190529.mat'; % 4 rats % remove jumping-out points, ind_approach = 3;
% file_input1 = 'group_gammaTFR_eachseq_lowfirEEG_20190519.mat'; % 4 rats % remove jumping-out points, ind_approach = 3; use EEG with low firing
% file_input1 = 'group_gammaTFR_eachseq_ds_stable_20190525.mat'; % 4 rats % remove jumping-out points, ind_approach = 3;
% file_input1 = 'group_gammaTFR_eachseq_ds_unstable_20190525.mat'; % 4 rats % remove jumping-out points, ind_approach = 3;
%file_input1 = strcat(data_folder,file_input1);

%file_analysis_out = 'data_gammaphase_20190628.mat';

load(file_analysis_name_ext)

numbins = 90; % Number of bins
bin_ang = 2*pi/numbins;
ang_bins = bin_ang/2:bin_ang:(2*pi-bin_ang/2);

directories_allData_v1

pos_prdc = nan(Nseq,2);
pos_real = nan(Nseq,2);
dist_prdc_real = nan(Nseq,1);
vel = nan(Nseq,1);
for nseq = 1:Nseq
    min0 = min(Pxn_all{nseq,4});
    max0 = max(Pxn_all{nseq,4});
    pos_prdc(nseq,:) = [min0,max0];
    
    min0 = min(Pxn_all{nseq,6});
    max0 = max(Pxn_all{nseq,6});
    pos_real(nseq,:) = [min0,max0];
    
    dist_prdc_real(nseq) = min(abs(angle(exp((Pxn_all{nseq,4}'-Pxn_all{nseq,6})*1i))));
    vel(nseq,1) = mean(Pxn_all{nseq,7});
end
ang_bins0 = mean(diff(Ang_RewardLoc_ontrack));

%%
ncell = 3;
nspk = 5;
nmove = 1;

ind_seq = {};
ind_seq{1} = find(data_info(:,3) == 2 & ... % 1(Pre-run); 2(Sample); 3(Test); 4(Post-test)
    data_info(:,4) <= 8 & ... % trial number
    pos_real(:,2) < data_info(:,8) &... % stop location
    data_info(:,10) >= ncell & data_info(:,11) >= nspk & ...  % limit active cells number and spikes
    data_info(:,6) == 1 &...  % only correct or error trials
    data_info(:,12) > 0 &... % slope    
    data_info(:,13) >= .6 &... % replay score
    para_all(:,1) <= 20 &... % max jump (in bins)
    para_all(:,4) >= nmove &... % total distance moved (in bins)
    dist_prdc_real <= 5*bin_ang); % minimal distance between predicted and real locations
title12{1} = 'Correct sample trials';

ind_seq{2} = find(data_info(:,3) == 2 & ... % 1(Pre-run); 2(Sample); 3(Test); 4(Post-test)
    data_info(:,4) <= 8 & ... % trial number
    pos_real(:,2) < data_info(:,8) &... % stop location
    data_info(:,10) >= ncell & data_info(:,11) >= nspk & ...  % limit active cells number and spikes
    data_info(:,6) == 0 &...  % only correct or error trials
    data_info(:,12) > 0 &... % slope    
    data_info(:,13) >= .6 &... % replay score    
    para_all(:,1) <= 20 &... % max jump (in bins)
    para_all(:,4) >= nmove &... % total distance moved (in bins)
    dist_prdc_real <= 5*bin_ang);
title12{2} = 'Error sample trials';

ind_seq{3} = find(data_info(:,3) == 3 & ... % 1(Pre-run); 2(Sample); 3(Test); 4(Post-test)
    data_info(:,4) <= 8 & ... % trial number
    pos_real(:,2) < data_info(:,8) &... % stop location
    data_info(:,10) >= ncell & data_info(:,11) >= nspk & ...  % limit active cells number and spikes
    data_info(:,6) == 1 &...  % only correct or error trials
    data_info(:,12) > 0 &... % slope    
    data_info(:,13) >= .6 &... % replay score    
    para_all(:,1) <= 20 &... % max jump (in bins)
    para_all(:,4) >= nmove &... % total distance moved (in bins)
    dist_prdc_real <= 5*bin_ang);
title12{3} = 'Correct test trials';

ind_seq{4} = find(data_info(:,3) == 3 & ... % 1(Pre-run); 2(Sample); 3(Test); 4(Post-test)
    data_info(:,4) <= 8 & ... % trial number
    pos_real(:,2) < data_info(:,8) &... % stop location
    data_info(:,10) >= ncell & data_info(:,11) >= nspk & ...  % limit active cells number and spikes
    data_info(:,6) == 0 &...  % only correct or error trials
    data_info(:,12) > 0 &... % slope    
    data_info(:,13) >= .6 &... % replay score    
    para_all(:,1) <= 20 &... % max jump (in bins)
    para_all(:,4) >= nmove &... % total distance moved (in bins)
    dist_prdc_real <= 5*bin_ang);
title12{4} = 'Error test trials';

% Downsample Method1: 
% downsample ind_seq so that sequence# in each category are same
% s = RandStream('mt19937ar','Seed',1); 
% min0 = Nseq;
% for itrial = 1:length(ind_seq)
%     min0 = min(min0,length(ind_seq{itrial}));
% end
% for itrial = 1:length(ind_seq)
%     if length(ind_seq{itrial}) > min0
%         ind = datasample(s,1:length(ind_seq{itrial}),min0,'Replace',false);
%         ind = sort(ind);
%         ind_new = ind_seq{itrial}(ind);
%         ind_seq_new{itrial} = ind_new;
%     elseif length(ind_seq{itrial}) == min0
%         ind_seq_new{itrial} = ind_seq{itrial};
%     end
% end


% % Downsample Method2: 
% % downsample ind_seq in each location area
% % so that sequence# within each location bin are same across trial types
% realloc_mean = nan(Nseq,1);
% predloc_max = nan(Nseq,1);
% for ii = 1:Nseq
%     realloc_mean(ii) = mean(Pxn_all{ii,6});
%     predloc_max(ii) = max(Pxn_all{ii,4});
% end
% dist_real_stop = realloc_mean-data_info(:,8);
% dist_pred_stop = predloc_max-data_info(:,8);
% 
% binw = mean(diff(Ang_RewardLoc_ontrack))*2;
% dist_bins = -binw*6:binw:binw*6;
% 
% ind_seq_new = cell(size(ind_seq));
% for i_r = 1:length(dist_bins)-1
%     ind_all = {};
%     len_all = [];
%     for itrial = 1:length(ind_seq)
%         dist_real_stop0 = dist_real_stop(ind_seq{itrial},:);
%         dist_pred_stop0 = dist_pred_stop(ind_seq{itrial},:);
%         ind_all{itrial} = find(dist_real_stop0 >= dist_bins(i_r) &...
%             dist_real_stop0 <= dist_bins(i_r+1) &...
%             dist_pred_stop0 >= dist_real_stop0);
%         len_all(itrial) = length(ind_all{itrial});
%     end
%     min0 = min(len_all);
%     
%     for itrial = 1:length(ind_seq)
%         if len_all(itrial) > min0
%             ind = datasample(1:len_all(itrial),min0,'Replace',false);
%             ind = sort(ind);
%             ind_new = ind_all{itrial}(ind);
%             ind_seq_new{itrial} = [ind_seq_new{itrial};ind_seq{itrial}(ind_new)];
%         elseif len_all(itrial) == min0
%             ind_seq_new{itrial} = [ind_seq_new{itrial};ind_seq{itrial}(ind_all{itrial})];
%         end
%     end
% end


% ind_seq = ind_seq_new;


%% group data
file_input1 = 'Data_gammaPhs.mat';
file_input2 = 'Data_gammaBandpass.mat';
file_input3 = 'Data_gammaPhs_cyclenum.mat';
file_input_cell = 'Cells_ds_ALLLaps_v2_vel_0.mat';
% file_input_cell = 'Cells_ALLLaps_v2_vel_0.mat';
TTList0 = 'TTList_dCA1_pyr_ds.txt';
dt = .04;
step = .01;

ind_seq0 = cat(1,ind_seq{:});
[ind_seq_2,I] = sort(ind_seq0);

Pxn0 = Pxn_all(ind_seq_2,:);
data_info0 = data_info(ind_seq_2,:);
ns = 0;
Data0 = [];
Data0_seqnum = [];
power_thr = NaN; % in z score; Set to NaN for no threshold
for nseq = 1:length(ind_seq_2)
    ns_new = data_info0(nseq,2);
    if ns_new ~= ns
        % a new session now
        ns = ns_new;
        path_ns = pathRats{ns};
        cd(path_ns);
        disp(['Currently in ' pathRats{ns}])
        
        load (file_input2)
        load (file_input3)
        load (file_input_cell)
        
        fid=fopen(TTList0);
        if (fid == -1)
            warning([ 'Could not open tfile ' TTList]);
        else
            % read the file names from the t-file list
            TT0 = ReadFileList(TTList0);
            numCells0 = length(TT0);
            if numCells0==1 && max(TT0{1}==-1)
                % no cells in ttlist
                numcells=0;
            else
                numcells=numCells0;
            end
        end
        fclose(fid);
        if size(spikes,1) ~= numcells
            warning(['cell number not matched']);
        end
        
        % use the channel with highest theta power for estimate of gamma
        % phase globally
        tetNew = CSClist_CA1_sort{ns}(1);
        
    end
    
    t_startstop = Pxn0{nseq,1}([1,end]);
    % Obtain peak slow/fast gamma power within a sequence event
    InRange = BP{1,1} >= t_startstop(1) & BP{1,1} <= t_startstop(2)+dt-step;
    SGpower_max = max(BP{2,4}(InRange));
    FGpower_max = max(BP{2,3}(InRange));
    
    Ncell = size(spikes,1);
    for nc = 1:Ncell
        ind = find(spikes{nc} >= t_startstop(1) & spikes{nc} <= t_startstop(2)+dt-step);
        if ~isempty(ind)
%             if TT0{nc,1}(4) ~= '_';
%                 tetNew = str2double(TT0{nc,1}(3:4));
%             elseif TT0{nc,1}(4) == '_';
%                 tetNew = str2double(TT0{nc,1}(3));
%             else
%                 disp('problem with something');
%             end
            nspikes = length(ind);
            for nspk0 = 1:nspikes
                [~,i] = min(abs(spikes{nc}(ind(nspk0))-Phase{1,1}));
                phsth_i = Phase{1,2}(i);
                phsf_i = Phase{1,3}(tetNew,i); phsfcyc_i = Phase{2,3}(tetNew,i);
                phss_i = Phase{1,4}(tetNew,i); phsscyc_i = Phase{2,4}(tetNew,i);
                temp = NaN(5,1);
                temp(1) = phsth_i;
                if FGpower_max >= power_thr
                    temp(2) = phsf_i;
                    temp(3) = phsfcyc_i;
                end
                if SGpower_max >= power_thr
                    temp(4) = phss_i;
                    temp(5) = phsscyc_i;
                end
                if isnan(power_thr)
                    temp(2) = phsf_i;
                    temp(3) = phsfcyc_i;
                    temp(4) = phss_i;
                    temp(5) = phsscyc_i;
                end
                Data0 = [Data0,temp];
                Data0_seqnum = [Data0_seqnum,ind_seq_2(nseq)];
            end
        end
    end
end
% Data:
% row1: theta phase
% row2: fast gamma phase
% row3: fast gamma cycle
% row4: slow gamma phase
% row5: slow gamma cycle

%% Re-assign data to different trial conditions (e.g. sample correct)
Ntrial = length(ind_seq);
for itrial = 1:Ntrial
    a = ind_seq{itrial};
    I = ismember(Data0_seqnum,a);
    Data_seqnum{itrial} = Data0_seqnum(I);
    Data{itrial} = Data0(:,I);
end

save(file_analysis_out_ext, 'Data','Data_seqnum','title12','-v7.3');

%% Ploting
xbin = 10;
xaxis0=-360+xbin/2:xbin:360-xbin/2;
ybin = 10;
yaxis0=0+ybin/2:ybin:720-ybin/2;
axis0={xaxis0,yaxis0};

cycLimit = [-1 1];
Ntrial = length(Data);
h1 = figure('Units','normalized','Position',[0 0 1 1]);
for itrial = 1:Ntrial
    in = ~isnan(Data{itrial}(4,:));
    ph_th=Data{itrial}(1,in);
    ph_sg=Data{itrial}(4,in);
    cyc_sg=Data{itrial}(5,in);
    ind_cycle = unique(cyc_sg);
    ind_cycle = ind_cycle(~isnan(ind_cycle) & ind_cycle>=cycLimit(1) & ind_cycle<=cycLimit(2));
    z_sg=[];
    for ifig=1:length(ind_cycle)
        in = cyc_sg==ind_cycle(ifig);
        z=(hist3([[ph_th(in)',ph_sg(in)'];[ph_th(in)',ph_sg(in)'+360]],axis0))';
        z=z./sum(sum(z));
        z_sg(:,:,ifig) = z;
    end

    xaxis0=-360:xbin:360;
    yaxis0=0:xbin:720;
    axis0={xaxis0,yaxis0};
    for ifig=1:length(ind_cycle)
        subplot(Ntrial,length(ind_cycle),ifig+(itrial-1)*length(ind_cycle))
        z=z_sg(:,:,ifig);
        imagesc(xaxis0,yaxis0,smooth2a(normz(z,0),5,5));
        set(gca,'YDir','normal','FontSize',12);

        xlim([-90,270])
        ylim([0,360])
        set(gca,'Xtick',[-90 270])
        if ifig == 1
            set(gca,'Ytick',0:180:360);
            ylabel(sprintf([title12{itrial},'\nSlow gamma phase']));
        else
            set(gca,'Ytick',[])
        end
        if itrial == Ntrial
            xlabel('Theta phase');
        end
        if itrial == 1
            title(strcat({'Cycle '}, num2str(ind_cycle(ifig))));
        end
        axis square
    end
end
saveas(h1,[FiguresDir,'\SupFigure5a'],'fig')
saveas(h1,[FiguresDir,'\SupFigure5a'],'epsc')
close(h1)

cycLimit = [-2 2];
Ntrial = length(Data);
h2 = figure('Units','normalized','Position',[0 0 1 1]);
for itrial = 1:Ntrial
    in = ~isnan(Data{itrial}(2,:));
    ph_th=Data{itrial}(1,in);
    ph_fg=Data{itrial}(2,in);
    cyc_fg=Data{itrial}(3,in);
    ind_cycle = unique(cyc_fg);
    ind_cycle = ind_cycle(~isnan(ind_cycle) & ind_cycle>=cycLimit(1) & ind_cycle<=cycLimit(2));
    z_fg=[];
    for ifig=1:length(ind_cycle)
        in = cyc_fg==ind_cycle(ifig);
        z=(hist3([[ph_th(in)',ph_fg(in)'];[ph_th(in)',ph_fg(in)'+360]],axis0))';
        z=z./sum(sum(z));
        z_fg(:,:,ifig) = z;
    end

    xaxis0=-360:xbin:360;
    yaxis0=0:xbin:720;
    axis0={xaxis0,yaxis0};
    for ifig=1:length(ind_cycle)
        subplot(Ntrial,length(ind_cycle),ifig+(itrial-1)*length(ind_cycle))
        z=z_fg(:,:,ifig);
        imagesc(xaxis0,yaxis0,smooth2a(normz(z,0),5,5));
        set(gca,'YDir','normal','FontSize',12);

        xlim([-90,270])
        ylim([0,360])
        set(gca,'Xtick',[-90 270])
        if ifig == 1
            set(gca,'Ytick',0:180:360);
            ylabel(sprintf([title12{itrial},'\nFast gamma phase']));
        else
            set(gca,'Ytick',[])
        end
        if itrial == Ntrial
            xlabel('Theta phase');
        end
        if itrial == 1
            title(strcat({'Cycle '}, num2str(ind_cycle(ifig))));
        end
        axis square
    end
end
saveas(h2,[FiguresDir,'\SupFigure5b'],'fig')
saveas(h2,[FiguresDir,'\SupFigure5b'],'epsc')
close(h2)
end