% Edit on 07/23/2019
% Modified from doall_get_gammaTFR_eachtrial.m and Stats_Pxn_gamma_v2.m
% Select continuous swipe ahead or behind sequences (Pheiffer and Foster 2013 Nature)
% Aim to check out gamma change with swipe ahead distance, and the
% prediction of reward location

%clear

parentfd = fileparts(mfilename('fullpath'));
% use the overall place field excluding pre-running trials as a decoder
% Downsample the place cells so that fields equally distrubuted on the track
file_input_speed = 'Data_angle_ontrack.mat';  % input the speed information
%data_folder = [parentfd,'\GroupData\'];
data_folder = 'E:\ColginLab\Data Analysis\GroupData\';

% file_input = 'BayesData_CircMap_v1_dt40ms_10ms_1cells_1spk_v2_ds.mat';
% file_input_TFR = 'Data_TFR_gamma_lowfir.mat';
% file_spike_input1 = 'Cells_ds_ALLLaps_v2_vel_0.mat';  % used to get all spikes
% file_output1 = 'group_gammaTFR_eachseq_lowfirEEG_20190519.mat'; % 4 rats % remove jumping-out points, ind_approach = 3; use EEG with low firing

file_input = 'BayesData_CircMap_v1_dt40ms_10ms_1cells_1spk_v2_ds.mat';
file_input_TFR = 'Data_TFR_gamma_SR.mat';
file_spike_input1 = 'Cells_ds_ALLLaps_v2_vel_0.mat';  % used to get all spikes
file_output1 = 'group_gammaTFR_eachseq_20190723.mat'; % 4 rats % remove jumping-out points, ind_approach = 3;

% file_input = 'BayesData_CircMap_v1_dt40ms_10ms_1cells_1spk_v2_ds_stable.mat';
% file_input_TFR = 'Data_TFR_gamma.mat';
% file_spike_input1 = 'Cells_ds_stable_ALLLaps_v2_vel_0.mat';  % used to get all spikes
% file_output1 = 'group_gammaTFR_eachseq_ds_stable_20190525.mat'; % 4 rats % remove jumping-out points, ind_approach = 3; use EEG with low firing

% file_input = 'BayesData_CircMap_v1_dt40ms_10ms_1cells_1spk_v2_ds_unstable.mat';
% file_input_TFR = 'Data_TFR_gamma.mat';
% file_spike_input1 = 'Cells_ds_unstable_ALLLaps_v2_vel_0.mat';  % used to get all spikes
% file_output1 = 'group_gammaTFR_eachseq_ds_unstable_20190525.mat'; % 4 rats % remove jumping-out points, ind_approach = 3; use EEG with low firing

file_output1 = strcat(data_folder,file_output1);

directories_allData_v1
% fig_folder = 'Figure_Pxn_gammaTFR_sequence_20190518\'; % 4 rats, zscore gamma over all sessions, remove jumping-out points
fig_folder = 'Figure_Pxn_gammaTFR_sequence_20190529\'; % 4 rats, zscore gamma over all sessions, remove jumping-out points

dt = .04;
step = .01;

% Default parameters
maxJump_thr = 30;  % in position bins
timeWin = 6;     % in time bin
timeStep = 1;     % in time bin
Distance_thr = 0; % in position bins
jump_prop_thr = 0;

Nseq = 0;
Fs = 2000;
Pxn_all = {};
TFRz_all = {};
data_info = [];
para_all = [];

sign_plot = 0;
ffa = figure('Units','normalized','Position',[0 0 0.35 1]);
for ns = 1:isession
    path_ns = pathRats{ns};
    cd(path_ns);
    disp(path_ns)
    
    csclist_ns = CSClist_CA1{ns};
    
    if exist(file_input,'file')>0 && ~isempty(csclist_ns)
        load(file_input);
        load(file_input_TFR);
        trackdata_ns = trackdata{ns};
        load(trackdata_ns,'Ind_rewardloc','Ind_rewardloc_sample','Ind_rewardloc_test','Ind_rewardloc_posttest','Diam_inner',...
            'Sign_correct_sample','Sign_correct_test','Sign_correct_posttest',...
            'Ang_RewardLoc_ontrack','n_prerunning');
        S1 = load(file_spike_input1);  % used to get all spikes
        spikes = S1.spikes(:,1);
        Ncell = size(spikes,1);
        
        ang_vel_limit = 5/(Diam_inner/2); % 5cm/s
        Ind_rewardloc_all = [{nan(n_prerunning,1)},Ind_rewardloc_sample,Ind_rewardloc_test,Ind_rewardloc_posttest];
        Sign_correct_all = [{nan(n_prerunning,1)},Sign_correct_sample,Sign_correct_test,Sign_correct_posttest];
        
        for nseq = 1:size(scores,2)
            score_nseg = scores{nseq};
            for nl = 1:size(score_nseg,1)
                t_nl = score_nseg{nl,2};
                pxn_nl = score_nseg{nl,3};
                
                % Select continuous swipe ahead sequences
                [onset, offset, para] = DetectSequenceEvents_cz(score_nseg{nl,3},...
                    maxJump_thr,timeWin,timeStep,Distance_thr,jump_prop_thr,true);
                Nseq0 = length(onset);
                
                % for each sequence
                for i = 1:Nseq0
                    ind1 = onset(i); ind2 = offset(i);
                    
                    % get Pxn
                    pxn_nseq = score_nseg{nl,3}(:,ind1:ind2);
                    % get TFR
                    ind11 = find(abs(TFR_gamma{2,1}-score_nseg{nl,6}(ind1)) < 1/Fs);
                    if ~isempty(ind11)
                        [~,ind11] = min(abs(TFR_gamma{2,1}-score_nseg{nl,6}(ind1)));
                    end
                    ind22 = find(abs(TFR_gamma{2,1}-score_nseg{nl,6}(ind2)-(dt-step)) < 1/Fs);
                    if ~isempty(ind22)
                        [~,ind22] = min(abs(TFR_gamma{2,1}-score_nseg{nl,6}(ind2)-(dt-step)));
                    end
                    
                    if isempty(ind11) || isempty(ind22) 
                        continue
                    else
                        Nseq = Nseq+1;
                        
                        % get number of active cells and spikes
                        nspk = nan(Ncell,1);
                        for nc = 1:Ncell
                            ind = find(spikes{nc,1} >= score_nseg{nl,6}(ind1) &...
                                spikes{nc,1} < score_nseg{nl,6}(ind2)+dt-step);
                            nspk(nc) = length(ind);
                        end
                        ind = find(nspk>0);
                        ncell_nseq = length(ind);
                        nspk_nseq = sum(nspk);
                        
                        % get sequence fit
                        bin2use = find(~isnan(sum(pxn_nseq)));
                        [~,calphase,~,~,slope] = ...
                            Cir_reg(pxn_nseq,score_nseg{nl,5},score_nseg{nl,6}(ind1:ind2)-score_nseg{nl,6}(ind1),bin2use);
                        [~,posbin_ind] = histc(calphase,score_nseg{nl,5});
                        replaySC = replayScore_cir(pxn_nseq,posbin_ind'+1);
                        
                        % get Pxn for this sequence
                        [~, decodedbin] = max(pxn_nseq);
                        Pxn_all(Nseq,:) = [{score_nseg{nl,6}(ind1:ind2)},...
                            {pxn_nseq},...
                            {decodedbin},...
                            {score_nseg{nl,5}(decodedbin)'},...
                            {score_nseg{nl,4}(ind1:ind2)},...
                            {score_nseg{nl,5}(score_nseg{nl,4}(ind1:ind2))},...
                            {score_nseg{nl,8}(1,ind1:ind2)},...
                            {score_nseg{nl,8}(4,ind1:ind2)},...
                            {posbin_ind'+1}];

                        % get TFR for each sequence
                        %tfr_nseq = TFR_gamma{2,2}(:,ind11:ind22);
                        TFR_gamma{1,2} = squeeze(TFR_gamma{1,2});
                        TFR_gamma{2,2} = squeeze(TFR_gamma{2,2});
                        tfr_nseq = TFR_gamma{2,2}(:,ind11:ind22);
                        TFRz_all(Nseq,:) = [{TFR_gamma{2,1}(1,ind11:ind22)},...
                            {tfr_nseq}];
                        
                        % get info for each sequence
                        
                        Sign_correct_nl = Sign_correct_all{nseq}(nl);
                        Ind_rewardloc_nl = Ind_rewardloc_all{nseq}(nl);
                        if ~isnan(Sign_correct_nl)
                            data_info(Nseq,:) = [Ind_Rat(ns),...
                                ns,...
                                nseq,...
                                nl,...
                                Ind_rewardloc,...
                                Sign_correct_nl,...
                                Ind_rewardloc_nl-Ind_rewardloc,...
                                Ang_RewardLoc_ontrack(Ind_rewardloc_nl),...
                                Ang_RewardLoc_ontrack(Ind_rewardloc),...
                                ncell_nseq,...
                                nspk_nseq,...
                                slope,...
                                replaySC];
                        else
                            data_info(Nseq,:) = [Ind_Rat(ns),...
                                ns,...
                                nseq,...
                                nl,...
                                Ind_rewardloc,...
                                Sign_correct_nl,...
                                nan,...
                                nan,...
                                Ang_RewardLoc_ontrack(Ind_rewardloc),...
                                ncell_nseq,...
                                nspk_nseq,...
                                slope,...
                                replaySC];
                        end
                        
                        % in sample lap, use correct sign in its paired test lap
                        if nseq == 2
                            Sign_correct_nl = Sign_correct_all{3}(nl);
                            data_info(Nseq,6) = Sign_correct_nl;
                        end
                        
                        % get para for each sequence
                        para_all(Nseq,:) = para(i,:);
                        
                        % plot each sequence
                        if sign_plot
                            % find reward location and stop location
                            ang_bins = score_nseg{nl,5};
                            [~,ind_stop] = min(abs(data_info(Nseq,8)-ang_bins));
                            [~,ind_reward] = min(abs(data_info(Nseq,9)-ang_bins));
                            
                            filename = strcat('Seq',num2str(Nseq),...
                                '-Seg',num2str(nseq),...
                                '-Err',num2str(data_info(Nseq,7)));
                            
                            subplot(7,1,1:3);
                            imagesc(pxn_nseq)
                            colorbar
                            axis xy
                            set(gca,'XTick',[1,size(pxn_nseq,2)]);
                            hold on
                            plot(score_nseg{nl,4}(ind1:ind2),'--','Color','w');
                            plot(ind_reward*ones(1,size(pxn_nseq,2)),'--','Color','r');
                            if Sign_correct_nl == 1 || isnan(Sign_correct_nl)
                                set(gca,'YTick',[ind_reward]);
                                set(gca,'YTickLabel','reward');
                            elseif ind_stop < ind_reward
                                set(gca,'YTick',[ind_stop,ind_reward]);
                                set(gca,'YTickLabel',{'stop','reward'});
                                plot(ind_stop*ones(1,size(pxn_nseq,2)),'--','Color','y');
                            elseif ind_stop > ind_reward
                                set(gca,'YTick',[ind_reward,ind_stop]);
                                set(gca,'YTickLabel',{'reward''stop'});
                                plot(ind_stop*ones(1,size(pxn_nseq,2)),'--','Color','y');
                            end
                            hold off
                            set(gca,'fontsize',14);
                            title(filename,'fontweight','normal');
                            
                            subplot(7,1,4:6);
                            imagesc(TFRz_all{Nseq,1},freq_TFR,tfr_nseq)
                            colorbar
                            axis xy
                            set(gca,'XTick',TFRz_all{Nseq,1}([1,end]));
                            set(gca,'YTick',[8,20:20:100]);
                            set(gca,'fontsize',14);
                            ylabel('Frequency (Hz)')
                            
                            subplot(7,1,7)
                            plot(1:size(pxn_nseq,2),Pxn_all{Nseq,7});
                            colorbar
%                             ylim([0,60])
                            xlim([1,size(pxn_nseq,2)])
                            set(gca,'XTick',[1,size(pxn_nseq,2)]);
                            set(gca,'fontsize',14);
                            xlabel('Time(s) or Time bins')
                            ylabel('Speed (cm/s)')
                            
                            bmpImage = strcat(parentfd,'\',fig_folder,filename,'.bmp');
                            f = getframe(gcf);
                            [pic, cmap] = frame2im(f);
                            imwrite(pic,bmpImage,'bmp');
                            clf
                        end
                    end
                end
            end
        end
        clear scores TFR_gamma score_nseg tfr_nseq S1
    end
end

save(file_output1,'-v7.3');