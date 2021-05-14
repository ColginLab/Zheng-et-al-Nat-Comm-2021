function BayesianDecodingExample(DataDir, FiguresDir)
%parentfd = fileparts(mfilename('fullpath'));
%Outdir = [parentfd,'\GroupData Figures'];
%Outdir = 'E:\ColginLab\MATLAB\Figures\GroupData';
%Outdir = 'E:\ColginLab\Figures\Figures3_4';
% use the overall place field excluding pre-ruuning trials as a decoder
file_input = 'BayesData_CircMap_v1_dt40ms_10ms_1cells_1spk_v2_ds.mat';
file_input_speed = 'Data_angle_ontrack.mat';  % input the speed information

directories_allData_v1

% Selected_session = [16, 25, 44, 24];
% Selected_Trial_all = {[7 3],[2 3],[6 3],[3 4];... % Test trial number
%                       [4 3],[5 3],[1 2],[3 4]};   % Sample trial number
% tWin_all = {[-4 0], [-2 0], [-3 0],[-2 0];...
%             [-4 0], [-2 0], [-4 0],[-2 0]}; % time window with respect to stop; in sec

Selected_session = 25;
Selected_Trial_all = {[2 3];... % Test trial number
                      [5 3]};   % Sample trial number
tWin_all = {[-2 0];...
            [-2 0]}; % time window with respect to stop; in sec       

FigLabels = {'\Fig3a', '\Fig4a'};
for ii = 1:length(Selected_session)
    path_ns = pathRats{Selected_session(ii)};
    str = strsplit(path_ns,'\');
    RatID = str{end-2};
    DateID = str{end-1};
    cd(path_ns);
    load (file_input);
    load(file_input_speed,'data_angle');

    trackdata_ns = trackdata{Selected_session(ii)};
    load(trackdata_ns);

    % Only plot sample-test laps for now
    scores_sample = scores{1,2};
    scores_test = scores{1,3};
    nlaps = size(scores_sample,1);
    %% Plot results
    scores_all = {scores_test,scores_sample};
    ts_stop_all = {ts_test_reward, ts_sample_reward};
    TrialID = {'test','sample'};
    for tt = 1:length(TrialID) % loop through sample and test
        Selected_Trial = Selected_Trial_all{tt,ii};
        scores0 = scores_all{tt};
        ts_stop0 = ts_stop_all{tt};
        tWin = tWin_all{tt,ii}; % time window with respect to stop; in sec
        h1 = figure('Units','normalized','Position',[0 0 .6 .5]);
        for ss = 1:length(Selected_Trial)
            nl = Selected_Trial(ss);
            pxn = scores0{nl,3};
            pos = scores0{nl,4};
            posAxis = scores0{nl,5};
            pxn_t = scores0{nl,6};

            InRange = find(pxn_t>=ts_stop0(nl)/1e6+tWin(1) & pxn_t<=ts_stop0(nl)/1e6)+tWin(2);
            subplot(length(Selected_Trial),1,ss)
            hold on
            imagesc(pxn_t(InRange)-pxn_t(InRange(end)),posAxis,pxn(:,InRange))
            axis xy tight
            plot(pxn_t(InRange)-pxn_t(InRange(end)),posAxis(pos(InRange)),'b-','LineWidth',1)
            plot(xlim,[mode(ang_sample_reward_ontrack) mode(ang_sample_reward_ontrack)],'g--','LineWidth',1)
            set(gca,'CLim',[0 .3])
            colormap hot
            cb = colorbar;
            if sign_correct_test(nl)
                title(['Correct ',TrialID{tt},' trial'])
            else
                title(['Error ',TrialID{tt},' trial'])
            end
            if ss == 1
                set(gca,'XTick',[])
            elseif ss == length(Selected_Trial)
                xlabel('Time before stop (s)')
                ylabel('Position (rad)')
                ylabel(cb,'Posterior probability')
            end
        end
        saveas(h1,[FiguresDir, FigLabels{tt}],'fig')
        saveas(h1,[FiguresDir, FigLabels{tt}],'epsc')
        close(h1)
    end
end
end