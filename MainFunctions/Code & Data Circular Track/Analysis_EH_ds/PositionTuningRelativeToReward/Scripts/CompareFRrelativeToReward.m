function CompareFRrelativeToReward(Path,trackdata,AnalysisFD,spktime_fd)

outdir = strcat(AnalysisFD,'\PositionTuningRelativeToReward');
if ~isdir(outdir)
    mkdir(outdir)
end
posbinsize = 4;
posaxis = 0:posbinsize:360-posbinsize;
FRkernel = gausskernel(5,2);
FR_preRun_session = zeros(length(Path),length(posaxis));
FR_AfterRewardRun_session = zeros(length(Path),length(posaxis));
FRdiff_session = zeros(length(Path),length(posaxis));
for rr = 1:length(Path)
    load(trackdata{rr});
    % Obtain reward position
    sign_correct_sample(isnan(sign_correct_sample)) = false;
    rewardPos = mode(ang_sample_reward_ontrack(logical(sign_correct_sample)));
    % Obtain position tunning
    eeg2posIndfile = strcat(Path{rr},'\',spktime_fd,'\eeg2pos.mat');
    posIndfile = strcat(Path{rr},'\Data_angle_ontrack.mat');
    spikeIndfile = strcat(Path{rr},'\',spktime_fd,'\spike_index.mat');
    [~,FR_preRun] = GenerateSpkMatrix_v6(eeg2posIndfile,posIndfile,posbinsize,false,spikeIndfile,1);    
    [~,FR_AfterRewardRun] = GenerateSpkMatrix_v6(eeg2posIndfile,posIndfile,posbinsize,false,spikeIndfile,[2 3 4]);
    % Smoothen position firing rate
    for ii = 1:size(FR_preRun,1)
        FR_preRun(ii,:) = conv_cir(FR_preRun(ii,:),FRkernel);
        FR_AfterRewardRun(ii,:) = conv_cir(FR_AfterRewardRun(ii,:),FRkernel);
    end
    % align to reward position
    shift = match(rewardPos/2/pi*360,posaxis)-length(posaxis)/2;
    FR_preRun = circshift(FR_preRun,[0,-shift]);
    FR_AfterRewardRun = circshift(FR_AfterRewardRun,[0,-shift]);
    
    rmv = max(FR_preRun,[],2) < 1 & max(FR_AfterRewardRun,[],2) < 1;
    FR_preRun_session(rr,:) = median(FR_preRun(~rmv,:));
    FR_AfterRewardRun_session(rr,:) = median(FR_AfterRewardRun(~rmv,:));
    FRdiff_session(rr,:) = median(FR_AfterRewardRun(~rmv,:)-FR_preRun(~rmv,:));
    
end
posaxis = posaxis - 180;

%% Export analysis results and all the necessary scripts
FunctionPath = mfilename('fullpath');
FunctionName = mfilename;
OutScript_dir = strcat(outdir,'\Scripts');
if isdir(OutScript_dir)
    rmdir(OutScript_dir,'s')
end
mkdir(OutScript_dir)
[fList] = matlab.codetools.requiredFilesAndProducts(FunctionPath);
for ff = 1:length(fList)
    [~,fname,ext] = fileparts(fList{ff});
    copyfile(fList{ff},strcat(OutScript_dir,'\',fname,ext));
end
save(strcat(outdir,'\FiringRelativeToReward.mat'),'FR_preRun_session','FR_AfterRewardRun_session','FRdiff_session','posaxis','FunctionName')
%% Make plots
h = figure;
set(h,'OuterPosition',[615,144,690,929])
subplot(2,1,1); hold on
CI_pre = bootci(5000,@median,FR_preRun_session);
CI_post = bootci(5000,@median,FR_AfterRewardRun_session);
ha = shadedErrorBar(posaxis,median(FR_preRun_session),[CI_pre(2,:)-median(FR_preRun_session); median(FR_preRun_session)-CI_pre(1,:)],'b');
ha2 = shadedErrorBar(posaxis,median(FR_AfterRewardRun_session)+.00001,[CI_post(2,:)-median(FR_AfterRewardRun_session); median(FR_AfterRewardRun_session)-CI_post(1,:)],'r');
legend([ha.mainLine ha2.mainLine],'Before reward','After reward')
ylabel('Median firing rate (Hz)')
axis tight square
subplot(2,1,2);
CI_FRdiff = bootci(5000,@median,FRdiff_session);
shadedErrorBar(posaxis,median(FRdiff_session),[CI_FRdiff(2,:)-median(FRdiff_session); median(FRdiff_session)-CI_FRdiff(1,:)]);
ylabel('After reward FR - Pre FR (Hz)')
xlabel('Angular position aligned to reward position (deg)')
axis tight square

saveas(h,strcat(outdir,'\FRcomparison'),'epsc');
saveas(h,strcat(outdir,'\FRcomparison'),'png');
close(h)

    