parentfd = fileparts(mfilename('fullpath'));
ind = strfind(parentfd,'MATLAB');
%Dir = [parentfd(1:ind-1),'Circular Track\'];
%DataDir = 'E:\ColginLab\Data\Circular Track\';

isession = 0;
Ind_Rat = [];
pathRats = {};
trackdata = {};
CSClist_CA1 = {};
CSClist_CA1_sort = {};
CSClist_CA1_lowfir = {};
CSClist_R2 = {};
rewardlast = [];
rewardcurrent = [];

%% Rat 139
isession_Rat139 = 0;

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-02-20-CT-2\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170220_CT_tracking_2.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[1,6,8,9,12,14,16,17];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[6,17,12,16,1,14,9,8];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[6,12,17];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 12;
rewardcurrent_Rat139(isession_Rat139,1) = 10;
% 36 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [1,0,0,1,0,0,0,1]

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-02-21-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170221_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[1,6,8,9,14,16,17];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[6,17,16,14,1,9,8];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[6,17];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 10;
rewardcurrent_Rat139(isession_Rat139,1) = 14;
% 32 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [0,0,1,0,0,0,0,1]

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-02-22-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170222_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[6,8,9,10,12,14,16,17];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[6,17,12,16,14,10,9,8];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[6,10,12,14];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 9;
rewardcurrent_Rat139(isession_Rat139,1) = 11;
% 15 CA1 place cells, 0 CA1 interneurons
% sign_correct_test = [0,0,1,1,1,1,1,0]

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-02-25-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170225_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[];
% csclist_Rat139_CA1{isession_Rat139,1}=[6,8,9,10,11,14,16,17];
% csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[6,11];
csclist_Rat139_R2{isession_Rat139,1} = [];
rewardlast_Rat139(isession_Rat139,1) = 9;
rewardcurrent_Rat139(isession_Rat139,1) = 12;
% 21 CA1 place cells, 0 CA1 interneurons
% sign_correct_test = [1,1,1,1,1,1,1,1]
% The EEG are fluctuating and sometimes flowing out of range.

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-02-26-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170226_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[];
% csclist_Rat139_CA1{isession_Rat139,1}=[6,8,9,10,11,14];
% csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[6,9];
csclist_Rat139_R2{isession_Rat139,1} = [];
rewardlast_Rat139(isession_Rat139,1) = 12;
rewardcurrent_Rat139(isession_Rat139,1) = 8;
% 21 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [1,1,0,0,0,0,0,1]
% The EEG are fluctuating and sometimes flowing out of range.

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-02-27-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170227_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[];
% csclist_Rat139_CA1{isession_Rat139,1}=[6,7,8,10,11,14,16,17];
% csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[6,7,16,17];
csclist_Rat139_R2{isession_Rat139,1} = [];
rewardlast_Rat139(isession_Rat139,1) = 8;
rewardcurrent_Rat139(isession_Rat139,1) = 13;
% 24 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [0,1,1,0,0,1,1,1]
% The EEG are fluctuating and sometimes flowing out of range.

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-02-28-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170228_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[6,7,8,9,10,11,14,16,17];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[6,17,14,16,11,10,9,7,8];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[7,9,14,16,17];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 13;
rewardcurrent_Rat139(isession_Rat139,1) = 10;
% 35 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [NaN,1,1,1,1,0,0,NaN]
% The computer tracking in post-trials was not good. He actually stopped at the correct location, but he didn't reach as far as the probe point (radias-3cm), so the tone didn't go off.

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-02-28-CT-2\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170228_CT_tracking_2.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[6,7,8,10,11,14,16,17];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[6,17,14,16,11,10,7,8];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[7,16,17];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 10;
rewardcurrent_Rat139(isession_Rat139,1) = 7;
% 32 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [NaN,0,1,1,1,1,1,0]

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-03-01-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170301_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[6,7,8,9,10,11,14,15,16,17];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[6,17,14,9,11,16,10,15,7,8];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[9,15,16,17];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 7;
rewardcurrent_Rat139(isession_Rat139,1) = 8;
% 46 CA1 place cells, 4 CA1 interneurons
% sign_correct_test = [1,1,1,1,0,0,0,1]

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-03-01-CT-2\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170301_CT_tracking_2.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[6,7,8,9,10,11,14,16];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[6,14,16,11,9,10,7,8];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[9,16];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 8;
rewardcurrent_Rat139(isession_Rat139,1) = 14;
% 41 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [1,1,1,1,1,1,1,1]

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-03-02-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170302_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[6,7,8,9,10,11,14,15,17];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[6,17,14,9,11,10,15,7,8];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[15,17];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 14;
rewardcurrent_Rat139(isession_Rat139,1) = 9;
% 42 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [0,0,0,1,1,1,0,1]

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-03-03-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170303_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[7,8,10,11,15];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[11,10,15,7,8];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[15];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 9;
rewardcurrent_Rat139(isession_Rat139,1) = 11;
% 32 CA1 place cells, 0 CA1 interneurons
% sign_correct_test = [0,1,1,1,0,1,1,0]

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-03-04-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170304_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[6,7,8,9,10,11,14,15,17];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[6,17,14,9,11,10,15,7,8];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[6,17];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 11;
rewardcurrent_Rat139(isession_Rat139,1) = 15;
% 48 CA1 place cells, 0 CA1 interneurons
% sign_correct_test = [1,0,0,1,1,1,0,0]

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-03-05-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170305_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[6,7,8,9,10,11,14,15];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[6,14,9,11,10,15,8,7];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[15,17];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 15;
rewardcurrent_Rat139(isession_Rat139,1) = 11;
% 42 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [1,1,1,0,1,1,1,1]
% There was no tone in sample5 maybe because he passed the location and then turned back, but he got a reward.

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-03-06-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170306_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[6,7,8,9,10,11,14,15];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[6,14,9,11,10,15,8,7];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[15];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 11;
rewardcurrent_Rat139(isession_Rat139,1) = 14;
% 41 CA1 place cells, 1 CA1 interneurons
% sign_correct_test = [0,1,1,1,1,1,1,1]

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-03-08-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170308_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[7,8,9,10,11,14,17];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[17,14,9,11,8,10,7];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[17];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 10;
rewardcurrent_Rat139(isession_Rat139,1) = 12;
% 34 CA1 place cells, 0 CA1 interneurons
% sign_correct_test = [1,1,0,1,1,1,1,1]

isession_Rat139 = isession_Rat139+1;
path_Rat139{isession_Rat139,1} = strcat(DataDir,'Rat139\2017-03-09-CT-1\');
trackdata_Rat139{isession_Rat139,1} = strcat(path_Rat139{isession_Rat139,1},'20170309_CT_tracking_1.mat');
csclist_Rat139_CA1{isession_Rat139,1}=[7,8,9,10,11,14,17];
csclist_Rat139_CA1_sort{isession_Rat139,1}=[17,14,9,11,8,10,7];
csclist_Rat139_CA1_lowfir{isession_Rat139,1}=[14,17];
csclist_Rat139_R2{isession_Rat139,1} = 'HS2R1.ncs';
rewardlast_Rat139(isession_Rat139,1) = 12;
rewardcurrent_Rat139(isession_Rat139,1) = 8;
% 36 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [1,0,0,1,NaN,0,1,0]

isession = isession+isession_Rat139;
Ind_Rat139 = ones(isession_Rat139,1)*139;
Ind_Rat = [Ind_Rat;Ind_Rat139];
pathRats = [pathRats;path_Rat139];
trackdata = [trackdata;trackdata_Rat139];
CSClist_CA1 = [CSClist_CA1;csclist_Rat139_CA1];
CSClist_CA1_sort = [CSClist_CA1_sort;csclist_Rat139_CA1_sort];
CSClist_CA1_lowfir = [CSClist_CA1_lowfir;csclist_Rat139_CA1_lowfir];
CSClist_R2 = [CSClist_R2;csclist_Rat139_R2];
rewardlast = [rewardlast;rewardlast_Rat139];
rewardcurrent = [rewardcurrent;rewardcurrent_Rat139];

%% Rat 150
isession_Rat150 = 0;

isession_Rat150 = isession_Rat150+1;
path_Rat150{isession_Rat150,1} = strcat(DataDir,'Rat150\2017-06-08-CT\');
trackdata_Rat150{isession_Rat150,1} = strcat(path_Rat150{isession_Rat150,1},'20170608_CT_tracking.mat');
csclist_Rat150_CA1{isession_Rat150,1}=[7,9,10,12];
csclist_Rat150_CA1_sort{isession_Rat150,1}=[10,12,7,9];
csclist_Rat150_CA1_lowfir{isession_Rat150,1}=[7];
csclist_Rat150_R2{isession_Rat150,1} = 'HS3R1.ncs';
rewardlast_Rat150(isession_Rat150,1) = 8;
rewardcurrent_Rat150(isession_Rat150,1) = 13;
% 20 CA1 place cells, 1 CA1 interneurons
% sign_correct_test = [0,0,1,1,1,1,1,1]

isession_Rat150 = isession_Rat150+1;
path_Rat150{isession_Rat150,1} = strcat(DataDir,'Rat150\2017-06-10-CT\');
trackdata_Rat150{isession_Rat150,1} = strcat(path_Rat150{isession_Rat150,1},'20170610_CT_tracking.mat');
csclist_Rat150_CA1{isession_Rat150,1}=[7,9,10,12];
csclist_Rat150_CA1_sort{isession_Rat150,1}=[10,12,7,9];
csclist_Rat150_CA1_lowfir{isession_Rat150,1}=[10];
csclist_Rat150_R2{isession_Rat150,1} = 'HS3R1.ncs';
rewardlast_Rat150(isession_Rat150,1) = 9;
rewardcurrent_Rat150(isession_Rat150,1) = 14;
% 24 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [0,0,1,1,1,1,1,1]

isession_Rat150 = isession_Rat150+1;
path_Rat150{isession_Rat150,1} = strcat(DataDir,'Rat150\2017-06-12-CT-1\');
trackdata_Rat150{isession_Rat150,1} = strcat(path_Rat150{isession_Rat150,1},'20170612_CT_tracking_1.mat');
csclist_Rat150_CA1{isession_Rat150,1}=[7,8,9,10,12];
csclist_Rat150_CA1_sort{isession_Rat150,1}=[10,12,7,9,8];
csclist_Rat150_CA1_lowfir{isession_Rat150,1}=[10];
csclist_Rat150_R2{isession_Rat150,1} = 'HS3R1.ncs';
rewardlast_Rat150(isession_Rat150,1) = 10;
rewardcurrent_Rat150(isession_Rat150,1) = 15;
% 41 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [1,1,1,1,1,1,1,1]
% The tracking file was missing.  The file "20170612_CT_tracking_1.mat" was created offline by code "code_create_offline_tracking.m"

isession_Rat150 = isession_Rat150+1;
path_Rat150{isession_Rat150,1} = strcat(DataDir,'Rat150\2017-06-12-CT-2\');
trackdata_Rat150{isession_Rat150,1} = strcat(path_Rat150{isession_Rat150,1},'20170612_CT_tracking_2.mat');
csclist_Rat150_CA1{isession_Rat150,1}=[7,8,9,10,12];
csclist_Rat150_CA1_sort{isession_Rat150,1}=[10,12,7,8,9];
csclist_Rat150_CA1_lowfir{isession_Rat150,1}=[10];
csclist_Rat150_R2{isession_Rat150,1} = 'HS3R1.ncs';
rewardlast_Rat150(isession_Rat150,1) = 15;
rewardcurrent_Rat150(isession_Rat150,1) = 12;
% 36 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [0,0,0,0,1,1,1,1]

isession_Rat150 = isession_Rat150+1;
path_Rat150{isession_Rat150,1} = strcat(DataDir,'Rat150\2017-06-13-CT-1\');
trackdata_Rat150{isession_Rat150,1} = strcat(path_Rat150{isession_Rat150,1},'20170613_CT_tracking_1.mat');
csclist_Rat150_CA1{isession_Rat150,1}=[7,8,9,10,12];
csclist_Rat150_CA1_sort{isession_Rat150,1}=[10,12,7,8,9];
csclist_Rat150_CA1_lowfir{isession_Rat150,1}=[10];
csclist_Rat150_R2{isession_Rat150,1} = 'HS3R1.ncs';
rewardlast_Rat150(isession_Rat150,1) = 12;
rewardcurrent_Rat150(isession_Rat150,1) = 10;
% 36 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [0,1,1,1,1,1,0,1]

isession_Rat150 = isession_Rat150+1;
path_Rat150{isession_Rat150,1} = strcat(DataDir,'Rat150\2017-06-13-CT-2\');
trackdata_Rat150{isession_Rat150,1} = strcat(path_Rat150{isession_Rat150,1},'20170613_CT_tracking_2.mat');
csclist_Rat150_CA1{isession_Rat150,1}=[7,8,9,10,12];
csclist_Rat150_CA1_sort{isession_Rat150,1}=[10,12,7,8,9];
csclist_Rat150_CA1_lowfir{isession_Rat150,1}=[10];
csclist_Rat150_R2{isession_Rat150,1} = 'HS3R1.ncs';
rewardlast_Rat150(isession_Rat150,1) = 10;
rewardcurrent_Rat150(isession_Rat150,1) = 7;
% 26 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [0,1,1,1,1,0,0,0]

isession = isession+isession_Rat150;
Ind_Rat150 = ones(isession_Rat150,1)*150;
Ind_Rat = [Ind_Rat;Ind_Rat150];
pathRats = [pathRats;path_Rat150];
trackdata = [trackdata;trackdata_Rat150];
CSClist_CA1 = [CSClist_CA1;csclist_Rat150_CA1];
CSClist_CA1_sort = [CSClist_CA1_sort;csclist_Rat150_CA1_sort];
CSClist_CA1_lowfir = [CSClist_CA1_lowfir;csclist_Rat150_CA1_lowfir];
CSClist_R2 = [CSClist_R2;csclist_Rat150_R2];
rewardlast = [rewardlast;rewardlast_Rat150];
rewardcurrent = [rewardcurrent;rewardcurrent_Rat150];


%% Rat 149
isession_Rat149 = 0;

isession_Rat149 = isession_Rat149+1;
path_Rat149{isession_Rat149,1} = strcat(DataDir,'Rat149\2017-08-08-CT-1\');
trackdata_Rat149{isession_Rat149,1} = strcat(path_Rat149{isession_Rat149,1},'20170808_CT_tracking_1.mat');
csclist_Rat149_CA1{isession_Rat149,1}=[1,2,3,6,7,8,9,10,14,16,17,18];
csclist_Rat149_CA1_sort{isession_Rat149,1}=[9,10,8,3,6,7,18,14,2,1,16,17];
csclist_Rat149_CA1_lowfir{isession_Rat149,1}=[6,7,9];
csclist_Rat149_R2{isession_Rat149,1} = 'CSC13.ncs';
rewardlast_Rat149(isession_Rat149,1) = 8;
rewardcurrent_Rat149(isession_Rat149,1) = 13;
% 59 CA1 place cells, 4 CA1 interneurons
% sign_correct_test = [0,1,1,1,1,1,1,1]

isession_Rat149 = isession_Rat149+1;
path_Rat149{isession_Rat149,1} = strcat(DataDir,'Rat149\2017-08-08-CT-2\');
trackdata_Rat149{isession_Rat149,1} = strcat(path_Rat149{isession_Rat149,1},'20170808_CT_tracking_2.mat');
csclist_Rat149_CA1{isession_Rat149,1}=[1,2,3,6,7,8,10,14,16,17,18];
csclist_Rat149_CA1_sort{isession_Rat149,1}=[8,10,3,6,7,14,18,2,1,16,17];
csclist_Rat149_CA1_lowfir{isession_Rat149,1}=[7,8];
csclist_Rat149_R2{isession_Rat149,1} = 'CSC13.ncs';
rewardlast_Rat149(isession_Rat149,1) = 13;
rewardcurrent_Rat149(isession_Rat149,1) = 9;
% 55 CA1 place cells, 5 CA1 interneurons
% sign_correct_test = [1,1,0,0,1,1,1,1]

isession_Rat149 = isession_Rat149+1;
path_Rat149{isession_Rat149,1} = strcat(DataDir,'Rat149\2017-08-11-CT-1\');
trackdata_Rat149{isession_Rat149,1} = strcat(path_Rat149{isession_Rat149,1},'20170811_CT_tracking_1.mat');
csclist_Rat149_CA1{isession_Rat149,1}=[1,2,3,6,7,8,9,10,14,16,17,18];
csclist_Rat149_CA1_sort{isession_Rat149,1}=[9,8,3,2,18,10,14,6,7,1,16,17];
csclist_Rat149_CA1_lowfir{isession_Rat149,1}=[2,9];
csclist_Rat149_R2{isession_Rat149,1} = 'CSC13.ncs';
rewardlast_Rat149(isession_Rat149,1) = 13;
rewardcurrent_Rat149(isession_Rat149,1) = 11;
% 61 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [1,0,0,1,1,1,1,0]

isession_Rat149 = isession_Rat149+1;
path_Rat149{isession_Rat149,1} = strcat(DataDir,'Rat149\2017-08-11-CT-2\');
trackdata_Rat149{isession_Rat149,1} = strcat(path_Rat149{isession_Rat149,1},'20170811_CT_tracking_2.mat');
csclist_Rat149_CA1{isession_Rat149,1}=[1,2,3,6,7,8,10,14,16,17,18];
csclist_Rat149_CA1_sort{isession_Rat149,1}=[3,2,8,14,10,18,6,7,1,16,17];
csclist_Rat149_CA1_lowfir{isession_Rat149,1}=[2,3,17];
csclist_Rat149_R2{isession_Rat149,1} = 'CSC13.ncs';
rewardlast_Rat149(isession_Rat149,1) = 11;
rewardcurrent_Rat149(isession_Rat149,1) = 9;
% 55 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [1,1,1,1,0,0,0,0]

isession_Rat149 = isession_Rat149+1;
path_Rat149{isession_Rat149,1} = strcat(DataDir,'Rat149\2017-08-12-CT-1\');
trackdata_Rat149{isession_Rat149,1} = strcat(path_Rat149{isession_Rat149,1},'20170812_CT_tracking_1.mat');
csclist_Rat149_CA1{isession_Rat149,1}=[1,3,6,7,8,9,10,14,16,17,18];
csclist_Rat149_CA1_sort{isession_Rat149,1}=[9,8,18,10,3,6,14,7,1,16,17];
csclist_Rat149_CA1_lowfir{isession_Rat149,1}=[6,9];
csclist_Rat149_R2{isession_Rat149,1} = 'CSC13.ncs';
rewardlast_Rat149(isession_Rat149,1) = 9;
rewardcurrent_Rat149(isession_Rat149,1) = 14;
% 56 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [1,0,1,0,1,1,1,0]

isession_Rat149 = isession_Rat149+1;
path_Rat149{isession_Rat149,1} = strcat(DataDir,'Rat149\2017-08-13-CT-1\');
trackdata_Rat149{isession_Rat149,1} = strcat(path_Rat149{isession_Rat149,1},'20170813_CT_tracking_1.mat');
csclist_Rat149_CA1{isession_Rat149,1}=[1,2,3,6,7,8,9,10,14,16,17,18];
csclist_Rat149_CA1_sort{isession_Rat149,1}=[9,8,10,3,18,14,6,2,7,1,16,17];
csclist_Rat149_CA1_lowfir{isession_Rat149,1}=[2,9];
csclist_Rat149_R2{isession_Rat149,1} = 'CSC13.ncs';
rewardlast_Rat149(isession_Rat149,1) = 14;
rewardcurrent_Rat149(isession_Rat149,1) = 10;
% 60 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [0,1,1,1,1,1,0,0]

isession_Rat149 = isession_Rat149+1;
path_Rat149{isession_Rat149,1} = strcat(DataDir,'Rat149\2017-08-13-CT-2\');
trackdata_Rat149{isession_Rat149,1} = strcat(path_Rat149{isession_Rat149,1},'20170813_CT_tracking_2.mat');
csclist_Rat149_CA1{isession_Rat149,1}=[1,3,6,7,8,9,10,14,16,17,18];
csclist_Rat149_CA1_sort{isession_Rat149,1}=[9,8,3,18,10,14,6,7,1,17,16];
csclist_Rat149_CA1_lowfir{isession_Rat149,1}=[6,9];
csclist_Rat149_R2{isession_Rat149,1} = 'CSC13.ncs';
rewardlast_Rat149(isession_Rat149,1) = 10;
rewardcurrent_Rat149(isession_Rat149,1) = 15;
% 56 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [0,0,0,1,1,1,1,1]

isession_Rat149 = isession_Rat149+1;
path_Rat149{isession_Rat149,1} = strcat(DataDir,'Rat149\2017-08-14-CT-1\');
trackdata_Rat149{isession_Rat149,1} = strcat(path_Rat149{isession_Rat149,1},'20170814_CT_tracking_1.mat');
csclist_Rat149_CA1{isession_Rat149,1}=[1,3,6,7,8,9,10,14,16,17,18];
csclist_Rat149_CA1_sort{isession_Rat149,1}=[9,8,3,10,18,6,14,7,1,17,16];
csclist_Rat149_CA1_lowfir{isession_Rat149,1}=[7,9,17];
csclist_Rat149_R2{isession_Rat149,1} = 'CSC13.ncs';
rewardlast_Rat149(isession_Rat149,1) = 15;
rewardcurrent_Rat149(isession_Rat149,1) = 12;
% 53 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [0,0,1,0,0,0,1,1]

isession_Rat149 = isession_Rat149+1;
path_Rat149{isession_Rat149,1} = strcat(DataDir,'Rat149\2017-08-14-CT-2\');
trackdata_Rat149{isession_Rat149,1} = strcat(path_Rat149{isession_Rat149,1},'20170814_CT_tracking_2.mat');
csclist_Rat149_CA1{isession_Rat149,1} = [1,3,6,7,8,9,10,16,17,18];
csclist_Rat149_CA1_sort{isession_Rat149,1} = [9,8,3,18,10,6,7,1,17,16];
csclist_Rat149_CA1_lowfir{isession_Rat149,1}=[7,9];
csclist_Rat149_R2{isession_Rat149,1} = 'CSC13.ncs';
rewardlast_Rat149(isession_Rat149,1) = 12;
rewardcurrent_Rat149(isession_Rat149,1) = 14;
% 46 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [1,1,1,1,1,1,0,1]

isession_Rat149 = isession_Rat149+1;
path_Rat149{isession_Rat149,1} = strcat(DataDir,'Rat149\2017-08-15-CT-1\');
trackdata_Rat149{isession_Rat149,1} = strcat(path_Rat149{isession_Rat149,1},'20170815_CT_tracking_1.mat');
csclist_Rat149_CA1{isession_Rat149,1} = [1,3,6,7,8,9,10,14,16,17,18];
csclist_Rat149_CA1_sort{isession_Rat149,1} = [9,8,10,18,3,6,14,7,1,16,17];
csclist_Rat149_CA1_lowfir{isession_Rat149,1}=[3,6,17];
csclist_Rat149_R2{isession_Rat149,1} = 'CSC13.ncs';
rewardlast_Rat149(isession_Rat149,1) = 14;
rewardcurrent_Rat149(isession_Rat149,1) = 8;
% 47 CA1 place cells, 4 CA1 interneurons
% sign_correct_test = [0,1,0,1,1,0,1,1]

isession = isession+isession_Rat149;
Ind_Rat149 = ones(isession_Rat149,1)*149;
Ind_Rat = [Ind_Rat;Ind_Rat149];
pathRats = [pathRats;path_Rat149];
trackdata = [trackdata;trackdata_Rat149];
CSClist_CA1 = [CSClist_CA1;csclist_Rat149_CA1];
CSClist_CA1_sort = [CSClist_CA1_sort;csclist_Rat149_CA1_sort];
CSClist_CA1_lowfir = [CSClist_CA1_lowfir;csclist_Rat149_CA1_lowfir];
CSClist_R2 = [CSClist_R2;csclist_Rat149_R2];
rewardlast = [rewardlast;rewardlast_Rat149];
rewardcurrent = [rewardcurrent;rewardcurrent_Rat149];

%% Rat 148
isession_Rat148 = 0;

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-11-15-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171115_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[2,3,4,7,8,9,10,11,12,13,15,16];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[10,12,7,11,8,15,13,2,3,9,4,16];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[7,8,10];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 7;
rewardcurrent_Rat148(isession_Rat148,1) = 13;
% 47 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [1,1,1,1,1,0,1,1]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-11-16-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171116_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[1,2,3,4,7,8,9,10,11,13,14,15,16];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[10,7,11,8,15,14,13,2,9,3,1,16,4];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[8,10,14];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 13;
rewardcurrent_Rat148(isession_Rat148,1) = 8;
% 44 CA1 place cells, 4 CA1 interneurons
% sign_correct_test = [1,0,1,1,1,1,0,0]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-11-17-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171117_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[1,2,3,4,5,7,8,9,10,11,12,13,16];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[10,7,12,11,8,5,13,2,9,3,1,16,4];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[5,8,10,12];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 8;
rewardcurrent_Rat148(isession_Rat148,1) = 11;
% 52 CA1 place cells, 5 CA1 interneurons
% sign_correct_test = [1,1,1,1,1,1,1,0]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-11-18-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171118_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[2,3,7,9,10,11,12,13,14,15,16];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[10,7,12,11,15,14,13,9,3,2,16];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[2,10,12,14];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 11;
rewardcurrent_Rat148(isession_Rat148,1) = 8;
% 37 CA1 place cells, 4 CA1 interneurons
% sign_correct_test = [0,0,0,0,1,1,0,0]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-11-21-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171121_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[1,2,3,4,5,7,8,9,10,11,12,13,15];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[10,12,7,15,5,11,8,13,9,2,3,1,4];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[4,5,12,13];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 8;
rewardcurrent_Rat148(isession_Rat148,1) = 14;
% 43 CA1 place cells, 5 CA1 interneurons
% sign_correct_test = [NaN,1,1,1,1,1,0,1]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-11-22-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171122_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[1,2,3,4,5,7,8,9,10,11,12,13,15];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[7,10,15,12,8,5,11,2,13,9,1,3,4];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[1,12,15];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 14;
rewardcurrent_Rat148(isession_Rat148,1) = 9;
% 51 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [NaN,0,1,0,0,1,1,0]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-11-24-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171124_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[2,3,4,5,7,8,9,10,11,12,13,14,15,16];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[10,14,12,7,15,5,11,8,13,2,9,3,16,4];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[12,14,15];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 9;
rewardcurrent_Rat148(isession_Rat148,1) = 15;
% 54 CA1 place cells, 4 CA1 interneurons
% sign_correct_test = [NaN,1,1,1,0,1,1,1]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-11-29-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171129_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[1,2,3,4,5,7,8,9,11,12,13,16];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[5,7,12,11,8,2,1,13,3,9,4,16];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[1,12];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 12;
rewardcurrent_Rat148(isession_Rat148,1) = 6;
% 51 CA1 place cells, 2 CA1 interneurons
% sign_correct_test = [1,0,1,0,0,1,0,1]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-11-30-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171130_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[1,2,3,4,5,7,8,9,11,12,13,14,15];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[5,14,7,15,12,11,8,13,2,1,3,9,4];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[12,14,15];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 6;
rewardcurrent_Rat148(isession_Rat148,1) = 11;
% 57 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [0,0,0,0,1,0,0,0]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-12-01-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171201_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[2,3,4,5,7,8,9,11,12,13,15];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[5,15,7,12,11,8,13,2,9,3,4];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[5,12,15];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 11;
rewardcurrent_Rat148(isession_Rat148,1) = 7;
% 62 CA1 place cells, 3 CA1 interneurons
% sign_correct_test = [1,0,0,0,0,1,0,0]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-12-02-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171202_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[1,2,3,4,5,7,8,9,10,11,12,13,15,16,18];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[10,7,15,12,5,11,8,2,13,18,3,1,9,16,4];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[5,12,15,18];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 7;
rewardcurrent_Rat148(isession_Rat148,1) = 12;
% 64 CA1 place cells, 4 CA1 interneurons
% sign_correct_test = [0,1,0,1,1,1,1,1]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-12-03-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171203_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[1,2,3,4,5,7,8,9,10,11,12,13,15];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[10,7,5,15,12,11,8,13,2,9,1,3,4];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[1,5,15];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 12;
rewardcurrent_Rat148(isession_Rat148,1) = 9;
% 66 CA1 place cells, 5 CA1 interneurons
% sign_correct_test = [0,0,0,0,1,1,0,1]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-12-09-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171209_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[1,2,3,4,7,8,9,10,11,12,13,15,16];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[10,7,12,15,11,8,13,2,1,16,3,9,4];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[1,4,15];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 14;
rewardcurrent_Rat148(isession_Rat148,1) = 6;
% 66 CA1 place cells, 5 CA1 interneurons
% sign_correct_test = [NaN,1,0,1,1,0,1,1]

isession_Rat148 = isession_Rat148+1;
path_Rat148{isession_Rat148,1} = strcat(DataDir,'Rat148\2017-12-29-CT\');
trackdata_Rat148{isession_Rat148,1} = strcat(path_Rat148{isession_Rat148,1},'20171229_CT_tracking.mat');
csclist_Rat148_CA1{isession_Rat148,1}=[1,2,3,5,7,8,9,12,13,14,15];
csclist_Rat148_CA1_sort{isession_Rat148,1}=[7,5,12,8,15,14,13,2,1,3,9];
csclist_Rat148_CA1_lowfir{isession_Rat148,1}=[1,5,15];
csclist_Rat148_R2{isession_Rat148,1} = 'HS2R1.ncs';
rewardlast_Rat148(isession_Rat148,1) = 10;
rewardcurrent_Rat148(isession_Rat148,1) = 13;
% % 44 CA1 place cells, 5 CA1 interneurons
% % sign_correct_test = [1,1,1,1,1,1,1,1]

isession = isession+isession_Rat148;
Ind_Rat148 = ones(isession_Rat148,1)*148;
Ind_Rat = [Ind_Rat;Ind_Rat148];
pathRats = [pathRats;path_Rat148];
trackdata = [trackdata;trackdata_Rat148];
CSClist_CA1 = [CSClist_CA1;csclist_Rat148_CA1];
CSClist_CA1_sort = [CSClist_CA1_sort;csclist_Rat148_CA1_sort];
CSClist_CA1_lowfir = [CSClist_CA1_lowfir;csclist_Rat148_CA1_lowfir];
CSClist_R2 = [CSClist_R2;csclist_Rat148_R2];
rewardlast = [rewardlast;rewardlast_Rat148];
rewardcurrent = [rewardcurrent;rewardcurrent_Rat148];