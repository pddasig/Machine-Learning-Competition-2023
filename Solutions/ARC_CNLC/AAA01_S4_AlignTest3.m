
clear all;
clc;

mflenSEP  = 41;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% step 01 :: processing the test dataset
Headline_str = 'DEPT,GR,RHOB,NPHI,RD';
load('Step01_LoadPreprocessing.mat','-mat');
test_csv_file = 'test_well_03.csv';

test = sub_read_the_Train_Data(test_csv_file,Headline_str);   
test(:,5) = log10(test(:,5));
test_scal = test;
test_scal(:,2:end) = (test(:,2:end) - Mean_all)./STD_all;
Test.depth = test_scal(:,1);
Test.gr    = test_scal(:,2);
Test.rhob  = test_scal(:,3);
Test.nphi  = test_scal(:,4);
Test.rd    = test_scal(:,5);

Test.gr_trend   = medfilt2(Test.gr,[mflenSEP 1]);
Test.rhob_trend = medfilt2(Test.rhob,[mflenSEP 1]);
Test.nphi_trend = medfilt2(Test.nphi,[mflenSEP 1]);
Test.rd_trend   = medfilt2(Test.rd,[mflenSEP 1]);
Test.gr_detail       = Test.gr-Test.gr_trend;
Test.rd_detail       = Test.rd-Test.rd_trend;
Test.nphi_detail     = Test.nphi-Test.nphi_trend;
Test.rhob_detail     = Test.rhob-Test.rhob_trend;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% step 03 :: Auto-Correction of Gamma
a0.detail = Test.gr_detail;  %%% reference log 
a0.trend  = Test.gr_trend;   %%% reference log 
a0.depth  = Test.depth;

% % b0.detail = Test.rd_detail; %%% to-be-alligned log
% % b0.trend  = Test.rd_trend;  %%% to-be-alligned log
% % b0.depth  = Test.depth;
% % s0_scope   = 1024;
% % s1_scope   = 256;
% % pad_numval = 999.25;
% % [bx_rd,bd_rd,bc_rd,ac_rd,PP0_rd,PP1_rd] = sxk_combo_PickPairs_Correction_RD(b0,a0,s0_scope,s1_scope,pad_numval);

load('test3_outcome.mat','-mat');

b0.detail = Test.nphi_detail; %%% to-be-alligned log
b0.trend  = Test.nphi_trend;  %%% to-be-alligned log
b0.depth  = Test.depth;
s0_scope   = 512;
s1_scope   = 128;
pad_numval = 999.25;
[bx_nphi0,bd_nphi,bc_nphi,ac_nphi,PP0_nphi,PP1_nphi] = sxk_combo_PickPairs_Correction_NPHI(b0,a0,s0_scope,s1_scope,pad_numval);
s0_scope   = 128;
s1_scope   = 16;
[bx_nphi,bd_nphi,bc_nphi,ac_nphi,PP0_nphi,PP1_nphi] = sxk_combo_PickPairs_Correction_NPHI(bx_nphi0,a0,s0_scope,s1_scope,pad_numval);


b0.detail  = Test.rhob_detail; %%% to-be-alligned log
b0.trend   = Test.rhob_trend;  %%% to-be-alligned log
b0.depth   = Test.depth;
s0_scope   = 512;
s1_scope   = 128;
pad_numval = 999.25;
[bx_rhob,bd_rhob,bc_rhob,ac_rhob,PP0_rhob,PP1_rhob] = sxk_combo_PickPairs_Correction_RHOB(b0,bx_nphi,s0_scope,s1_scope,pad_numval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhob_dept_pred = bx_rhob.depth_pred;
nphi_dept_pred = bx_nphi.depth_pred;
rd_dept_pred   = bx_rd.depth_pred;

rhob_pred = bx_rhob.log_pred*STD_all(2) + Mean_all(2);
nphi_pred = bx_nphi.log_pred*STD_all(3) + Mean_all(3);
rd_pred = bx_rd.log_pred*STD_all(4) + Mean_all(4);
rd_pred = 10.^(rd_pred);


save test3_outcome_sub34;

final_result = [rhob_pred nphi_pred rd_pred rhob_dept_pred nphi_dept_pred rd_dept_pred];
writematrix(final_result, 'test3-add.csv');



