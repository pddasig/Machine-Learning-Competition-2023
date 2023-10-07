function  [bx,bd,bc,ac,PPairs_s0,PPairs_s1] = sxk_combo_PickPairs_Correction_RD(b0,a0,s0_scope,s1_scope,pad_numval)

jud_alog  = a0.trend ~= pad_numval & a0.depth ~= pad_numval;
aidx_abs = [1:1:length(a0.trend)].';
aidx_eff = aidx_abs(jud_alog);

jud_blog  = b0.trend ~= pad_numval & b0.depth ~= pad_numval;
bidx_abs = [1:1:length(b0.trend)].';
bidx_eff = bidx_abs(jud_blog);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrad_bg      = 32;
nrad_scope   = 3;
[pick3,pick2,a0.pick1] = sub_hierarchy_Peaks(diff(a0.detail),nrad_scope,nrad_bg,4,2,1);
[pick3,pick2,b0.pick1] = sub_hierarchy_Peaks(diff(b0.detail),nrad_scope,nrad_bg,4,2,1);

A0picks   = unique([a0.pick1(:);b0.pick1(:)]);
A0picks   = A0picks(A0picks>1 & A0picks<length(a0.trend));
A0picks(a0.trend(A0picks) == pad_numval) = [];

B0picks   = unique([a0.pick1(:);b0.pick1(:)]);
B0picks   = B0picks(B0picks>1 & B0picks<length(b0.trend));
B0picks(b0.trend(B0picks) == pad_numval) = [];

[AS0,BS0] = sxk_DTW_newobj_Shift(a0.trend,b0.trend,A0picks,B0picks,s0_scope,s1_scope); %%% s1_scope as the max shift 
PPairs_s0 = [AS0,BS0];

figure(4001)
clf;
subplot(211)
hold on;
plot(a0.depth(jud_alog),a0.trend(jud_alog));
plot(b0.depth(jud_blog),b0.trend(jud_blog));
plot(a0.depth(AS0),a0.trend(AS0),'r^');
plot(b0.depth(BS0),b0.trend(BS0),'b^');
xlabel('Depth,ft');
subplot(212);
plot(aidx_eff,a0.trend(jud_alog));
hold on;
plot(bidx_eff+(AS0-BS0),b0.trend(jud_blog));
title(num2str((AS0-BS0)));
xlabel('Index');
title('Shift Correction (S0)');

[a1.trend, b1.trend, a1.depth,b1.depth] = sxk_ShiftLogs(a0.trend,b0.trend,a0.depth,b0.depth,AS0,BS0,pad_numval);
[a1.detail,b1.detail,a1.depth,b1.depth] = sxk_ShiftLogs(a0.detail,b0.detail,a0.depth,b0.depth,AS0,BS0,pad_numval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Step 02 :: Pop-up PickPairs Moderate Scale
%%%%%% shift the logs for this step
nrad_bg      = 32;
nrad_scope   = 3;
[a1r,b1r]    = sxk_Effective_Logs_Picks(a1,b1,pad_numval,nrad_bg,nrad_scope);

[PPairs_s1,flags_s1,a1r,b1r] = sxk_PopUp_PickPairs_V3_OBJ_RD(a1,b1,s0_scope,s1_scope,pad_numval);




% suba0 = a1r.trend + a1r.detail;
% subb0 = b1r.trend + b1r.detail;
% sub_Apicks = a1r.pick1;
% sub_Bpicks = b1r.pick1;
% PPairs_s1    = sxk_PopUp_PickPairs_V2(suba0,subb0,sub_Apicks,sub_Bpicks,s0_scope,s1_scope,pad_numval);

% % % % %%%%% to test and display the pickpairs, use the following function 
% % % % PPairs_s1 = [4 16;
% % % %             188 186;
% % % %              909 909;
% % % %              613 612;
% % % %              1258 1257;
% % % %              2650 2648;
% % % %              2390 2456;
% % % %              3965 3965;
% % % %              4526 4527;
% % % %              4838 4838;
% % % %              4745 4743;
% % % %              4774 4770;
% % % %              4813 4810;
% % % %              5317 5315;
% % % %              5529 5526;
% % % %              5718 5717;
% % % %              5846 5843;
% % % %              5893 5894;
% % % %              6000 5999;
% % % %              6026 6027;
% % % %              6050 6050;
% % % %              6085 6086;
% % % %              6370 6369;
% % % %              6661 6663;
% % % %              6814 6819;
% % % %              6868 6869;
% % % %              6696 6701];         
PPairs_s1 = sortrows(PPairs_s1,1);
sxk_DispCheck_PairPairs(a1r,b1r,PPairs_s1,pad_numval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Step 03 :: Correction ::: Shift & Stretch
[bc.trend, bc.depori,ac.trend, ac.depth] = sxk_pw_Alignment_FreeStyle(PPairs_s1, a1.trend,  b1.trend,  a1.depth, b1.depth, pad_numval);
[bc.detail,bc.depori,ac.detail,ac.depth] = sxk_pw_Alignment_FreeStyle(PPairs_s1, a1.detail, b1.detail, a1.depth, b1.depth, pad_numval);
bc.depth = ac.depth;

bd.depth_pred = sxk_pw_Alignment_depthpred(PPairs_s1,a1.depth, b1.depth, pad_numval); %%%% where the b-log samples locate
bd.trend      = b1.trend;
bd.detail     = b1.detail;
bd.depth      = b1.depth; %%%% the original depth 

jud_b1eff = b1.depth~=pad_numval;
jud_aceff = ac.trend~=pad_numval & ac.depth ~=pad_numval;
jud_bceff = bc.trend~=pad_numval & bc.depth ~=pad_numval;


bc_log = bc.trend; %%% 
b1_log = b1.trend;
ac_log = ac.trend;

figure(4005)
clf;
subplot(211)
plot(ac.depth(jud_aceff),ac_log(jud_aceff),'b-');
hold on;
plot(ac.depth(jud_bceff & jud_aceff),bc_log(jud_bceff & jud_aceff),'r-');
ylabel('Value')
xlabel('Depth,ft');
box on;
title('after alignment');
legend('Reference Log','Target Log');
set(gca,'fontsize',14);
box on;
subplot(212)
plot(b1.depth(jud_b1eff),b1_log(jud_b1eff),'b-');
hold on;
plot(bc.depori(jud_bceff),bc_log(jud_bceff),'r--');
ylabel('Value')
xlabel('Depth,ft');
box on;
title('check target log');
legend('Original','Original-rep');
set(gca,'fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% %%%%%% %%%%%% determine the outputs
jud_b1eff      = b1.depth ~= pad_numval & b1.trend ~= pad_numval;
depth_pred0    = bd.depth_pred(jud_b1eff);
depth_pred1    = sxk_nearestfilling_depth(depth_pred0,pad_numval);

jud_aceff    = ac.depth ~= pad_numval & ac.trend ~= pad_numval;
trend_pred0  = bc.trend(jud_aceff);
detail_pred0 = bc.detail(jud_aceff);
trend_pred1  = sxk_nearestfilling(trend_pred0,pad_numval);
detail_pred1 = sxk_nearestfilling(detail_pred0,pad_numval);

bx.depth_pred = depth_pred1;
bx.log_pred   = trend_pred1 + detail_pred1;

bx.depth      = ac.depth(jud_aceff); %%%% this is only for connections !!!
bx.trend      = trend_pred1; %%%% this is only for connections !!!
bx.detail     = detail_pred1; %%%% this is only for connections !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% the following is used for checking the final results
jud_effbd      = bd.trend~=pad_numval & bd.depth~=pad_numval;
depth_pred_bd  = bd.depth_pred(jud_effbd);
trend_pred_bd  = bd.trend(jud_effbd);
detail_pred_bd = bd.detail(jud_effbd);

figure(4006)
clf;
subplot(211)
plot(a0.depth(jud_alog),a0.trend(jud_alog)+a0.detail(jud_alog),'b-');
hold on;
plot(a0.depth(jud_alog),trend_pred1+detail_pred1,'r-');
subplot(212)
plot(a0.depth(jud_alog),trend_pred1+detail_pred1,'r-');
hold on;
plot(depth_pred_bd,trend_pred_bd+detail_pred_bd,'g--');

return;

end
