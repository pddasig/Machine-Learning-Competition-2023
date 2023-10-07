function  [pickpairs_out,flags_out,a1r,b1r] = sxk_PopUp_PickPairs_V3_OBJ_NPHI(a1,b1,s0_scope,s1_scope,pad_numval)

flags_out = [];

a1log_ori = a1.trend + a1.detail;
b1log_ori = b1.trend + b1.detail;
nrad_bg      = 32;
nrad_scope   = 3;
[a1r,b1r]    = sxk_Effective_Logs_Picks(a1,b1,pad_numval,nrad_bg,nrad_scope);
Apicks = a1r.pick1;
Bpicks = b1r.pick1;
[a0,b0,rtidx0,d0,r0] = sxk_DTW_newobj_Shift_iter_dist(a1log_ori,b1log_ori,Apicks,Bpicks,256,256,pad_numval,[]);

mflen_hor = 5;
wrad_ver  = 2;
[ridx_in] = sxk_DTWroute_enhance(d0,r0,Apicks,mflen_hor,wrad_ver,s1_scope*4);
ridx_in = sortrows(ridx_in,3);
ridx_in = flipud(ridx_in);
Apick_opt = Apicks(ridx_in(:,1));
Bpick_opt = Bpicks(ridx_in(:,2));
Value_opt = ridx_in(:,3);

figure(111)
clf;
imagesc(d0);
hold on;
plot(r0(:,2),r0(:,1),'r.');
plot(ridx_in(:,2),ridx_in(:,1),'w*');
plot(ridx_in(:,2),ridx_in(:,1),'wo');

P11 = [Apicks(1),  Bpicks(1)];
P22 = [Apicks(end),Bpicks(end)];
p1_refine     = 96;
jud_rm = abs(Apick_opt-Apicks(1))<2*p1_refine | abs(Apick_opt-Apicks(end))<2.*p1_refine; % 
Apick_opt(jud_rm,:) = [];
Bpick_opt(jud_rm,:) = [];
pickpairs_out = [P11; [Apick_opt Bpick_opt];P22];

% % % % %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%
% % % tensile_pool  = [0.80:0.01:1.20];
% % % p1_refine     = 96;
% % % a1idxvec = [1:1:length(a1log_ori)];
% % % juda1eff = a1log_ori~=pad_numval & a1log_ori~=pad_numval*2;
% % % a1idx_head = min(a1idxvec(juda1eff));
% % % a1idx_tail = max(a1idxvec(juda1eff));
% % % [pickA4,pickB4] = sxk_RefineCorrection_GR_NPHI_tail(a1,b1,a1idx_tail-560,s1_scope,p1_refine,tensile_pool,pad_numval);
% % % P44 = [pickA4,pickB4(1)];
% % % 
% % % flags      = [];
% % % pickpairs  = [];
% % % for ip = 1:1:length(Apick_opt)
% % %     pickAi = Apick_opt(ip);
% % %     [pickAii,pickBii,flagval] = sxk_RefineCorrection_GR_NPHI(a1,b1,pickAi,s1_scope,p1_refine,tensile_pool,pad_numval);
% % %     
% % %     pickpairs = [pickpairs; [pickAii pickBii(end)]];
% % %     flags     = [flags;flagval];    
% % % end
% % % jud_standout = (flags>=prctile(flags,50) & flags<0.90) & isinf(abs(flags))==0  ;
% % % pickpairs_out    = pickpairs(jud_standout,:);
% % % flags_out        = flags(jud_standout,:);  
% % % 
% % % 
% % % pickpairs_out = [ P11; P22;pickpairs_out; P44];

return;

end
