function  [pickpairs_out,flags_out,a1r,b1r] = sxk_PopUp_PickPairs_V3_OBJ_RHOB(a1,b1,s0_scope,s1_scope,pad_numval)


a1log_ori = a1.trend + a1.detail;
b1log_ori = b1.trend + b1.detail;
nrad_bg      = 32;
nrad_scope   = 3;
[a1r,b1r]    = sxk_Effective_Logs_Picks(a1,b1,pad_numval,nrad_bg,nrad_scope);
Apicks = a1r.pick1;
Bpicks = b1r.pick1;
[a0,b0,rtidx0,d0,r0] = sxk_DTW_newobj_Shift_iter(a1log_ori,b1log_ori,Apicks,Bpicks,2.*s1_scope,round(s1_scope),pad_numval,[]);

mflen_hor = 5;
wrad_ver  = 2;
[ridx_in] = sxk_DTWroute_enhance(d0,r0,Apicks,mflen_hor,wrad_ver,s1_scope./2);
Apick_opt = Apicks(ridx_in(:,1));
figure(111)
clf;
imagesc(d0);
hold on;
plot(r0(:,2),r0(:,1),'r.');
plot(ridx_in(:,2),ridx_in(:,1),'w*');
plot(ridx_in(:,2),ridx_in(:,1),'wo');


%%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%
tensile_pool  = [0.85:0.01:1.15];
p1_refine     = 64;

a1idxvec = [1:1:length(a1log_ori)];
juda1eff = a1log_ori~=pad_numval & a1log_ori~=pad_numval*2;
a1idx_head = min(a1idxvec(juda1eff));
a1idx_tail = max(a1idxvec(juda1eff));
[pickA1,pickB1,flag1] = sxk_RefineCorrection_NPHI_RHOB_head(a1,b1,a1idx_head,s1_scope,p1_refine,tensile_pool,pad_numval);
[pickA2,pickB2,flag2] = sxk_RefineCorrection_NPHI_RHOB_tail(a1,b1,a1idx_tail,s1_scope,p1_refine,tensile_pool,pad_numval);
    
jud_rm = abs(Apick_opt-pickA1)<3.*p1_refine | abs(Apick_opt-pickA2)<3.*p1_refine;
Apick_opt(jud_rm,:) = [];

flags      = [];
pickpairs  = [];
for ip = 1:1:length(Apick_opt)
    pickAi = Apick_opt(ip);
    [pickAii,pickBii,flagval] = sxk_RefineCorrection_NPHI_RHOB(a1,b1,pickAi,s1_scope,p1_refine,tensile_pool,pad_numval);
    
    pickpairs = [pickpairs; [pickAii pickBii(end)]];
    flags     = [flags;flagval];    
end
jud_standout = (flags>=0.30 & flags<10) & isinf(abs(flags))==0  ;
pickpairs_out    = pickpairs(jud_standout,:);
flags_out        = flags(jud_standout,:);  

P11 = [pickA1,pickB1(end)];
P22 = [pickA2,pickB2(end)];
pickpairs_out  = [pickpairs_out ;  P11;  P22];
flags_out      = [flags_out;       flag1(end);  flag2(end)];
mat_out       = sortrows([pickpairs_out flags_out],1);
pickpairs_out = mat_out(:,1:2);
flags_out     = mat_out(:,3);
 

return;

end
