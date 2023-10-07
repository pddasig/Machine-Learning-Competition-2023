function  [pickA0,pickB0,pickB0_retune] = sub_Step01_FindA0B0(loga0_trend,logb0_trend,PickSets,L0_search_nrad,L0_retune_nrad,tensile_pool)

% L0_search_nrad  = 8000;
% L0_retune_nrad  = 256;
% tensile_pool    = [0.90:0.02:1.10];  %%% if you wanna test tensile factors
% L0_retune_step  = 4;                 %%% only useful for point-wise retuning

Apidx_L1 = PickSets.Apidx_L1;
Apidx_L2 = PickSets.Apidx_L2;
Apidx_L3 = PickSets.Apidx_L3;
Bpidx_L1 = PickSets.Bpidx_L1;
Bpidx_L2 = PickSets.Bpidx_L2;
Bpidx_L3 = PickSets.Bpidx_L3;

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Step 1-1: Find the First Pair(L0)
loga00    = loga0_trend;
logb00    = logb0_trend;
ndepth_aa = length(loga00);
ndepth_bb = length(logb00);

[~,idx] = min(abs(Apidx_L3-0.5.*ndepth_aa));
pickA0  = Apidx_L3(idx);
pickB0s = Bpidx_L2(Bpidx_L2>=(pickA0-L0_search_nrad) & Bpidx_L2<=(pickA0+L0_search_nrad));

minrad = min([512 floor(0.5.*ndepth_aa)]);
if pickA0 < minrad || pickA0 > (ndepth_aa-minrad)
    [~,idx] = min(abs(Apidx_L2-0.5.*ndepth_aa));
    pickA0  = Apidx_L2(idx);
    pickB0s = Bpidx_L1(Bpidx_L1>=(pickA0-L0_search_nrad) & Bpidx_L1<=(pickA0+L0_search_nrad));
end

%%%% how many pairs to be tested
obj2d_later   = sub_FindB0_ObjFun_LaterPart_Search(loga00,logb00,pickA0,pickB0s,tensile_pool);
obj2d_front   = sub_FindB0_ObjFun_FrontPart_Search(loga00,logb00,pickA0,pickB0s,tensile_pool);
[~,tmpcidx] = max(max(obj2d_later+obj2d_front,[],1));
pickB0      = pickB0s(tmpcidx);

logbshift =  pickA0 - pickB0;
figure(10001)
clf;
subplot(211);
plot([1:1:ndepth_aa],loga00,'b-');
hold on;
plot(pickA0,loga00(pickA0),'r*');
plot([1:1:ndepth_bb],logb00,'m-');
subplot(212);
plot([1:1:ndepth_aa],loga00,'b-');
hold on;
plot([1:1:ndepth_bb]+logbshift,logb00,'m-');
plot(pickA0,loga00(pickA0),'r*');
plot(pickA0,loga00(pickA0),'ro');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Step 1-2: Retune the First Pair(L0)
nrad_retune = L0_retune_nrad*2; %% log segments for comparsion

idxaa = pickA0 + [-nrad_retune:1:nrad_retune];
idxbb = pickB0 + [-nrad_retune:1:nrad_retune];
judee = idxaa>=1 & idxaa<=ndepth_aa & idxbb>=1 & idxbb<=ndepth_bb;
idxaae = idxaa(judee);
idxbbe = idxbb(judee);
[~,cenidx] = min(abs(idxaae-pickA0));

lw_retune = max([min(idxaae-pickA0) min(idxbbe-pickB0) -L0_retune_nrad]);
up_retune = min([max(idxaae-pickA0) max(idxbbe-pickB0) +L0_retune_nrad]);
% % picksRE   = cenidx + [lw_retune:L0_retune_step:up_retune];  %%% point-wise
picksRE_abs   = Bpidx_L1(Bpidx_L1>=(pickB0+lw_retune) & Bpidx_L1<=(pickB0 +up_retune)); %%% peak-wise
picksRE_rel   = cenidx + picksRE_abs-pickB0;

logaa_tar = loga00(idxaae);
logbb_tar = logb00(idxbbe);
obj2d_laterRE   = sub_FindB0_ObjFun_LaterPart_retune(logaa_tar,logbb_tar,cenidx,picksRE_rel,tensile_pool);
obj2d_frontRE   = sub_FindB0_ObjFun_FrontPart_retune(logaa_tar,logbb_tar,cenidx,picksRE_rel,tensile_pool);
[~,tmpcidx] = max(max(obj2d_laterRE+obj2d_frontRE,[],1));

pickB0_retune = picksRE_abs(tmpcidx);
logbshift_retune     =  pickA0 - pickB0_retune;

figure(10002)
clf;
subplot(211);
plot([1:1:ndepth_aa],loga00,'b-');
hold on;
plot(pickA0,loga00(pickA0),'r*');
plot([1:1:ndepth_bb],logb00,'m-');
subplot(212);
plot([1:1:ndepth_aa],loga00,'b-');
hold on;
plot([1:1:ndepth_bb]+logbshift_retune,logb00,'m-');
plot(pickA0,loga00(pickA0),'r*');
plot(pickA0,loga00(pickA0),'ro');

return;

end
