function  pair_dtw = sxk_additional_DTW(bc,ac,pad_numval)

jud_ac = ac.trend~=pad_numval;
jud_bc = bc.trend~=pad_numval;
jud_comm  = jud_ac & jud_bc;

idxvec_all  = [1:1:length(ac.trend)];
idxvec_comm = idxvec_all(jud_comm);
hidx_comm   = min(idxvec_comm)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ac_act       = ac.trend(jud_comm);
bc_act       = bc.trend(jud_comm);
ac_depth_act = ac.depth(jud_comm);
bc_depth_act = bc.depth(jud_comm);

% load test;
[route0,acc,bcc]   = dtw_dist(ac_act(:),bc_act(:),0); %%% most common dtw form

min_dist  = 3;
pair_rel  = sxk_RouteNodes(route0,min_dist,acc,bcc);
pair_dtw  = pair_rel + hidx_comm;


acc_depth = ac_depth_act(route0(:,1));
bcc_depth = bc_depth_act(route0(:,2));
figure(111)
clf;
subplot(211)
plot(ac_depth_act,ac_act);
hold on;
plot(bc_depth_act,bc_act);
subplot(212)
plot(acc_depth,acc,'b-'); %%% repeat the original curve
hold on
plot(bcc_depth,bcc,'r-');

figure(222)
clf;
subplot(211);
plot(ac_act);
hold on;
plot(pair_rel(:,1),ac_act(pair_rel(:,1)),'r*');
subplot(212)
plot(bc_act);
hold on;
plot(pair_rel(:,2),bc_act(pair_rel(:,2)),'b*');

return;

end
