function  [pair_abs] = sxk_RouteNodes(route0,min_dist,acc,bcc)

np = size(route0,1);

idx4front = repmat([1:1:np-1].',[1 2]);
idx4later = repmat([2:1:np].',[1 2]);

jud_flat = diff(route0,1,1)==0;
tmp_zero  = [1==0 1==0];
jud_hnode = diff([tmp_zero;jud_flat],1,1)==1; %%% head node for a flat segment
hnode_c1  = idx4front(jud_hnode(:,1),1);  %% head index for the first curve
hnode_c2  = idx4front(jud_hnode(:,2),2);  %% tail index for the first curve

jud_tnode = diff([jud_flat;tmp_zero],1,1)==(-1); %%% head node for a flat segment
tnode_c1  = idx4later(jud_tnode(:,1),1);  %% head index for the first curve
tnode_c2  = idx4later(jud_tnode(:,2),2);  %% tail index for the first curve

acc_hts = [hnode_c1(:).'; tnode_c1(:).'];  %%% index for a prolonged log curve
bcc_hts = [hnode_c2(:).'; tnode_c2(:).'];  %%% index for a prolonged log curve
jud_begin = any(acc_hts==1,1);
jud_close = any(acc_hts==np,1);
acc_hts = acc_hts(:,(~jud_begin) & (~jud_close));
jud_begin2 = any(bcc_hts==1,1);
jud_close2 = any(bcc_hts==np,1);
bcc_hts = bcc_hts(:,(~jud_begin2) & (~jud_close2));

acc_hts1 = [];
for ia = 1:1:size(acc_hts,2)-1
    tmp = [acc_hts(2,ia);acc_hts(1,ia+1)];    
    acc_hts1 = [acc_hts1 tmp];
end

bcc_hts1 = [];
for ib = 1:1:size(bcc_hts,2)-1
    tmp = [bcc_hts(2,ib);bcc_hts(1,ib+1)];    
    bcc_hts1 = [bcc_hts1 tmp];
end

acc_hts = acc_hts1;
bcc_hts = bcc_hts1;
% % % acc_hts_abs = route0(acc_hts,1);
% % % bcc_hts_abs = route0(bcc_hts,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1_acc    = acc_hts(1);
p2_acc    = [];
p1_bcc    = bcc_hts(1);
p2_bcc    = [];

segidx_acc = [];
segidx_bcc = [];
for ip   = 2:1:min([size(acc_hts,2) size(bcc_hts,2)])    
    
    p2_acc_test = max([acc_hts(2,ip) bcc_hts(2,ip)]);
    p2_bcc_test = max([acc_hts(2,ip) bcc_hts(2,ip)]);    
    pt_acc = [p1_acc:1:p2_acc_test];
    pt_bcc = [p1_bcc:1:p2_bcc_test];    
    eff_anum = sum(diff(acc(pt_acc))~=0);
    eff_bnum = sum(diff(bcc(pt_bcc))~=0);    

    if all([eff_anum eff_bnum]>=min_dist) == 1
        p2_acc = p2_acc_test;
        segidx_acc = [segidx_acc; [p1_acc p2_acc]];
        p1_acc = p2_acc + 1;
        p2_acc = [];
                
        p2_bcc = p2_bcc_test;
        segidx_bcc = [segidx_bcc; [p1_bcc p2_bcc]];
        p1_bcc = p2_bcc + 1;
        p2_bcc = [];        
    end
end

if isempty(p2_acc)==1 || p2_acc~=acc_hts(end)           
    segidx_acc(end)   = acc_hts(end);  %%% index for prolonged logs
end
if isempty(p2_bcc)==1 || p2_bcc~=bcc_hts(end)         
    segidx_bcc(end)     = bcc_hts(end);  %%% index for prolonged logs
end

segidx_acc_abs = [route0(segidx_acc(:,1),1) route0(segidx_acc(:,2),1)];
segidx_bcc_abs = [route0(segidx_bcc(:,1),2) route0(segidx_bcc(:,2),2)];

pair_abs = sortrows([segidx_acc_abs(:) segidx_bcc_abs(:)]);
jud_near = [0==1;any(diff(pair_abs,1,1)<=1,2)];
pair_abs = pair_abs(~jud_near,:);

% % figure(111)
% % clf;
% % subplot(211)
% % plot(acc);
% % hold on;
% % plot(acc_hts(1,:),acc(acc_hts(1,:)),'r.');
% % plot(acc_hts(2,:),acc(acc_hts(2,:)),'r>');
% % subplot(212)
% % plot(bcc);
% % hold on;
% % plot(bcc_hts(1,:),bcc(bcc_hts(1,:)),'r.');
% % plot(bcc_hts(2,:),bcc(bcc_hts(2,:)),'r>');


return;

end
