function  [blog_out,bdep_ori,alog_out,adep_out] = sxk_pw_Alignment_FreeStyle_v2(pickpairs,suba0,subb0,dep_suba0,dep_subb0,pad_numval)

pickpairs = sortrows(pickpairs,1);
npair = size(pickpairs,1);

ndepa0 = length(suba0);
ndepb0 = length(subb0);

blog_aligned = ones(ndepa0,1).*pad_numval;
bdep_aligned = ones(ndepa0,1).*pad_numval;
alog_aligned = ones(ndepa0,1).*pad_numval;
adep_aligned = ones(ndepa0,1).*pad_numval;
for ip = 1:1:npair-1      
    a1   = [pickpairs(ip,1):1:pickpairs(ip+1,1)];    %%% reference index
    b1   = [pickpairs(ip,2):1:pickpairs(ip+1,2)];    %%% target index 
    
    alog = suba0(a1); 
    adep = dep_suba0(a1);    
    blog = subb0(b1);
%     bdep = dep_subb0(b1);   
    bdep = b1;
    
    [blog_al,bdep_al] = sub_allignB2A_all(alog,adep,blog,bdep);
% % %     tlog_al = sub_allignB2A(rlog,tlog);

    blog_aligned(a1) = blog_al;
    bdep_aligned(a1) = bdep_al; %%% the depth is original depth: where the data coming from
    
    alog_aligned(a1) = alog;
    adep_aligned(a1) = adep;    
end




judeffbody = alog_aligned~=pad_numval;
blog_aligned = blog_aligned(judeffbody);
bdep_aligned = bdep_aligned(judeffbody); %%%% this is different 
alog_aligned = alog_aligned(judeffbody);
adep_aligned = adep_aligned(judeffbody);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% process the head
na_head = pickpairs(1,1)-1;
nb_head = pickpairs(1,2)-1;
nc_head = max([na_head nb_head]);
head_tmpl = ones(nc_head,1).*pad_numval;

blog_head = head_tmpl;
blog_head(end-nb_head+1:end) = subb0(1:1:nb_head);
bdep_head = head_tmpl;
bdep_head(end-nb_head+1:end) = [1:1:nb_head]; %%% %%% changed
alog_head = head_tmpl;
alog_head(end-na_head+1:end) = suba0(1:1:na_head);
adep_head = head_tmpl;
adep_head(end-na_head+1:end) = dep_suba0(1:1:na_head);

%%%%%% process the tail
na_tail = ndepa0-pickpairs(end,1);
nb_tail = ndepb0-pickpairs(end,2);
nc_tail = max([na_tail nb_tail]);
tail_tmpl = ones(nc_tail,1).*pad_numval;

blog_tail = tail_tmpl;
blog_tail(1:nb_tail) = subb0(pickpairs(end,2)+1:1:end);
bdep_tail = tail_tmpl;
bdep_tail(1:nb_tail) = [pickpairs(end,2)+1:1:length(dep_subb0)]; %%% changed

alog_tail = tail_tmpl;
alog_tail(1:na_tail) = suba0(pickpairs(end,1)+1:1:end);
adep_tail = tail_tmpl;
adep_tail(1:na_tail) = dep_suba0(pickpairs(end,1)+1:1:end);

blog_out = [blog_head;blog_aligned;blog_tail];
bdep_ori = [bdep_head;bdep_aligned;bdep_tail];

alog_out = [alog_head;alog_aligned;alog_tail];
adep_out = [adep_head;adep_aligned;adep_tail];

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%% check the allignment outputs
% % % jud_aeff = alog_out~=pad_numval;
% % % jud_beff = blog_out~=pad_numval;
% % % figure(4004)
% % % clf;
% % % subplot(211)
% % % plot(adep_out(jud_aeff),alog_out(jud_aeff),'b-');
% % % hold on;
% % % plot(adep_out(jud_aeff & jud_beff),blog_out(jud_aeff & jud_beff),'r-');
% % % ylabel('Value')
% % % xlabel('Depth,ft');
% % % box on;
% % % title('after alignment');
% % % legend('Reference Log','Target Log');
% % % set(gca,'fontsize',14)
% % % 
% % % jud_b0eff = subb0~=pad_numval;
% % % 
% % % subplot(212)
% % % plot(dep_subb0(jud_b0eff),subb0(jud_b0eff),'b-');
% % % hold on;
% % % plot(bdep_ori(jud_beff),blog_out(jud_beff),'r--');
% % % ylabel('Value')
% % % xlabel('Depth,ft');
% % % box on;
% % % title('check target log');
% % % legend('Original','Original-rep');
% % % set(gca,'fontsize',14)

return;

end
