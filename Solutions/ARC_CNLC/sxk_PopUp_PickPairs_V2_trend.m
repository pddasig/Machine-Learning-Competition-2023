function  pickpairs = sxk_PopUp_PickPairs_V2_trend(a1,b1,s0_scope,s1_scope,pad_numval)

a1log_ori = a1.trend ;
b1log_ori = b1.trend ;

nrad_bg      = 32;
nrad_scope   = 3;
[a1r,b1r]    = sxk_Effective_Logs_Picks(a1,b1,pad_numval,nrad_bg,nrad_scope);
Apicks = a1r.pick1;
Bpicks = b1r.pick1;

a1.idxvec = [1:1:length(a1log_ori)].';
judeff    = a1.trend~=pad_numval; %%%% based on a0. %%%% & b1.trend~=pad_numval
a1.hidx0     = min(a1.idxvec(judeff));
a1.tidx0     = max(a1.idxvec(judeff));
ht_tmpl0  = [a1.hidx0 a1.hidx0; a1.tidx0 a1.tidx0];

b1.idxvec = [1:1:length(b1log_ori)].';
judeff    = b1.trend~=pad_numval; %%%% based on a0. %%%% & b1.trend~=pad_numval
b1.hidx0     = min(b1.idxvec(judeff));
b1.tidx0     = max(b1.idxvec(judeff));


[a0,b0,rtidx0,d0,r0] = sxk_DTW_newobj_Shift_iter(a1log_ori,b1log_ori,Apicks,Bpicks,2.*s1_scope,round(s1_scope),pad_numval,[]);

jud_clearA = (Apicks >= (a0-s1_scope) & Apicks <= (a0+s1_scope));
jud_clearB = (Bpicks >= (b0-s1_scope) & Bpicks <= (b0+s1_scope));
Apicks(jud_clearA) = [];
Bpicks(jud_clearB) = [];

a0_idx   = [1:1:length(a1log_ori)];
b0_idx   = [1:1:length(b1log_ori)];
judeffa0 = a1log_ori < pad_numval;
judeffb0 = b1log_ori < pad_numval;
figure(4002)
clf;
set(gcf,'position',[-2235 699 2498 420]);
figure(4002)
subplot(1,10,1:2)
cla;
imagesc(d0);
hold on;
plot(r0(:,2),r0(:,1),'r.');
plot(rtidx0(2),rtidx0(1),'w*','markersize',8);
plot(rtidx0(2),rtidx0(1),'wo','markersize',8);
% % caxis([-0.5 0.5])
subplot(1,10,3:10)
hold on
plot(a0_idx(judeffa0),a1log_ori(judeffa0));
plot(b0_idx(judeffb0),b1log_ori(judeffb0));
plot(a0,a1log_ori(a0),'r^','markersize',8);
plot(b0,b1log_ori(b0),'b^','markersize',8);
grid on;
grid minor;
box on;
xlabel('depth index')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntimes = max([round(length(a1log_ori)./(1.*s1_scope)),16]);

a1_iter = a1;
b1_iter = b1;
a1_hidx0_iter = a1.hidx0;
a1_tidx0_iter = a1.tidx0;
b1_hidx0_iter = b1.hidx0;
b1_tidx0_iter = b1.tidx0;
b1chkidx_iter = round(b0_idx);
pickpairs = [a0,b0];

for it = 1:1:ntimes      
    if length(Apicks)<s1_scope || length(Bpicks)<s1_scope
        break;
    end    
    
    ailog_iter = a1_iter.trend ;
    bilog_iter = b1_iter.trend ;    
    [ai,bi,rtidx,di,ri] = sxk_DTW_newobj_Shift_iter(ailog_iter,bilog_iter,Apicks,Bpicks,2.*s1_scope,round(s1_scope),pad_numval,[]);
   
    jud_clearA = ~(Apicks >= (ai-1.*s1_scope) & Apicks <= (ai+1.*s1_scope));
    jud_clearB = ~(Bpicks >= (bi-1.*s1_scope) & Bpicks <= (bi+1.*s1_scope));
    Apicks = Apicks(jud_clearA);
    Bpicks = Bpicks(jud_clearB);    
    jud_eff  = all(sign(pickpairs(:,1)-ai) == sign(pickpairs(:,2)-bi)) & ailog_iter(ai)~=pad_numval & bilog_iter(bi)~=pad_numval;    
    if ~jud_eff 
        continue;
    end 
    
    jud_goodenough = abs(ai-bi)<(0.5.*s1_scope) & ai>a1_hidx0_iter & ai<a1_tidx0_iter ...
                                                & bi>b1_hidx0_iter & bi<b1_tidx0_iter; %%% to correction and update the iteration 
    
    if jud_goodenough==1     
        a1i_idxvec = [1:1:length(ailog_iter)].';        
        judeff     = a1_iter.trend~=pad_numval; 
        hidx_ii    = min(a1i_idxvec(judeff));
        a1p_now    = ai - (hidx_ii-a1.hidx0);   %%%% find the index on original a1            
        b1p_now    = b1chkidx_iter(bi);      %%%% find the index on original b1      
        
        pickpairs    = sortrows([pickpairs;[a1p_now b1p_now]],1);         
        picks_corr   = sortrows(pickpairs,1); %%% sortrows([ht_tmpl0;pickpairs],1)
       
        figure(4002)
        subplot(1,10,1:2)
        cla;
        imagesc(di);
        hold on;
        plot(ri(:,2),ri(:,1),'r.');
        plot(rtidx(2),rtidx(1),'w*','markersize',8);
        plot(rtidx(2),rtidx(1),'wo','markersize',8);
        xlim([1 size(di,2)])
        ylim([1 size(di,1)])
        subplot(1,10,3:10)
        plot(a1p_now,a1log_ori(a1p_now),'r*','markersize',8);
        plot(a1p_now,a1log_ori(a1p_now),'ro','markersize',8);
        hold on
        plot(b1p_now,b1log_ori(b1p_now),'b*','markersize',8);
        plot(b1p_now,b1log_ori(b1p_now),'bo','markersize',8); 
       
        [b1i.trend,b1chkidx_iter, a1i.trend,a1i.depth] = sxk_pw_Alignment_FreeStyle_v2(picks_corr,a1.trend,b1.trend,a1.depth,b1.depth,pad_numval);
        [b1i.detail,b1chkidx_iter,a1i.detail,a1i.depth] = sxk_pw_Alignment_FreeStyle_v2(picks_corr,a1.detail,b1.detail,a1.depth,b1.depth,pad_numval);
        b1i.depth = a1i.depth;               
        b1chkidx_iter = round(b1chkidx_iter);
        b1chkidx_iter(b1chkidx_iter<1) = 1;
        b1chkidx_iter(b1chkidx_iter>length(b1log_ori)) = length(b1log_ori);
       
        [a1ir,b1ir]    = sxk_Effective_Logs_Picks(a1i,b1i,pad_numval,nrad_bg,nrad_scope);
        a1_iter        = a1ir;
        b1_iter        = b1ir;
        Apicks = a1ir.pick1;
        Bpicks = b1ir.pick1;
        
        for ip = 1:1:size(pickpairs,1)
            aip = pickpairs(ip,1);
            bip = pickpairs(ip,1);
            
            jud_clearA = ~(Apicks >= (aip-1.*s1_scope) & Apicks <= (aip+1.*s1_scope));
            jud_clearB = ~(Bpicks >= (bip-1.*s1_scope) & Bpicks <= (bip+1.*s1_scope));
            Apicks = Apicks(jud_clearA);
            Bpicks = Bpicks(jud_clearB);
        end      
        
        a1idxvec = [1:1:length(a1_iter.trend)].';
        judeff    = a1.trend~=pad_numval; %%%% based on a0. %%%% & b1.trend~=pad_numval
        a1_hidx0_iter     = min(a1idxvec(judeff));
        a1_tidx0_iter     = max(a1idxvec(judeff));      
        b1idxvec     = [1:1:length(b1_iter.trend)].';
        judeff       = b1.trend~=pad_numval; %%%% based on a0. %%%% & b1.trend~=pad_numval
        b1_hidx0_iter  = min(b1idxvec(judeff));
        b1_tidx0_iter  = max(b1idxvec(judeff));        
    end
      
end
pickpairs = sortrows(pickpairs,1);

return;

end
