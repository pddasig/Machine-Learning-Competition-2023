function  [pickAi,pickBi,isend] = sub_Find_L1_peakpairs_DW(loga0_trend,logb0_trend,PickSets,pickA0,pickB0,L1_search_nrad,proc_warplen)

Apidx_L1 = PickSets.Apidx_L1;
Apidx_L2 = PickSets.Apidx_L2;
Apidx_L3 = PickSets.Apidx_L3;
Bpidx_L1 = PickSets.Bpidx_L1;
Bpidx_L2 = PickSets.Bpidx_L2;
Bpidx_L3 = PickSets.Bpidx_L3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pickAis    = Apidx_L1(Apidx_L1 >= (pickA0+L1_search_nrad));
[~,idx]    = min(abs(pickAis-(pickA0+proc_warplen)));
pickAi     =  pickAis(idx);
pickBi_ref = pickB0 + (pickAi-pickA0);
pickBis    = Bpidx_L1(Bpidx_L1 >= (pickBi_ref-L1_search_nrad) & Bpidx_L1 <= (pickBi_ref+L1_search_nrad));

if isempty(pickBis)==1
    pickAi = [];
    pickBi = [];
    isend  = 1;
    return;
end

loga00 = loga0_trend(pickA0:pickAi);
loga11 = (loga00-mean(loga00))./std(loga00); %%%% useful for short segments
% % % loga11 = loga00;

if pickAi==max(Apidx_L1) || length(loga11)<=64
    isend = 1;
else
    isend = 0;
end

ntry = length(pickBis);
obj_ei = zeros(ntry,1);
obj_si = zeros(ntry,1);  
for it = 1:1:ntry    
    logb00 = logb0_trend(pickB0:pickBis(it));    
    rate = length(logb00)./length(loga00);
    if rate<0.50 || rate>2.0
        continue;
    end    
    logb1 = sub_allignB2A(loga00,logb00);
    logb1(isnan(logb1)) = 0;        
    
    logb11 = (logb1-mean(logb1))./std(logb00);%%%% useful for short segments
% %     logb11 = logb1;
    
    nrad_corr = 16;
    [obj_simi,obj_area] = sub_Search_objFun_segments(loga11,-logb11,nrad_corr);        
    objvec_simi   = (abs(obj_simi));
    objval_si     = norm(objvec_simi,2);                
    objvec_area   = obj_area; 
    objval_ei     = norm(objvec_area,2);     
    obj_ei(it) = objval_ei;
    obj_si(it) = objval_si;  
    
%     figure(2001)
%     clf;
%     subplot(211)
%     plot(loga00);
%     hold on;
%     plot(logb1);
%     title(it)
%     subplot(212)
%     plot(obj_ei./40);
%     hold on;
%     plot(obj_si./20);    
end
obj_ei1 = (obj_ei-min(obj_ei(:)))./(max(obj_ei(:))-min(obj_ei(:)));
obj_ei1 = smooth(medfilt2(obj_ei1,[3 1]));
obj_si1 = (obj_si-min(obj_si(:)))./(max(obj_si(:))-min(obj_si(:)));
obj_com = (obj_ei1.^2).*obj_si1;

[~,tmpcidx] = max(obj_com,[],1);
pickBi = pickBis(tmpcidx);

logb00 = logb0_trend(pickB0:pickBi);  
if isempty(logb00)==0    
    logb11 = sub_allignB2A(loga00,logb00);
    figure(2001)
    clf;
    subplot(211)
    plot(loga00);
    hold on;
    plot(logb00);
    title('before allignment')
    subplot(212)
    plot(loga00);
    hold on;
    plot(logb11);
    title('after allignment')
end

return;

end
