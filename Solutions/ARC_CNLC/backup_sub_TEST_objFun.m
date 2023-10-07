function [obj_simi,obj_area,loc_refCO,loc_rawCO] = sub_TEST_objFun(loga0,logb0,nrad_corr)

% % % nrad_corr = 8;

loga = loga0;
logb = logb0;
ndep_log = length(loga);

idxvec = [1:1:ndep_log];
idxmat = repmat(idxvec,[3 1]) + repmat([-1:1:1].',[1 ndep_log]);
judeff = idxmat>=1 & idxmat<=ndep_log;
vaamat = zeros(3,ndep_log);
vaamat(judeff) = loga(idxmat(judeff));
vbbmat = zeros(3,ndep_log);
vbbmat(judeff) = logb(idxmat(judeff));
vaamin = min(vaamat,[],1);
vaamax = max(vaamat,[],1);
vbbmin = min(vbbmat,[],1);
vbbmax = max(vbbmat,[],1);
vccmin = max([vaamin;vbbmin]);
vccmax = min([vaamax;vbbmax]);

is_rawCO  = (vccmin)<=(vccmax);% is_rawCO(1) = 0;% is_rawCO(end) = 0;
loc_rawCO = idxvec(is_rawCO);
dif_rawC0 = abs(loga(is_rawCO)-logb(is_rawCO));
[loc_refCO,records]= sub_mergeCOpoints(loc_rawCO,dif_rawC0,nrad_corr,nrad_corr*8,3);

% % figure(2002)
% % clf;
% % plot(loga,'r-');
% % hold on;
% % plot(logb,'b-');
% % plot(loc_rawCO,loga(is_rawCO),'go');
% % plot(loc_rawCO,logb(is_rawCO),'go');
% % hold on;
% % plot(loc_refCO,loga(loc_refCO),'r*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_CO = length(loc_refCO)-1;
obj_simi = zeros(1,ndep_log);
obj_area = zeros(1,ndep_log);

% % % thrd_bigCOlen = 256;
% % % num_bigCO     = sum(diff(loc_refCO,1,2)>=thrd_bigCOlen);
% % % num_avgCO     = round(ndep_log./thrd_bigCOlen);

% % % nrad_corr = 8;
for in = 1:1:num_CO    
    idx01 = loc_refCO(in);
    idx02 = loc_refCO(in+1);
    
    len_ii = idx02-idx01+1;
    loga_ii = loga(idx01:idx02);
    logb_ii = logb(idx01:idx02);      
    idxinner = [1:1:len_ii];
    idxmat   = repmat(idxinner,[nrad_corr.*2+1 1]) + repmat([-nrad_corr:1:nrad_corr].',[1 len_ii]);
    judmat   = idxmat>=1 & idxmat<=len_ii;
    aamat = zeros(nrad_corr.*2+1,len_ii);
    aamat(judmat) = loga_ii(idxmat(judmat))+0.1; %%%%%% !!!!!!
    bbmat = zeros(nrad_corr.*2+1,len_ii);
    bbmat(judmat) = logb_ii(idxmat(judmat))+0.1; %%%%%% !!!!!!    
    
    area_pos = median(abs(aamat-bbmat));
    area_neg = median(abs(aamat)+abs(bbmat));
    area_vec = min([area_pos;area_neg],[],1) ;   %%% area_vec = max([area_pos;area_neg],[],1);
    area_avg = 0.5.*(min(abs(aamat),[],1) + min(abs(bbmat),[],1));
    area_vec = (area_vec- area_avg);
        
%     aarms    = sqrt(mean(aamat.^2,1));
%     bbrms    = sqrt(mean(bbmat.^2,1));
%     area_vec = (aarms+bbrms);
%     area_vec(area_vec<0) = 0;
    
    simi_vec  = mean(aamat.*bbmat)./(mean(aamat.^2+bbmat.^2));    
    obj_simi(idx01:idx02) = (simi_vec);
    obj_area(idx01:idx02) = area_vec;
    
% % %     figure(2003)
% % %     clf;   
% % %     plot(loga_ii,'r-');
% % %     hold on;
% % %     plot(logb_ii,'b-');
% % %     plot((simi_vec),'k-','linewidth',1.5);    
% % %     plot(area_pos,'g-','linewidth',1.5);  
% % %     plot(area_neg,'c-','linewidth',1.5);
% % %     plot(area_vec,'k--','linewidth',1.5);    
end
% obj_comb = abs(obj_simi).^0.5.*obj_area;
% norm(obj_comb)


% figure(2002)
% clf;
% subplot(211)
% plot(loga,'r-');
% hold on;
% plot(logb,'b-');
% plot(loc_rawCO,loga(is_rawCO),'go');
% plot(loc_refCO,loga(loc_refCO),'r*');
% subplot(212)
% plot(obj_simi);
% hold on;
% plot(obj_area);
% plot(loc_rawCO,obj_area(is_rawCO),'go');
% plot(loc_refCO,obj_area(loc_refCO),'r*');

return;

end
