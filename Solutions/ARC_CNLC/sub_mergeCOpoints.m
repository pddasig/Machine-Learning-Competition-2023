function  [output,record_cell]= sub_mergeCOpoints(loc_rawCO,dif_rawC0,min_nrad,max_nrad,maxiter)
%%%%%% Developer : Xuekai Sun  @ CNLC(China National Logging Corporation)
%%%%%% Date      : Apr-21,2023
%%%%% Copyright (c) 2023, CNLC(China National Logging Corporation)
%%%%% this program is designed for SPWLA PDDA Tests.

locs_prev = loc_rawCO;
difs_prev = dif_rawC0;
nrad_prev = min_nrad;

niter       = 0;
record_cell = [];
while niter<=maxiter    
    nrad_now = 2.^niter.*nrad_prev;
    if nrad_now>=max_nrad
        break;
    end
        
    refine_locs = [];
    refine_difs = [];
    while isempty(locs_prev)==0    
        locii = locs_prev(1);
        loc01 = locii-nrad_now;
        loc02 = locii+nrad_now;

        jud_merge = locs_prev>=loc01 & locs_prev<=loc02;
        loc_ten   = locs_prev(jud_merge);
        dif_ten   = difs_prev(jud_merge);
        [~,idx] = min(abs(dif_ten));    
        refine_locs = [refine_locs;loc_ten(idx)];
        refine_difs = [refine_difs;dif_ten(idx)];
        
        locs_prev(jud_merge) = [];
        difs_prev(jud_merge) = [];
    end
    
    niter = niter + 1;
    record_cell{1,niter} = [refine_locs(:).'; refine_difs(:).'];
    
    locs_prev = refine_locs;
    difs_prev = refine_difs;
    
% % %     figure(111)
% % %     clf;
% % %     plot(loc_rawCO,dif_rawC0,'ko');
% % %     hold on;
% % %     plot(refine_locs,refine_difs,'r.');
    
end

output = record_cell{1,end}(1,:);


return;

end
