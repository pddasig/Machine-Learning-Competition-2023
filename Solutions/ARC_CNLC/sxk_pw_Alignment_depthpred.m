function bdepth_pred = sxk_pw_Alignment_depthpred(pickpairs,dep_suba0,dep_subb0,pad_numval)

pickpairs = sortrows(pickpairs,1);
npair = size(pickpairs,1);

ndepb0 = length(dep_subb0);
bdepth_pred = zeros(ndepb0,1); %%% find where the b sample locates after alignments
for ip = 1:1:npair-1   
    a1   = [pickpairs(ip,1):1:pickpairs(ip+1,1)];    %%% reference index
    b1   = [pickpairs(ip,2):1:pickpairs(ip+1,2)];    %%% target index 
    
    adepii = dep_suba0(a1);    
    bdepii = dep_subb0(b1); 
    bdep_predii = linspace(adepii(1),adepii(end),length(bdepii));        
    bdepth_pred(b1) = bdep_predii;
end

%%%% %%%% %%%% if any head to be processed
jud_padval  = dep_subb0 == pad_numval;
head_shift  = dep_suba0(pickpairs(1,1)) - dep_subb0(pickpairs(1,2));
bdepth_pred(1:pickpairs(1,2)) = dep_subb0(1:pickpairs(1,2)) + head_shift;

tail_shift  = dep_suba0(pickpairs(end,1)) - dep_subb0(pickpairs(end,2));
bdepth_pred(pickpairs(end,2):end) = dep_subb0(pickpairs(end,2):end) + tail_shift;

bdepth_pred(bdepth_pred==jud_padval) = pad_numval;

return;

end
