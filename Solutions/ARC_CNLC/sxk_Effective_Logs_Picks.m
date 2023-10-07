function  [a1r,b1r] = sxk_Effective_Logs_Picks(a1,b1,pad_numval,nrad_bg,nrad_scope)

idxori        = [1:1:length(a1.trend)].';
% % jud_effvals   = a1.trend~=pad_numval & b1.trend~=pad_numval;
jud_effvals_aa   = a1.trend~=pad_numval;
jud_effvals_bb   = b1.trend~=pad_numval;
idxori_aa        = idxori(jud_effvals_aa);
idxori_bb        = idxori(jud_effvals_bb);

a1r.trend     = a1.trend; %%% the length is not changed 
a1r.detail    = a1.detail;
a1r.depth     = a1.depth;
a1r.idxori    = idxori_aa;

b1r.trend     = b1.trend;%%% the length is not changed 
b1r.detail    = b1.detail;
b1r.depth     = b1.depth;
b1r.idxori    = idxori_bb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % nrad_bg   = 32;
% % % nrad_scope   = 3;
[~,~,AL1_d,~] = sub_hierarchy_Peaks(a1r.detail,nrad_scope,nrad_bg,4,2,1);
[~,~,BL1_d,~] = sub_hierarchy_Peaks(b1r.detail,nrad_scope,nrad_bg,4,2,1);
[~,~,AL1_t,~] = sub_hierarchy_Peaks(diff(a1r.detail+a1r.trend),nrad_scope,nrad_bg,4,2,1);
[~,~,BL1_t,~] = sub_hierarchy_Peaks(diff(b1r.detail+b1r.trend),nrad_scope,nrad_bg,4,2,1);

%%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%
AL0 = unique([AL1_d(:);AL1_t(:)]);
BL0 = unique([BL1_d(:);BL1_t(:)]);

% [aa,bb]    = meshgrid(AL0,BL0);
% ddmat      = abs(aa-bb);
% jud_b2a    = any(ddmat<8,2);
% jud_a2b    = any(ddmat<8,1);
% bidx_compl = BL0(~jud_b2a);
% aidx_compl = AL0(~jud_a2b);
% AL1        = unique([AL0(:);bidx_compl(:)]);
% BL1        = unique([BL0(:);aidx_compl(:)]);

AL1        = unique([min(idxori_aa); AL0(:); BL0(:); max(idxori_aa)]);
BL1        = unique([min(idxori_bb); BL0(:); AL0(:); max(idxori_bb)]);


AL1(AL1<min(idxori_aa)) = [];
AL1(AL1>max(idxori_aa)) = [];
BL1(BL1<min(idxori_bb)) = [];
BL1(BL1>max(idxori_bb)) = [];

a1r.pick1 = AL1;
b1r.pick1 = BL1;

return;

end
