function  pairs = sxk_effDTWpairs(b2a_route,flat_min_space)

np      = length(b2a_route);

fidxmat = [1:1:np-1].';
judflat = diff(b2a_route,1,1)==0;
jud_ar_pidx = (judflat(:,1)==1 & judflat(:,2)==0);
jud_br_pidx = (judflat(:,2)==1 & judflat(:,1)==0);
sub_ar_pidx   = fidxmat(jud_ar_pidx);
sub_br_pidx   = fidxmat(jud_br_pidx);

ar_pidx = b2a_route(sub_ar_pidx,1);
br_pidx = b2a_route(sub_br_pidx,2);
jud_ar_edge = [diff(ar_pidx(:),1,1)>1; 0==1];
jud_br_edge = [diff(br_pidx(:),1,1)>1; 0==1];

sub_ar_pidx = sub_ar_pidx(jud_ar_edge);
sub_br_pidx = sub_br_pidx(jud_br_edge);

pidx = unique([sub_br_pidx(:);sub_ar_pidx(:)]);
pairs = b2a_route(pidx,:);

jud_eff = all([diff(pairs,1,1)>=flat_min_space; [1 1]==1],2);
pairs = pairs(jud_eff,:);


% % % ar_pidx1 = b2a_route(sub_ar_pidx,1);
% % % br_pidx1 = b2a_route(sub_br_pidx,2);

return;

end
