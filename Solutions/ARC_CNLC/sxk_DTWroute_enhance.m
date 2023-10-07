function  [route_out] = sxk_DTWroute_enhance(d,route,Apicks,mflen_hor,wrad_ver,s1_scope)

[M,N] = size(d);

bkval = 4;
dr       = imrotate(d,45);
dr(dr==0) = bkval;
[ridx_min,cidx_min] = Get_RotationBack_Coordinates(d,45);
dr1 = medfilt2(dr,[1 mflen_hor]);
dr2 = sub_meanAMP(dr,wrad_ver);

drr = (dr1+dr2)./2;

dr3     = imrotate(drr,-45);
dr_back  = dr3(ridx_min:ridx_min+M-1,cidx_min:cidx_min+N-1);
dr_back(dr_back==0) = bkval;

%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%
idx1d       = (route(:,2)-1).*M + route(:,1);
route_value = dr_back(idx1d);
route_mat = [route route_value(:)];
Apick_route = Apicks(route_mat(:,1));

%%%%% iteration begins now
value_eff  = dr_back(dr_back<bkval);
value_base = median(value_eff(:));
jud_route  = route_mat(:,3)<=value_base;


route_mat = route_mat(jud_route,:);
Apick_route = Apick_route(jud_route,:);
route_out = [];
while isempty(route_mat)==0    
    [~,idx] = min(route_mat(:,3));
    ppi = route_mat(idx,:);
    aai = Apick_route(idx);
    route_out = [route_out; ppi];    
    
    Apick_min =  aai - s1_scope*1;
    Apick_max =  aai + s1_scope*1;    
    jud_clear =  Apick_route>=Apick_min & Apick_route<=Apick_max;
    
    route_mat(jud_clear,:) = [];
    Apick_route(jud_clear) = []; 
end




return;

end


function rms = sub_meanAMP(datin,nrad)

[nt,nshot] = size(datin);

tidx_vec = [1:1:nt];
tidx_mat = repmat(tidx_vec,[2.*nrad+1 1]) + repmat([-nrad:1:nrad].',[1 nt]);
tjud_mat = tidx_mat>=1 & tidx_mat<=nt;
dat_mat = zeros(2.*nrad+1,nt);

rms = zeros(nt,nshot);
for is = 1:1:nshot
    datii = datin(:,is);
    dat_mat(tjud_mat) = datii(tidx_mat(tjud_mat));
    rms(:,is) = mean(dat_mat);
end




return;

end
