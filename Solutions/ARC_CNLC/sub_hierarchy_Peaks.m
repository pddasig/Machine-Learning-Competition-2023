function  [pidx_L3,pidx_L2,pidx_L1,pidx0] = sub_hierarchy_Peaks(loga0_detail,nrad_scope,nrad_bg,L3,L2,L1)

% % % nrad_scope   = 3;
% % % nrad_bg = 16;

login   = abs(loga0_detail);
ndepth  = length(login);
idxvec0 = [1:1:ndepth];
idxmat  = repmat(idxvec0,[nrad_scope*2+1 1]) + repmat([-nrad_scope:1:nrad_scope].',[1 ndepth]);
judmat  = idxmat>=1 & idxmat<=ndepth;
bgmat   = zeros(nrad_scope*2+1,ndepth);
bgmat(judmat) = login(idxmat(judmat));
ispeaks       = all(login(:).'>=bgmat,1) & login(:).'> abs(median(loga0_detail));

pidx00   = idxvec0(ispeaks);
pamp00   = login(ispeaks);
pidx0     = [];
peak_amps = [];
while isempty(pidx00)==0    
    idxii      = pidx00(1);
    judnearby  = abs(pidx00-idxii)<=(nrad_scope+1);    
    idxii_can  = pidx00(judnearby);    
    ampii_can  = pamp00(judnearby);
    [~,idx]    = max(ampii_can);
    
    pidx0     = [pidx0;     idxii_can(idx)];
    peak_amps = [peak_amps; ampii_can(idx)];
    pidx00(judnearby) = [];  
    pamp00(judnearby) = [];  
end
peak_amps_ref = medfilt2(peak_amps,[nrad_bg*2+1 1]);

isL1 = peak_amps > L1.*peak_amps_ref;
isL2 = peak_amps > L2.*peak_amps_ref;
isL3 = peak_amps > L3.*peak_amps_ref & (peak_amps > L3.*median(peak_amps_ref));

% figure(111)
% clf;
% plot(login);
% hold on;
% plot(pidx0(isL1),peak_amps(isL1),'r*');


pidx_L1 = pidx0(isL1);
% pidx_L1(pidx_L1<=(pidx_L1(1)-nrad_bg)) = [];
% pidx_L1(pidx_L1>=(pidx_L1(end)-nrad_bg)) = [];

pidx_L2 = pidx0(isL2);
% pidx_L2(pidx_L2<=(pidx_L2(1)-nrad_bg)) = [];
% pidx_L2(pidx_L2>=(pidx_L2(end)-nrad_bg)) = [];

pidx_L3 = pidx0(isL3);
% pidx_L3(pidx_L3<=(pidx_L3(1)-nrad_bg)) = [];
% pidx_L3(pidx_L3>=(pidx_L3(end)-nrad_bg)) = [];


return;


end
