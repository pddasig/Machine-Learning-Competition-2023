function rms = sub_rmsAMP(datin,nrad)

[nt,nshot] = size(datin);

tidx_vec = [1:1:nt];
tidx_mat = repmat(tidx_vec,[2.*nrad+1 1]) + repmat([-nrad:1:nrad].',[1 nt]);
tjud_mat = tidx_mat>=1 & tidx_mat<=nt;
dat_mat = zeros(2.*nrad+1,nt);

rms = zeros(nt,nshot);
for is = 1:1:nshot
    datii = datin(:,is);
    dat_mat(tjud_mat) = datii(tidx_mat(tjud_mat));
    rms(:,is) = sqrt(mean(dat_mat.^2));
end




return;

end