function  [ridx_tar,cidx_tar] = sxk_getAnchorPoint_row(d,w,nrad,ar,at)


[M,N] = size(d);

% % % % nrad = 3;
ridx01    = w(:,1) - nrad; 
ridx01(ridx01<1) = 1;
ridx02    = w(:,1) + nrad; 
ridx02(ridx02>M) = M;
idx1d_cen = w(:,1) + (w(:,2)-1).*M;
idx1d_lw = ridx01 + (w(:,2)-1).*M;
idx1d_up = ridx02 + (w(:,2)-1).*M;
tmp = d;
tmp(abs(tmp)<1e-6)= 1e-4;
tmp(isnan(tmp))   = 1e-4;
tmp_up = tmp;
tmp_up(idx1d_up) = 0;
jud_up = cumprod(tmp_up,1) ==0;
tmp_dw = tmp;
tmp_dw(idx1d_lw) = 0;
jud_dw = cumprod(tmp_dw,1) ==0;
jud_cc = (jud_up + jud_dw)==1;

d(~jud_cc) = max(d(:));

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% 
%%%% respect peaks of the target log
d1 = d;
d1 = d1.*repmat(at(:).',[M 1]);
d1(~jud_cc)  = max(d1(:));
cidx_log1     = sub_rmsAMP(min(d1,[],1),3);
[~,cidx_tar1] = min(cidx_log1(:),[],1);
ridx_log1  = d(:,cidx_tar1);
[d_min1,ridx_tar1] = min(ridx_log1(:),[],1); 

    ridx_tar = ridx_tar1;
    cidx_tar = cidx_tar1;


% % % %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% 
% % % %%%% respect peaks of the reference log
% d2 = d;
% d2 = d2.*repmat(ar(:),[1 N]);
% d2(~jud_cc)  = max(d2(:));
% 
% ridx_log2     = sub_rmsAMP(min(d2,[],2),3);
% [~,ridx_tar2] = min(ridx_log2(:),[],1);
% cidx_log2     = d(ridx_tar2,:);
% [d_min2,cidx_tar2] = min(cidx_log2(:),[],1); 
% 
% if d_min1<=d_min2
%     ridx_tar = ridx_tar1;
%     cidx_tar = cidx_tar1;
% else
%     ridx_tar = ridx_tar2;
%     cidx_tar = cidx_tar2;
% end

return;

end


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