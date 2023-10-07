function [r0,t0,rtidx,d,w1] = sxk_DTW_newobj_Shift_iter(r,t,rpidx,tpidx,scope_nrad,max_shift,pad_numval,picks_now)

M = length(rpidx);
N = length(tpidx);

[tt,rr] = meshgrid(tpidx,rpidx);
mask    = abs(rr-tt)>= (1.0*max_shift);
bkval   = 4;
%%%%%% %%%%%% %%%%%% %%%%%% %%%%%% find the half 
nrad_ref = min([round(length(r)./4) scope_nrad]);
nrad_tar = min([round(length(t)./4) scope_nrad]);
nrad01 = min([nrad_ref nrad_tar]);

[aamat]    = sub_DatMat(r,nrad01);
[bbmat]    = sub_DatMat(t,nrad01);
aamat      = aamat(:,rpidx);
bbmat      = bbmat(:,tpidx);
aamat1 = (aamat - mean(aamat,1))./std(aamat,[],1);
bbmat1 = (bbmat - mean(bbmat,1))./std(bbmat,[],1);
obj_corr       = 1-abs(corr(aamat1,bbmat1));

jud_amat_padval = any(aamat==pad_numval,1);
jud_bmat_padval = any(bbmat==pad_numval,1);
obj_corr(jud_amat_padval,:) = 1;
obj_corr(:,jud_bmat_padval) = 1;


obj_area  = zeros(M,N);
for im = 1:1:M   
    avec = aamat1(:,im);  
    avec0 = aamat(:,im); 
    for in = 1:1:N
        bvec0 = bbmat(:,in); 
        bvec  = bbmat1(:,in);        
        area_all = abs(avec)+abs(bvec);
        area_avg = mean(area_all);   
        area_pos = (abs(avec-bvec));
        area_neg = (abs(avec+bvec));        
        area_vec = max([area_pos area_neg],[],2); 
        area_vec = area_vec - area_avg;            
        area_val = 1-sum(area_vec(area_vec>0))./sum(area_all);   
        
        obj_area(im,in) = area_val;
        obj_dist(im,in) = sqrt(mean(avec0-bvec0).^2);
    end
end

% % figure(111)
% % clf;
% % subplot(221)
% % imagesc(obj_corr); colorbar;
% % subplot(222)
% % imagesc(obj_area); colorbar;
% % subplot(223)
% % imagesc(obj_dist.^0.5);colorbar;
% % subplot(224)
% % imagesc(obj_corr+obj_dist+obj_area);colorbar;


d = (obj_corr+obj_area);
% % % d = (obj_dist).^0.5.*obj_corr+obj_area;

D=zeros(size(d));
D(1,1)=d(1,1); 
for m=2:M
    D(m,1)=d(m,1)+D(m-1,1);
end
for n=2:N
    D(1,n)=d(1,n)+D(1,n-1);
end
for m=2:M
    for n=2:N
        D(m,n)=d(m,n)+min(D(m-1,n),min(D(m-1,n-1),D(m,n-1))); % this double MIn construction improves in 10-fold the Speed-up. Thanks Sven Mensing
    end
end
 
Dist=D(M,N);
n=N;
m=M;
k=1;
w=[M N];
while ((n+m)~=2)
    if (n-1)==0
        m=m-1;
    elseif (m-1)==0
        n=n-1;
    else 
      [values,number]=min([D(m-1,n),1.*D(m,n-1),D(m-1,n-1)]); %%% 2.0 is a penality factor
% % %       [values,number]=min([D(m-1,n),inf,D(m-1,n-1)]);
      switch number
      case 1
        m=m-1;
      case 2
        n=n-1;
      case 3
        m=m-1;
        n=n-1;
      end
  end
    k=k+1;
    w=[m n; w]; % this replace the above sentence. Thanks Pau Mic
end


[~,~,inner_index] = sub_chopHeadTail(rpidx(w(:,1)),tpidx(w(:,2)));
w1 = w(inner_index,:);



d = d - median(d,1);
d = d - median(d,2);
d(mask==1) = bkval;

% % % % for ip = 1:1:size(picks_now,1)
% % % %     ridxii = picks_now(ip,1);
% % % %     tidxii = picks_now(ip,2);
% % % %     jud_ridx = rpidx>=(ridxii-scope_nrad) & rpidx<=(ridxii+scope_nrad);
% % % %     jud_tidx = tpidx>=(tidxii-scope_nrad) & tpidx<=(tidxii+scope_nrad); 
% % % %     
% % % %     tmp1 = d(jud_ridx,:);
% % % %     jud1 = tmp1 < bkval;
% % % %     tmp1(jud1) = tmp1(jud1)./2+bkval/2;
% % % %     d(jud_ridx,:) = tmp1;
% % % %     
% % % %     tmp2 = d(:,jud_tidx);
% % % %     jud2 = tmp2 < bkval;
% % % %     tmp2(jud2) = tmp2(jud2)./2+bkval/2;
% % % %     d(:,jud_tidx) = tmp2;
% % % % end


nrad     = 5;
dr       = imrotate(d,45);
dr(dr==0) = bkval;
% % % % dr = medfilt2(dr,[3 1]);
[ridx_min,cidx_min] = Get_RotationBack_Coordinates(d,45);
dlog1    = medfilt2(min(dr,[],2),[3 1]);
dlog2    = medfilt2(sum(dr,2)./(sum(dr~=0,2)+1),[nrad*2+1 1]);
dlog     = dlog1 + dlog2;


[~,ridx] = min(dlog,[],1);
ridx01   = max([1 ridx-nrad]);
ridx02   = min([length(dlog) ridx+nrad]);
dr1      = zeros(size(dr));
dr1(ridx01:ridx02,:) = dr(ridx01:ridx02,:);
dr1      = imrotate(dr1,-45);
dr_back  = dr1(ridx_min:ridx_min+M-1,cidx_min:cidx_min+N-1);
dr_back(dr_back==0) = bkval;


dr_log = min(dr_back,[],2);
dr_log = medfilt2(dr_log,[3 1]);
[~,r00] = min(dr_log);
dc_log  = min(dr_back,[],1);
dc_log = medfilt2(dc_log,[1 3]);
[~,t00] = min(dc_log);

r0 = rpidx(r00);
t0 = tpidx(t00);

rtidx = [r00 t00];




% % % % L    = min([M,N]);
% % % % nrad     = 3;
% % % % search_mat = [];
% % % % 
% % % % %%%% downward processing 
% % % % dw_limit = max([(M-2*nrad)-min(w1(:,1)), 0]);
% % % % diag_dw = [0:nrad:dw_limit];
% % % % for id = 1:1:length(diag_dw)    
% % % %     ridx = [1:1:L].' + diag_dw(id);
% % % %     cidx = [1:1:L].';
% % % %     jud_ridx = ridx>=1 & ridx<=M;    
% % % %     wi     = [ridx(jud_ridx) cidx(jud_ridx)];    
% % % %     [rtar,ctar,dtar] = sxk_getAnchorPoint_row_v0(d,wi,nrad);      
% % % %     rct_curr = [rtar,ctar,dtar];
% % % %     search_mat = [search_mat;rct_curr];       
% % % % 
% % % % %     figure(221)
% % % % %     clf;
% % % % %     imagesc(d);
% % % % %     hold on;
% % % % %     plot(cidx,ridx,'b.');
% % % % %     plot(ctar,rtar,'r*');
% % % % %     title(num2str(dtar*1e4));
% % % %     
% % % % end
% % % % 
% % % % %%%% upward processing 
% % % % up_limit = max([(M-2*nrad)-min(w1(:,1)), 0]);
% % % % diag_up  = [0:nrad:up_limit];
% % % % for iu = 1:1:length(diag_up)    
% % % %     ridx = [1:1:L].' ;
% % % %     cidx = [1:1:L].' + diag_up(iu);
% % % %     jud_cidx = cidx>=1 & cidx<=N;    
% % % %     wi     = [ridx(jud_cidx),cidx(jud_cidx)];
% % % %     
% % % %     [rtar,ctar,dtar] = sxk_getAnchorPoint_row_v0(d,wi,nrad);     
% % % %     rct_curr = [rtar,ctar,dtar];
% % % %     search_mat = [search_mat;rct_curr];    
% % % % %     figure(221)
% % % % %     clf;
% % % % %     imagesc(d);
% % % % %     hold on;
% % % % %     plot(cidx,ridx,'b.');
% % % % %     plot(ctar,rtar,'r*');
% % % % %     title(num2str(dtar*1e4));
% % % %     
% % % % end
% % % % [~,idx] = min(search_mat(:,3));
% % % % r00     = search_mat(idx,1);
% % % % t00     = search_mat(idx,2);
% % % figure(222)
% % % clf;
% % % imagesc(d);
% % % hold on;
% % % plot(t00,r00,'r*');
% % % 
% % % r0 = rpidx(r00);
% % % t0 = tpidx(t00);


return;

end

function [datmat] = sub_DatMat(suba0,nrad)

% nrad    = 16;
suba0   = suba0(:);
idxveca = [1:1:length(suba0)];
idxmat  = repmat(idxveca,[2.*nrad+1 1]) + repmat([-nrad:1:nrad].',[1 length(suba0)]);
judmat  = idxmat>=1 & idxmat<=length(suba0);
datmat = zeros(2.*nrad+1,length(suba0));
datmat(judmat) = suba0(idxmat(judmat));

return;

end

function  [ridx_tar,cidx_tar,d_min] = sxk_getAnchorPoint_row_v0(d,w,nrad)


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
d1(~jud_cc)  = max(d1(:));
cidx_log0    = medfilt2(min(d1,[],1),[1 3]);
cidx_log1    = sub_rmsAMP(cidx_log0-2.*min(cidx_log0),3);
cidx_log2    = sub_rmsAMP(cidx_log0-2.*min(cidx_log0),20);
[~,cidx_tar1] = min(cidx_log1(:) + cidx_log2(:),[],1);
ridx_log1  = d1(:,cidx_tar1);
[dmin1,ridx_tar1] = min(ridx_log1(:),[],1); 

% % dmin2 = sqrt(mean(cidx_log1));
d_min  = dmin1;

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
