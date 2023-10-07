function [r1,t1] = DTW_newobj_v1_sametype(r,t,rpidx,tpidx,scope_nrad,pad_numval)
% [Dist,D,k,w,rw,tw]=dtw(r,t,pflag)
% Dynamic Time Warping Algorithm
% Dist is unnormalized distance between t and r
% D is the accumulated distance matrix
% k is the normalizing factor
% w is the optimal path
% t is the vector you are testing against
% r is the vector you are testing
% rw is the warped r vector
% tw is the warped t vector
% pflag  plot flag: 1 (yes), 0(no)
%
% Version comments:
% rw, tw and pflag added by Pau Mic
 
% [row,M]=size(r); if (row > M) M=row; r=r'; end
% [row,N]=size(t); if (row > N) N=row; t=t'; end

M = length(rpidx);
N = length(tpidx);

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

% % % t1 = 1 - abs(t - medfilt2(t,[32 1]))./max(t);
% % % r1 = 1 - abs(r - medfilt2(r,[32 1]))./max(r);
% % % [at,ar] = meshgrid(t1(tpidx),r1(rpidx));
% % % aa = (abs(at.*ar));

d = (obj_corr+obj_area);

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
% ar = (max(r) - r(rpidx))./max(r);
% at = (max(t) - t(tpidx))./max(t);

ar = (max(r) - r(rpidx));
at = (max(t) - t(tpidx));

[ridx_tar,cidx_tar] = sxk_getAnchorPoint_row(d,w1,10,ar,at);
r1 = rpidx(ridx_tar);
t1 = tpidx(cidx_tar);


% % % head_pair = [rpidx(w1(1,1)) tpidx(w1(1,2))];
% % % tail_pair = [rpidx(w1(end,1)) tpidx(w1(end,2))];



figure(4002)
subplot(1,10,1:2)
cla;
imagesc(d);
hold on;
plot(w1(:,2),w1(:,1),'r.');
plot(cidx_tar,ridx_tar,'K*');
try
    xlim([1 size(d,2)])
    ylim([1 size(d,1)])
end

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
d1 = d1.*repmat(at(:).',[M 1]); %%%
d1(~jud_cc)  = max(d1(:));
cidx_log1     = sub_rmsAMP(min(d1,[],1)-min(d1(:)),3);
[~,cidx_tar1] = min(cidx_log1(:),[],1);
ridx_log1  = d(:,cidx_tar1);
[d_min1,ridx_tar1] = min(ridx_log1(:),[],1); 

% % %     ridx_tar = ridx_tar1;
% % %     cidx_tar = cidx_tar1;


% % % %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% 
% % % %%%% respect peaks of the reference log
d2 = d;
d2 = d2.*repmat(ar(:),[1 N]);
d2(~jud_cc)  = max(d2(:));

ridx_log2     = sub_rmsAMP(min(d2,[],2)-min(d2(:)),3);
[~,ridx_tar2] = min(ridx_log2(:),[],1);
cidx_log2     = d(ridx_tar2,:);
[d_min2,cidx_tar2] = min(cidx_log2(:),[],1); 

if d_min1<=d_min2
    ridx_tar = ridx_tar1;
    cidx_tar = cidx_tar1;
else
    ridx_tar = ridx_tar2;
    cidx_tar = cidx_tar2;
end

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
