function  [d] = sub_calculate_d(r,t,nrad)

M = length(r);
N = length(t);
% % nrad    = 16;

[aamat] = sub_DatMat(r,nrad);
[bbmat] = sub_DatMat(t,nrad);
d0       = 1-abs(corr(aamat,bbmat));
d1 = zeros(M,N);
for im = 1:1:M
    avec = aamat(:,im);
    for in = 1:1:N
        bvec = bbmat(:,in);
        d1(im,in) = 1-abs(sum(avec.*bvec)./(sum(avec.^2+bvec.^2)+0.001));       
    end
end
d = d1;
d(d>d0) = d0(d>d0);

return;

end
