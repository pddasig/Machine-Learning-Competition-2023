function [ridx_min,cidx_min] = Get_RotationBack_Coordinates(dat0,angleii)

[ny,nx] = size(dat0);

dat0(1,1,:)   = 1;
dat0(1,nx,:)  = 1;
dat0(ny,1,:)  = 1;
dat0(ny,nx,:) = 1;
dati       = imrotate(dat0, angleii);
dati_bb    = imrotate(dati,-angleii);
[n1_bb,n2_bb] = size(dati_bb);
extr_area = all(dati_bb~=0,3);
col_vec = [1:1:n2_bb];
row_vec = [1:1:n1_bb];
[cc,rr] = meshgrid(col_vec,row_vec);
ridx_min = min(rr(extr_area));
cidx_min = min(cc(extr_area));




return;

end