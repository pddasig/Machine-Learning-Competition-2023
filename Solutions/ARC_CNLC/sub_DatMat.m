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