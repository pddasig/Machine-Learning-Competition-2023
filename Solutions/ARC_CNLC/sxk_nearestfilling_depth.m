function  [rhob_fill] = sxk_nearestfilling_depth(rhob_pred,pad_numval)

jud_ineff = rhob_pred==pad_numval;

idxvec = [1:1:length(rhob_pred)];
idxvec_good = idxvec(~jud_ineff);
idxvec_inef = idxvec(jud_ineff);
rhob_good = rhob_pred(idxvec_good);


[bb,gg] = meshgrid(idxvec_inef,idxvec_good);
[~,idx] = min(abs(bb-gg),[],1);

rhob_fill = rhob_pred;
rhob_fill(idxvec_inef) = rhob_good(idx);


return;

end