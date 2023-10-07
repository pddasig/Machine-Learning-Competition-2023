function  hidx = sxk_get_hidx(gr_ac,pad_numval)


idxvec   = [1:1:length(gr_ac.trend)];
jud_grac = gr_ac.trend~=pad_numval;
idxeff   = idxvec(jud_grac);

hidx = min(idxeff)-1;


return;

end
