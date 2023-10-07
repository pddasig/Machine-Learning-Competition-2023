function  [bc,ac] = sxk_combo_Shift_Stretch(b0,a0,PPairs_s0,PPairs_s1,pad_numval)

AS0 = PPairs_s0(1);
BS0 = PPairs_s0(2);

[a1.trend, b1.trend, a1.depth,b1.depth] = sxk_ShiftLogs(a0.trend,b0.trend,a0.depth,b0.depth,AS0,BS0,pad_numval);
[a1.detail,b1.detail,a1.depth,b1.depth] = sxk_ShiftLogs(a0.detail,b0.detail,a0.depth,b0.depth,AS0,BS0,pad_numval);

%%% the below is an example
[bc.trend, bc.depori,ac.trend, ac.depth] = sxk_pw_Alignment_FreeStyle(PPairs_s1, a1.trend,  b1.trend, a1.depth,  b1.depth, pad_numval);
[bc.detail,bc.depori,ac.detail,ac.depth] = sxk_pw_Alignment_FreeStyle(PPairs_s1, a1.detail, b1.detail, a1.depth, b1.depth, pad_numval);

%%%%%%%%% bc's depth after alignment (or b2a depth)
dd_alog = mean(diff(a0.depth(a0.depth~=pad_numval)));  %%% this depth log are original, only be shifted
dd_blog = mean(diff(b0.depth(b0.depth~=pad_numval)));  %%% this depth log are original, only be shifted
idxvec   = [1:1:length(ac.trend)];
jud_alog = ac.trend~=pad_numval;
jud_blog = bc.trend~=pad_numval;
aidx_vec = idxvec(jud_alog);
bidx_vec = idxvec(jud_blog);
aidx01 = min(aidx_vec);
aidx02 = max(aidx_vec);
bidx01 = min(bidx_vec);
bidx02 = max(bidx_vec);

b2a_hidx = [bidx01:1:aidx01];
b2a_head_depth = [-length(b2a_hidx)+1:1:0].*dd_blog + ac.depth(aidx01); %%% use dd_blog because of the freestyle head

b2a_tidx       = [aidx02:1:bidx02];
b2a_tail_depth = [0:1:length(b2a_tidx)-1].*dd_blog + ac.depth(aidx02); %%% use dd_blog because of the freestyle head

b2a_depth = ac.depth;
b2a_depth(b2a_hidx) = b2a_head_depth(:);
b2a_depth(b2a_tidx) = b2a_tail_depth(:);

%%%%%% remove ineffective depth for blog
rm_hidx = [aidx01:1:bidx01-1];
rm_tidx = [bidx02+1:1:aidx02];
b2a_depth(rm_hidx) = pad_numval;
b2a_depth(rm_tidx) = pad_numval;
bc.depth = b2a_depth;

return;

end
