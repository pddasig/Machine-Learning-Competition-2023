function  [rw_idx,tw_idx,inner_index] = sub_chopHeadTail(rw_idx0,tw_idx0)

npair0  = length(rw_idx0);
later_idx = [2:1:npair0];
front_idx = [1:1:npair0-1];
ridx01 = min(later_idx(diff(rw_idx0)>0));
ridx02 = max(front_idx(diff(rw_idx0)>0));
inner_index = [ridx01:ridx02];
rw_idx = rw_idx0(inner_index);
tw_idx = tw_idx0(inner_index);


return;

end
