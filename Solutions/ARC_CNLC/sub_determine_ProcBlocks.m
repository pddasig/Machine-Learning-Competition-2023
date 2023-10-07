function  [proc_idx01,proc_idx02,nblock] = sub_determine_ProcBlocks(nshot_tar,Proc_RTlen,Proc_RTpad)

% Proc_RTlen   = 256; 
% Proc_RTpad   = 64; 
if Proc_RTlen > nshot_tar
    Proc_RTlen = nshot_tar;
    Proc_RTpad = ceil(nshot_tar./8); 
end       
proc_idx01 = [];
proc_idx02 = [];
curr_idx01 = 1;
flag       = 1;
while flag==1
    curr_idx02 = curr_idx01+Proc_RTlen-1;
    curr_idx02(curr_idx02>=nshot_tar) = nshot_tar;
    curr_idx01 = curr_idx02-Proc_RTlen+1;    
    proc_idx01 = [proc_idx01 curr_idx01];
    proc_idx02 = [proc_idx02 curr_idx02];    
    if curr_idx02==nshot_tar        
           break;    
    end
    curr_idx01 =   curr_idx02+1-Proc_RTpad;       
end
nblock = length(proc_idx01);


return;

end
