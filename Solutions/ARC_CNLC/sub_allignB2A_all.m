function  [blog_al,bdep_al] = sub_allignB2A_all(alog,adep,blog,bdep)

    alog = alog(:);
    blog = blog(:);
    adep = adep(:);
    bdep = bdep(:);
    
    na = length(alog);
    nb = length(blog);    
    
    try 
        bdep_al = linspace(bdep(1),bdep(end),na);    
        blog_al = interp1(bdep,blog,bdep_al,'spline');
    catch
        bdep_al = [];
        blog_al = [];        
    end
    
   
%     figure(111)
%     clf;
%     plot(bdep,blog,'b-');
%     hold on;
%     plot(bdep_al,blog_al,'r--');
    
    
return;
end
