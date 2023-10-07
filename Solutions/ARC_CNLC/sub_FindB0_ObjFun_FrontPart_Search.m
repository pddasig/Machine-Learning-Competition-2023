function  obj2d_com = sub_FindB0_ObjFun_FrontPart_Search(loga00,logb00,pickA0,pickB0s,tensile_pool)

% % ndepth_aa = length(loga00);
% % ndepth_bb = length(logb00);
ntry    = length(pickB0s);
nten = length(tensile_pool);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loga0     = loga00(1:pickA0); 
ndepth_a1 = length(loga0);

obj2d_ei = zeros(nten,ntry);
obj2d_si = zeros(nten,ntry);
for itop = 1:1:ntry     
    logb0 = zeros(ndepth_a1,1);  
    B2i   = pickB0s(itop);
    B1i   = max([B2i - ndepth_a1+1,1]);
    logb0(1:B2i-B1i+1) = logb00(B1i:B2i);
    
    logb1 = flipud(logb0(:)); %%% flipud is used for the convenience of objective function calculation
    loga1 = flipud(loga0(:)); 
    judeff = logb1~=0 & loga1~=0;
    logb1  = logb1(judeff);
    loga1  = loga1(judeff);    
    
    if sum(judeff)<=512  %%% designed for searching mode
        obj2d_com(:,itop) = zeros(nten,1);
        continue;
    end
            
    tmp_ei = zeros(nten,1);
    tmp_si = zeros(nten,1);    
   parfor itf = 1:1:nten                   
        [loga11,logb11] = sub_warpLogsTest(loga1,logb1,tensile_pool(itf));   
        logb11(isnan(logb11)) = 0;       
               
        nrad_corr = 16;
        [obj_simi,obj_area] = sub_Search_objFun(loga11,-logb11,nrad_corr); 
        objvec_simi   = abs(obj_simi);
        objval_si     = sqrt(mean(objvec_simi.^2));           
        objvec_area   = obj_area; 
        objval_ei     = norm(objvec_area)./(norm(loga11,2) + norm(logb11,2)); %%% 1-norm(objvec_area,2)
       
        tmp_ei(itf) = objval_ei;
        tmp_si(itf) = objval_si;        
    end                
    obj2d_ei(:,itop) = tmp_ei;
    obj2d_si(:,itop) = tmp_si;       
end
% obj2d_com = obj2d_si;
obj2d_com = obj2d_ei.*2.*obj2d_si;

% figure(1032)
% clf;
% subplot(131)
% imagesc(obj2d_ei);
% subplot(132)
% imagesc(obj2d_si);
% subplot(133)
% imagesc(obj2d_com);


return;

end
