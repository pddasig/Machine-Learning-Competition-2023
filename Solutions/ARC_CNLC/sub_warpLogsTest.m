function  [loga0,logb0,idxaa_back] = sub_warpLogsTest(loga00,logb00,tenfac)

idxaa = [1:1:length(loga00)];
idxbb = [1:1:length(logb00)];

ten_idxbb = (idxbb-idxbb(1)).*tenfac + idxbb(1);

tmp_depmin = max([min(idxaa) min(ten_idxbb)]);
tmp_depmax = min([max(idxaa) max(ten_idxbb)]);
jud_depthcoverage = idxaa>=tmp_depmin & idxaa<=tmp_depmax;
dep_coverage      = idxaa(jud_depthcoverage); %%% covered 
    
logb0   = interp1(ten_idxbb,logb00,dep_coverage,'linear');
loga0   = loga00(jud_depthcoverage);  

logb0 = logb0(:);
loga0 = loga0(:);

idxaa_back = idxaa(jud_depthcoverage);


% figure(111)
% clf;
% subplot(211)
% plot(idxbb,logb00,'b-');
% hold on;
% plot(dep_coverage,logb0,'r--');
% 
% subplot(212)
% plot(dep_coverage,loga0,'b-');
% hold on;
% plot(dep_coverage,logb0,'r-');


return;

end
