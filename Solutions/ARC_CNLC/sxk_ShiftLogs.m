function  [loga2,logb2,depa2,depb2] = sxk_ShiftLogs(loga0,logb0,depa0,depb0,Pa0,Pb0,pad_numval)


num_shift = Pa0-Pb0; %%%% 

pad_zeros = ones(abs(num_shift),1).*pad_numval;

switch sign(num_shift)
    case -1        
        loga1 = [pad_zeros;loga0(:);];      
        depa1 = [pad_zeros;depa0(:);];
        logb1 = logb0;
        depb1 = depb0;        
    case +1        
        logb1 = [pad_zeros;logb0(:);];      
        depb1 = [pad_zeros;depb0(:);];
        loga1 = loga0;
        depa1 = depa0;          
    case 0
        logb1 = logb0;
        depb1 = depb0;  
        loga1 = loga0;
        depa1 = depa0;    
end

numa1 = length(loga1);
numb1 = length(logb1);
numab = max([numa1 numb1]);

tmp_log = ones(numab,1)*pad_numval;
loga2 = tmp_log;
loga2(1:numa1) = loga1;
logb2 = tmp_log;
logb2(1:numb1) = logb1;

tmp_dep = ones(numab,1)*pad_numval;
depa2 = tmp_dep;
depa2(1:numa1) = depa1;
depb2 = tmp_dep;
depb2(1:numb1) = depb1;


% figure(101)
% clf;
% subplot(211)
% plot(loga0);
% hold on;
% plot(logb0);
% subplot(212)
% plot(loga1);
% hold on;
% plot(logb1);



return;

end
