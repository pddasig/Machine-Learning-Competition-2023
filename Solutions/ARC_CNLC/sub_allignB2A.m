function  [logb1] = sub_allignB2A(loga0,logb0)

loga0 = loga0(:);
logb0 = logb0(:);

aidxvec = [1:1:length(loga0)];
bidxvec = [1:1:length(logb0)];
ddb = (aidxvec(end)-aidxvec(1))./(bidxvec(end)-bidxvec(1));

bidxvec_al = bidxvec(1) + (bidxvec-1).*ddb;
logb1 = interp1(bidxvec_al,logb0,aidxvec,'line');
logb1 = logb1(:);

% figure(111)
% clf;
% subplot(211)
% plot(aidxvec,loga0,'b-');
% hold on
% plot(bidxvec_al,logb0,'r--');
% 
% subplot(212)
% plot(bidxvec_al,logb0,'b-');
% hold on
% plot(aidxvec,logb1,'r--');


return;

end