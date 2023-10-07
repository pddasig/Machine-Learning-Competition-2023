function  sxk_DispCheck_PairPairs(a1r,b1r,PPairs_s1,pad_numval)


idxvec_a1r = [1:1:length(a1r.trend)];
judeff_a1r = a1r.trend~=pad_numval & a1r.detail~=pad_numval;
idxvec_b1r = [1:1:length(b1r.trend)];
judeff_b1r = b1r.trend~=pad_numval & b1r.detail~=pad_numval;
idxvec_a1r_act = idxvec_a1r(judeff_a1r);
idxvec_b1r_act = idxvec_b1r(judeff_b1r);

% switch mode
%     case 'detail'
%         alog_act = a1r.detail(judeff_a1r);
%         blog_act = b1r.detail(judeff_b1r);  
%     case 'trend'
%         alog_act = a1r.trend(judeff_a1r);
%         blog_act = b1r.trend(judeff_b1r); 
%     case 'all'
%         alog_act = a1r.detail(judeff_a1r) + a1r.trend(judeff_a1r);
%         blog_act = b1r.detail(judeff_b1r) + b1r.trend(judeff_b1r);  
% end

figure(55555)
clf;
subplot(211)
alog_act = a1r.trend(judeff_a1r);
blog_act = b1r.trend(judeff_b1r); 

plot(idxvec_a1r_act,alog_act,'-','linewidth',1.5);
hold on;
plot(idxvec_b1r_act,blog_act,'-');
plot(PPairs_s1(:,1),a1r.trend(PPairs_s1(:,1)),'r*');
plot(PPairs_s1(:,1),a1r.trend(PPairs_s1(:,1)),'ro');
plot(PPairs_s1(:,2),b1r.trend(PPairs_s1(:,2)),'b*');
plot(PPairs_s1(:,2),b1r.trend(PPairs_s1(:,2)),'bo');
xlabel('Depth Index');
ylabel('');
grid on;
grid minor;
box on;
title('Trend Log');
subplot(212)
alog_act = a1r.detail(judeff_a1r) + a1r.trend(judeff_a1r);
blog_act = b1r.detail(judeff_b1r) + b1r.trend(judeff_b1r); 
plot(idxvec_a1r_act,alog_act,'-','linewidth',1.5);
hold on;
plot(idxvec_b1r_act,blog_act,'-');
plot(PPairs_s1(:,1),a1r.trend(PPairs_s1(:,1)) + a1r.detail(PPairs_s1(:,1)),'r*');
plot(PPairs_s1(:,1),a1r.trend(PPairs_s1(:,1)) + a1r.detail(PPairs_s1(:,1)),'ro');
plot(PPairs_s1(:,2),b1r.trend(PPairs_s1(:,2)) + b1r.detail(PPairs_s1(:,2)),'b*');
plot(PPairs_s1(:,2),b1r.trend(PPairs_s1(:,2)) + b1r.detail(PPairs_s1(:,2)),'bo');
xlabel('Depth Index');
ylabel('');
grid on;
grid minor;
box on;
title('Original Log');


return;

end
