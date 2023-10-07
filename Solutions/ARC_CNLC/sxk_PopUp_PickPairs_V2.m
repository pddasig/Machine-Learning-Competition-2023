function  pickpairs = sxk_PopUp_PickPairs_V2(a1,b1,sub_Apicks,sub_Bpicks,s0_scope,s1_scope,pad_numval)


nrad_bg      = 32;
nrad_scope   = 3;
[a1r,b1r]    = sxk_Effective_Logs_Picks(a1,b1,pad_numval,nrad_bg,nrad_scope);



idxvec   = [1:1:length(suba0)].';
judeff   = suba0~=pad_numval & subb0~=pad_numval;
idx_head = min(idxvec(judeff));
idx_tail = max(idxvec(judeff));

% % % % % [a0,b0]  = DTW_newobj_v1_sametype(suba0,subb0,sub_Apicks,sub_Bpicks,s0_scope,s1_scope,pad_numval);
[a0,b0] = sxk_DTW_newobj_Shift_iter(suba0,subb0,sub_Apicks,sub_Bpicks,s1_scope,round(s1_scope./2),pad_numval);

pickpairs = [a0,b0];
jud_clearA = ~(sub_Apicks >= (a0-s1_scope) & sub_Apicks <= (a0+s1_scope));
jud_clearB = ~(sub_Bpicks >= (b0-s1_scope) & sub_Bpicks <= (b0+s1_scope));
sub_Apicks = sub_Apicks(jud_clearA);
sub_Bpicks = sub_Bpicks(jud_clearB);

a0_idx   = [1:1:length(suba0)];
b0_idx   = [1:1:length(subb0)];
judeffa0 = suba0~=pad_numval;
judeffb0 = subb0~=pad_numval;

figure(4002)
clf;
set(gcf,'position',[-2235 699 2498 420]);
subplot(1,10,3:10)
hold on
plot(a0_idx(judeffa0),suba0(judeffa0));
plot(b0_idx(judeffb0),subb0(judeffb0));
plot(a0,suba0(a0),'r^');
plot(b0,subb0(b0),'b^');
grid on;
grid minor;
box on;
xlabel('depth index')

ht_tmpl = [idx_head idx_head; idx_tail idx_tail];
ntimes = max([round(length(suba0)./(1.*s1_scope)),20]);
for it = 1:1:ntimes      
    if length(sub_Apicks)<5 || length(sub_Bpicks)<5
        break;
    end     
% %     [ai,bi]  = DTW_newobj_v1_sametype(suba0,subb0,sub_Apicks,sub_Bpicks,s1_scope,pad_numval);
    [ai,bi,rtidx,di,ri] = sxk_DTW_newobj_Shift_iter(suba0,subb0,sub_Apicks,sub_Bpicks,s1_scope,round(s1_scope./2),pad_numval);
    new_pair = [ai bi]; 

    jud_eff  = all(sign(pickpairs(:,1)-ai) == sign(pickpairs(:,2)-bi)) & suba0(ai)~=pad_numval & subb0(bi)~=pad_numval;    
    if ~jud_eff 
        jud_clearA = ~(sub_Apicks >= (ai-1.*s1_scope) & sub_Apicks <= (ai+1.*s1_scope));
        jud_clearB = ~(sub_Bpicks >= (bi-1.*s1_scope) & sub_Bpicks <= (bi+1.*s1_scope));
        sub_Apicks = sub_Apicks(jud_clearA);
        sub_Bpicks = sub_Bpicks(jud_clearB);
        continue;
    end 

    jud_goodenough = abs(ai-bi)<1.*s1_scope;
    if jud_goodenough==1
        tmp_pair = sortrows([a0 b0; ai bi; ht_tmpl]);
        pickpairs = [pickpairs;new_pair]; 

        figure(4002)
        subplot(1,10,1:2)
        cla;
        imagesc(di);
        hold on;
        plot(ri(:,2),ri(:,1),'r.');
        plot(rtidx(2),rtidx(1),'K*');
        xlim([1 size(di,2)])
        ylim([1 size(di,1)])
        subplot(1,10,3:10)
        plot(ai,suba0(ai),'r*');
        plot(ai,suba0(ai),'ro');
        hold on
        plot(bi,subb0(bi),'b*');
        plot(bi,subb0(bi),'bo');    


        

% % %         [subb1, ~,suba1, ~] = sxk_pw_Alignment_FreeStyle(tmp_pair,suba0,subb0,idxvec,idxvec,pad_numval);
%         ai = [];        bi = [];
%         suba0 = suba1;
%         subb0 = subb1;       
    end    

        jud_clearA = ~(sub_Apicks >= (ai-1.*s1_scope) & sub_Apicks <= (ai+1.*s1_scope));
        jud_clearB = ~(sub_Bpicks >= (bi-1.*s1_scope) & sub_Bpicks <= (bi+1.*s1_scope));
        sub_Apicks = sub_Apicks(jud_clearA);
        sub_Bpicks = sub_Bpicks(jud_clearB);
    

end
pickpairs = sortrows(pickpairs,1);

return;

end
