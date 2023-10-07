function  csvdata = sub_read_the_Train_Data(file01,Headline_str)
%%%%%% Developer : Xuekai Sun  @ CNLC(China National Logging Corporation)
%%%%%% Date      : Apr-15,2023
%%%%% Copyright (c) 2023, CNLC(China National Logging Corporation)
%%%%% this program is designed for SPWLA PDDA Tests.

% % % Headline_str = 'DEPT,GR,RHOB,NPHI,RD';
ncol = 5;
nline_per_block = 200;

fid = fopen(file01,'rt');
tmp = fgetl(fid);
if strcmpi(tmp,Headline_str )==0
    disp('Wrong Data Format.Please check the csv file!');
    return;    
else
    csvfmt_spec = repmat('%f,',[1 5]);
    csvfmt_spec(end) = [];
end

csvdata = [];
while (~feof(fid))      
    tmp_dat = fscanf(fid,[csvfmt_spec '\n'],ncol.*nline_per_block);
    tmp_dat = reshape(tmp_dat,[ncol length(tmp_dat(:))./ncol]).';    
    csvdata  = [csvdata; tmp_dat];  
end

return;

end
