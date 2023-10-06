import pandas as pd
import glob
import numpy as np

train_flist = sorted(glob.glob("/h1/bhoon/SPWLA_proj/data/train/*.csv"))
test_flist = sorted(glob.glob("/h1/bhoon/SPWLA_proj/data/test/*.csv"))

total_df = []

train_df = []
for i in train_flist:
    df0=pd.read_csv(i)
    df0['well']=i.split('_')[-1].split('.')[0] 
    df0['log_RD']=np.log10(df0['RD']) 
    train_df.append(df0.copy()) 
    total_df.append(df0.copy()) 
train_df_c=pd.concat(train_df)

test_df = []
for i in test_flist:
    df0=pd.read_csv(i)
    df0['well']="1"+i.split('_')[-1].split('.')[0]
    df0['log_RD']=np.log10(df0['RD']) 
    test_df.append(df0.copy()) 
    total_df.append(df0.copy()) 
test_df_c=pd.concat(test_df)

total_df_c=pd.concat(total_df)

save_path = "/h1/bhoon/SPWLA_proj/total.csv"
total_df_c.to_csv(save_path,index=False)
