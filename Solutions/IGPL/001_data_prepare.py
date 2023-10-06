# Library Import 
import glob
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy import interpolate
from tqdm import tqdm

# Hyperparameter or What we need to consider
train_set_idx = [0,1,2,3,4,5,7,8] # training well index
valid_set_idx = [6]               # validation well index # We will test only [0, 4, 8] index for K-fold 

num_train_sample = 5000           # number of train sample from train wells 
num_valid_sample = 1000           # number of train sample from valid well

train_well_log_length = 3000      # patch length
valid_well_log_length = 3000      # patch length

anchor_range = (0.499,0.501)      # anchor points portion
#########################################################

# Data path & Data path list
DATA_files_path = "/h1/bhoon/SPWLA_proj/data/train/*.csv"
train_valid_flist = sorted(glob.glob(DATA_files_path))
# Save path
save_train_path = "/h1/bhoon/SPWLA_proj/new_train/"
save_valid_path = "/h1/bhoon/SPWLA_proj/new_valid/"

# Train dataset(DataFrame) merge
train_df = []
for i in train_set_idx:
    fname = train_valid_flist[i]
    df0 = pd.read_csv(fname)
    df0['log_RD']=np.log10(df0['RD'])
    train_df.append(df0.copy())

# Validation dataset(DataFrame) merge
valid_df = []
for i in valid_set_idx:
    fname = train_valid_flist[i]
    df0 = pd.read_csv(fname)
    df0['log_RD']=np.log10(df0['RD'])
    valid_df.append(df0.copy())

# Random Streching - **Important** - FOR SUPERVISED LEARNING
def random_stretch(curve,depth0,ank_p=0.2): 
    length=curve.shape[0] # Well log length 
    depth=depth0.copy()   # True depth
    n=int(length*ank_p)   # number of anchor point
    anchor_points0=np.sort(np.random.choice(np.arange(0,length),n,replace=False)) # anchor points1
    anchor_points1=np.sort(np.random.choice(np.arange(0,length),n,replace=False)) # anchor points2
    depth1=np.interp(depth,depth[anchor_points0],depth[anchor_points1],left=-99999,right=99999) # misaligned depth
    good=(-99999<depth1)*(depth1<99999) # Inside of data coverage 
    depth1[depth1==99999]=depth1[good][-1]-depth[good][-1]+depth[depth1==99999] # Consider out of anchor point coverage, Continuty
    depth1[depth1==-99999]=depth1[good][0]-depth[good][0]+depth[depth1==-99999] # Consider out of anchor point coverage, Continuty
    f=interpolate.interp1d(depth1,curve,fill_value=(curve[0], curve[-1]), bounds_error=False,kind=1) # Interpolation for getting misaligned well log data
    new_curve = f(depth0) # misaligned data
    if new_curve[-1]==0: # Temporary hardcoding to prevent problems in interpolation
        new_curve[-1]=new_curve[-2] # The problem was that the last value would be zero.
    #new_curve += np.random.normal(np.mean(new_curve),np.std(new_curve)**2,new_curve.shape)*0.03 # Process of adding random noise
    new_depth=np.interp(depth0, depth1, depth0,left=depth0[0],right=depth0[-1]) # correct depth for misaligned well log data 
    return new_curve,new_depth 

def get_0_frame(leng=train_well_log_length): # Function that produces a zero-filled DataFrame
    np_arr = np.zeros((leng,14))
    cols   = ["DEPT","GR","RHOB","NPHI","log_RD",
              "new_RHOB_depth","new_RHOB","final_RHOB",
              "new_NPHI_depth","new_NPHI","final_NPHI",
              "new_log_RD_depth","new_log_RD","final_log_RD"]
    return pd.DataFrame(np_arr,columns=cols)


df = train_df.copy()
# minimun_length = 3000 <<< can consider different patch length
for i in tqdm(range(num_train_sample)):

    # Randomly choosing well index
    random_well_idx = np.random.randint(0,len(df))
    random_well_data = df[random_well_idx]

    # random_length = np.random.randint(minimun_length,len(random_well_data)-1)
    random_length = train_well_log_length
    while True:
        random_start_depth_idx = np.random.randint(0,len(random_well_data)-1)
        if random_start_depth_idx+random_length <= len(random_well_data):
            extracted_well_log_sample = random_well_data.iloc[random_start_depth_idx:random_start_depth_idx+random_length]
            break
    
    # Log sample - patch
    tmp_df = extracted_well_log_sample[["DEPT","GR","RHOB","NPHI","log_RD"]]

    ank_p = np.random.uniform(anchor_range[0],anchor_range[-1])
    new_RHOB,new_RHOB_depth = random_stretch(tmp_df["RHOB"].values,tmp_df["DEPT"].values,ank_p=ank_p)
    tmp_df["new_RHOB_depth"] = new_RHOB_depth
    tmp_df["new_RHOB"] = new_RHOB
    f=interpolate.interp1d(new_RHOB_depth,new_RHOB,fill_value=(new_RHOB[0], new_RHOB[-1]), bounds_error=False,kind=1) 
    tmp_df["final_RHOB"]=f(tmp_df["DEPT"].values) 

    ank_p = np.random.uniform(anchor_range[0],anchor_range[-1])
    new_NPHI,new_NPHI_depth  = random_stretch(tmp_df["NPHI"].values,tmp_df["DEPT"].values,ank_p=ank_p)
    tmp_df["new_NPHI_depth"] = new_NPHI_depth
    tmp_df["new_NPHI"] = new_NPHI
    f=interpolate.interp1d(new_NPHI_depth,new_NPHI,fill_value=(new_NPHI[0], new_NPHI[-1]), bounds_error=False,kind=1) 
    tmp_df["final_NPHI"]=f(tmp_df["DEPT"].values)

    ank_p = np.random.uniform(anchor_range[0],anchor_range[-1])
    new_log_RD,new_log_RD_depth = random_stretch(tmp_df["log_RD"].values,tmp_df["DEPT"].values,ank_p=ank_p)
    tmp_df["new_log_RD_depth"]  = new_log_RD_depth
    tmp_df["new_log_RD"] = new_log_RD
    f=interpolate.interp1d(new_log_RD_depth,new_log_RD,fill_value=(new_log_RD[0], new_log_RD[-1]), bounds_error=False,kind=1) 
    tmp_df["final_log_RD"]=f(tmp_df["DEPT"].values)

    mk_df = get_0_frame(train_well_log_length)
    mk_df[:len(tmp_df)] = tmp_df[:len(tmp_df)].values

    fname = save_train_path+f"{str(i).zfill(5)}.csv"
    mk_df.to_csv(fname,index=False)

# minimun_length = 3000
df = valid_df.copy()
for i in tqdm(range(num_valid_sample)):

    # Randomly choosing well index
    random_well_idx = np.random.randint(0,len(df))
    random_well_data = df[random_well_idx]

    # random_length = np.random.randint(minimun_length,len(random_well_data)-1)
    random_length = valid_well_log_length
    while True:
        random_start_depth_idx = np.random.randint(0,len(random_well_data)-1)
        if random_start_depth_idx+random_length <= len(random_well_data):
            extracted_well_log_sample = random_well_data.iloc[random_start_depth_idx:random_start_depth_idx+random_length]
            break

    # Log sample - patch
    tmp_df = extracted_well_log_sample[["DEPT","GR","RHOB","NPHI","log_RD"]]

    ank_p = np.random.uniform(anchor_range[0],anchor_range[-1])
    new_RHOB,new_RHOB_depth = random_stretch(tmp_df["RHOB"].values,tmp_df["DEPT"].values,ank_p=ank_p)
    tmp_df["new_RHOB_depth"] = new_RHOB_depth
    tmp_df["new_RHOB"] = new_RHOB
    f=interpolate.interp1d(new_RHOB_depth,new_RHOB,fill_value=(new_RHOB[0], new_RHOB[-1]), bounds_error=False,kind=1)
    tmp_df["final_RHOB"]=f(tmp_df["DEPT"].values)

    ank_p = np.random.uniform(anchor_range[0],anchor_range[-1])
    new_NPHI,new_NPHI_depth  = random_stretch(tmp_df["NPHI"].values,tmp_df["DEPT"].values,ank_p=ank_p)
    tmp_df["new_NPHI_depth"] = new_NPHI_depth
    tmp_df["new_NPHI"] = new_NPHI
    f=interpolate.interp1d(new_NPHI_depth,new_NPHI,fill_value=(new_NPHI[0], new_NPHI[-1]), bounds_error=False,kind=1)
    tmp_df["final_NPHI"]=f(tmp_df["DEPT"].values)

    ank_p = np.random.uniform(anchor_range[0],anchor_range[-1])
    new_log_RD,new_log_RD_depth = random_stretch(tmp_df["log_RD"].values,tmp_df["DEPT"].values,ank_p=ank_p)
    tmp_df["new_log_RD_depth"]  = new_log_RD_depth
    tmp_df["new_log_RD"] = new_log_RD
    f=interpolate.interp1d(new_log_RD_depth,new_log_RD,fill_value=(new_log_RD[0], new_log_RD[-1]), bounds_error=False,kind=1)
    tmp_df["final_log_RD"]=f(tmp_df["DEPT"].values)

    mk_df = get_0_frame(valid_well_log_length)
    mk_df[:len(tmp_df)] = tmp_df[:len(tmp_df)].values

    fname = save_valid_path+f"{str(i).zfill(5)}.csv"
    mk_df.to_csv(fname,index=False)
    
    
