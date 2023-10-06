import glob
import pandas as pd
import numpy as np 
from tqdm import tqdm
import matplotlib.pyplot as plt
import tensorflow as tf
from dtwalign import dtw
from scipy import interpolate

#submission_path  = "/h1/bhoon/SPWLA_proj/submission/"
submission_path  = "./"
test_files_path  = "/h1/bhoon/SPWLA_proj/data/test/*.csv"
saved_model_path = "/h1/bhoon/SPWLA_proj/best_model.h5"
scaler = "standard"
total_data_describe = "/h1/bhoon/SPWLA_proj/total.csv"
total_df_c = pd.read_csv(total_data_describe)
DEPT_mu,    GR_mu,    RHOB_mu,    NPHI_mu,    LOG_RD_mu    = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['mean']
DEPT_sigma, GR_sigma, RHOB_sigma, NPHI_sigma, LOG_RD_sigma = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['std']
DEPT_max,   GR_max,   RHOB_max,   NPHI_max,   LOG_RD_max   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['max']
DEPT_min,   GR_min,   RHOB_min,   NPHI_min,   LOG_RD_min   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['min']

# custom loss
def custom_loss(y_true, y_pred):
    return tf.math.reduce_mean(tf.square(y_true - y_pred))


data = []
X = []
for fname in tqdm(sorted(glob.glob(test_files_path))):
    df0=pd.read_csv(fname)
    df0['log_RD']=np.log10(df0['RD'])
    scaled_df0=df0.copy()

    #total_df_c = pd.read_csv(fname)
    #DEPT_mu,    GR_mu,    RHOB_mu,    NPHI_mu,    LOG_RD_mu    = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['mean']
    #DEPT_sigma, GR_sigma, RHOB_sigma, NPHI_sigma, LOG_RD_sigma = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['std']
    #DEPT_max,   GR_max,   RHOB_max,   NPHI_max,   LOG_RD_max   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['max']
    #DEPT_min,   GR_min,   RHOB_min,   NPHI_min,   LOG_RD_min   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['min']

    # Depth >> minmax or max scaling

    #scaled_df0["DEPT"]-=DEPT_mu
    #scaled_df0["DEPT"]/=DEPT_sigma

    #scaled_df0["DEPT"]/=DEPT_max

    scaled_df0["DEPT"]-=DEPT_min
    scaled_df0["DEPT"]/=(DEPT_max-DEPT_min)

    # standard scaling
    if scaler == "standard":
#        scaled_df0["DEPT"]-=DEPT_mu
#        scaled_df0["DEPT"]/=DEPT_sigma
        scaled_df0["GR"]-=GR_mu
        scaled_df0["GR"]/=GR_sigma
        scaled_df0["RHOB"]-=RHOB_mu
        scaled_df0["RHOB"]/=RHOB_sigma
        scaled_df0["NPHI"]-=NPHI_mu
        scaled_df0["NPHI"]/=NPHI_sigma
        scaled_df0["log_RD"]-=LOG_RD_mu
        scaled_df0["log_RD"]/=LOG_RD_sigma

    # max scaling
    elif scaler == "max":
#        scaled_df0["DEPT"]/=DEPT_max
        scaled_df0["GR"]/=GR_max
        scaled_df0["RHOB"]/=RHOB_max
        scaled_df0["NPHI"]/=NPHI_max
        scaled_df0["log_RD"]/=LOG_RD_max

    # minmax scaling
    elif scaler == "minmax":
#        scaled_df0["DEPT"]-=DEPT_min
#        scaled_df0["DEPT"]/=(DEPT_max-DEPT_min)
        scaled_df0["GR"]-=GR_min
        scaled_df0["GR"]/=(GR_max-GR_min)
        scaled_df0["RHOB"]-=RHOB_min
        scaled_df0["RHOB"]/=(RHOB_max-RHOB_min)
        scaled_df0["NPHI"]-=NPHI_min
        scaled_df0["NPHI"]/=(NPHI_max-NPHI_min)
        scaled_df0["log_RD"]-=LOG_RD_min
        scaled_df0["log_RD"]/=(LOG_RD_max-LOG_RD_min)

    X.append(np.array(scaled_df0[["DEPT","GR","RHOB","NPHI","log_RD"]]))
    data.append(df0)

# Model load
#model = tf.keras.models.load_model(saved_model_path)
model = tf.keras.models.load_model(saved_model_path,custom_objects={"custom_loss":custom_loss})

for idx in range(3):

    Xdat = X[idx]
    dat  = data[idx]

    dept       = Xdat[:,0]
    gr         = Xdat[:,1]
    mis_rhob   = Xdat[:,2]
    mis_nphi   = Xdat[:,3]
    mis_log_rd = Xdat[:,4]
    
    pred_input = np.expand_dims(Xdat[:,1:5],0)
    pred_input = np.expand_dims(pred_input,-2)
    pred = np.squeeze(model.predict(pred_input))

    pred_rhob  = pred[:,0]
    pred_nphi  = pred[:,1]
    pred_log_rd= pred[:,2]
    
    alignment_OBE_rhob   = dtw(mis_rhob,pred_rhob)
    alignment_OBE_nphi   = dtw(mis_nphi,pred_nphi)
    alignment_OBE_log_rd = dtw(mis_log_rd,pred_log_rd)
    
    path_rhob   = alignment_OBE_rhob.get_warping_path(target="reference")
    path_nphi   = alignment_OBE_nphi.get_warping_path(target="reference")
    path_log_rd = alignment_OBE_log_rd.get_warping_path(target="reference")
    
    dep = dat["DEPT"][path_rhob]
    mis = dat["RHOB"].values
    f=interpolate.interp1d(dep,mis,fill_value=(mis[0], mis[-1]), bounds_error=False,kind=1)
    final_rhob=f(dat["DEPT"])
    
    dep = dat["DEPT"][path_nphi]
    mis = dat["NPHI"].values
    f=interpolate.interp1d(dep,mis,fill_value=(mis[0], mis[-1]), bounds_error=False,kind=1)
    final_nphi=f(dat["DEPT"])
    
    dep = dat["DEPT"][path_log_rd]
    mis = dat["log_RD"].values
    f=interpolate.interp1d(dep,mis,fill_value=(mis[0], mis[-1]), bounds_error=False,kind=1)
    final_log_rd=f(dat["DEPT"])

    dat["RHOB_pred"]      = final_rhob
    dat["NPHI_pred"]      = final_nphi
    dat["RD_pred"]        = 10**final_log_rd
    
    dat["RHOB_dept_pred"] = dat["DEPT"].values[path_rhob]
    dat["NPHI_dept_pred"] = dat["DEPT"].values[path_nphi]
    dat["RD_dept_pred"]   = dat["DEPT"].values[path_log_rd]

    dat[dat.columns[:-1]].to_csv(submission_path+f"aligned_well_{str(idx+1).zfill(2)}.csv",index=False)

