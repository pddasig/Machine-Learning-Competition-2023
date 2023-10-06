import glob
import pandas as pd
import numpy as np 
from tqdm import tqdm
import matplotlib.pyplot as plt
import tensorflow as tf
from dtwalign import dtw
from scipy import interpolate

# Hyperparameter or What we need to consider
scaler = "standard" # standard, max, minmax
display = "on" # on, off
convolve = False #True #False
conv_window = 30 #3

# custom loss
def custom_loss(y_true, y_pred):
    return tf.math.reduce_mean(tf.square(y_true - y_pred))

########################################################
sample_idx = np.arange(0,100,10) #[10,20,30,40,50,60,70,80,90]
valid_files_path="/h1/bhoon/SPWLA_proj/new_valid/*.csv"
saved_model_path ="/h1/bhoon/SPWLA_proj/best_value_model.h5"

total_data_describe = "/h1/bhoon/SPWLA_proj/total.csv"
total_df_c = pd.read_csv(total_data_describe)
DEPT_mu,    GR_mu,    RHOB_mu,    NPHI_mu,    LOG_RD_mu    = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['mean']
DEPT_sigma, GR_sigma, RHOB_sigma, NPHI_sigma, LOG_RD_sigma = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['std']
DEPT_max,   GR_max,   RHOB_max,   NPHI_max,   LOG_RD_max   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['max']
DEPT_min,   GR_min,   RHOB_min,   NPHI_min,   LOG_RD_min   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['min']

input_columns = ["DEPT","GR",'new_RHOB','new_NPHI','new_log_RD']
output_columns = ["RHOB","NPHI","log_RD"]

data = []
X = []
y = []
for fname in tqdm(sorted(glob.glob(valid_files_path))):
    df0=pd.read_csv(fname)
    scaled_df0=df0.copy()

    # Depth >> minmax or max scaling

    #scaled_df0["DEPT"]-=DEPT_mu
    #scaled_df0["DEPT"]/=DEPT_sigma

    #scaled_df0["DEPT"]/=DEPT_max
    # total_df_c = pd.read_csv(fname)
    # DEPT_mu,    GR_mu,    RHOB_mu,    NPHI_mu,    LOG_RD_mu    = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['mean']
    # DEPT_sigma, GR_sigma, RHOB_sigma, NPHI_sigma, LOG_RD_sigma = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['std']
    # DEPT_max,   GR_max,   RHOB_max,   NPHI_max,   LOG_RD_max   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['max']
    # DEPT_min,   GR_min,   RHOB_min,   NPHI_min,   LOG_RD_min   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['min']

    #scaled_df0["DEPT"]-=DEPT_min
    #scaled_df0["DEPT"]/=(DEPT_max-DEPT_min)

    if convolve:
        scaled_df0["GR"]         =np.convolve(scaled_df0["GR"].values,np.ones(conv_window)/conv_window,mode='same')         
        scaled_df0["RHOB"]       =np.convolve(scaled_df0["RHOB"].values,np.ones(conv_window)/conv_window,mode='same')
        scaled_df0["NPHI"]       =np.convolve(scaled_df0["NPHI"].values,np.ones(conv_window)/conv_window,mode='same')
        scaled_df0["log_RD"]     =np.convolve(scaled_df0["log_RD"].values,np.ones(conv_window)/conv_window,mode='same')
        scaled_df0['new_RHOB']   =np.convolve(scaled_df0['new_RHOB'].values,np.ones(conv_window)/conv_window,mode='same')
        scaled_df0['new_NPHI']   =np.convolve(scaled_df0['new_NPHI'].values,np.ones(conv_window)/conv_window,mode='same')
        scaled_df0['new_log_RD'] =np.convolve(scaled_df0['new_log_RD'].values,np.ones(conv_window)/conv_window,mode='same')

    # standard scaling
    if scaler == "standard":
        #scaled_df0["DEPT"]-=DEPT_mu
        #scaled_df0["DEPT"]/=DEPT_sigma
        scaled_df0["GR"]-=GR_mu
        scaled_df0["GR"]/=GR_sigma
        scaled_df0["RHOB"]-=RHOB_mu
        scaled_df0["RHOB"]/=RHOB_sigma
        scaled_df0["NPHI"]-=NPHI_mu
        scaled_df0["NPHI"]/=NPHI_sigma
        scaled_df0["log_RD"]-=LOG_RD_mu
        scaled_df0["log_RD"]/=LOG_RD_sigma
        scaled_df0['new_RHOB']-=RHOB_mu
        scaled_df0['new_RHOB']/=RHOB_sigma
        scaled_df0['new_NPHI']-=NPHI_mu
        scaled_df0['new_NPHI']/=NPHI_sigma
        scaled_df0['new_log_RD']-=LOG_RD_mu
        scaled_df0['new_log_RD']/=LOG_RD_sigma

    # max scaling
    elif scaler == "max":
        #scaled_df0["DEPT"]/=DEPT_max
        scaled_df0["GR"]/=GR_max
        scaled_df0["RHOB"]/=RHOB_max
        scaled_df0["NPHI"]/=NPHI_max
        scaled_df0["log_RD"]/=LOG_RD_max
        scaled_df0['new_RHOB']/=RHOB_max
        scaled_df0['new_NPHI']/=NPHI_max
        scaled_df0['new_log_RD']/=LOG_RD_max

    # minmax scaling
    elif scaler == "minmax":
        #scaled_df0["DEPT"]-=DEPT_min
        #scaled_df0["DEPT"]/=(DEPT_max-DEPT_min)
        scaled_df0["GR"]-=GR_min
        scaled_df0["GR"]/=(GR_max-GR_min)
        scaled_df0["RHOB"]-=RHOB_min
        scaled_df0["RHOB"]/=(RHOB_max-RHOB_min)
        scaled_df0["NPHI"]-=NPHI_min
        scaled_df0["NPHI"]/=(NPHI_max-NPHI_min)
        scaled_df0["log_RD"]-=LOG_RD_min
        scaled_df0["log_RD"]/=(LOG_RD_max-LOG_RD_min)
        scaled_df0['new_RHOB']-=RHOB_min
        scaled_df0['new_RHOB']/=(RHOB_max-RHOB_min)
        scaled_df0['new_NPHI']-=NPHI_min
        scaled_df0['new_NPHI']/=(NPHI_max-NPHI_min)
        scaled_df0['new_log_RD']-=LOG_RD_min
        scaled_df0['new_log_RD']/=(LOG_RD_max-LOG_RD_min)

    X.append(np.array(scaled_df0[input_columns]))
    y.append(np.array(scaled_df0[output_columns]))

    data.append(df0)

# Model load
model = tf.keras.models.load_model(saved_model_path,custom_objects={"custom_loss":custom_loss})

for idx in sample_idx:
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

    alignment_OBE_rhob   = dtw(mis_rhob,pred_rhob    )
    alignment_OBE_nphi   = dtw(mis_nphi,pred_nphi    )
    alignment_OBE_log_rd = dtw(mis_log_rd,pred_log_rd)
    
    path_rhob   = alignment_OBE_rhob.get_warping_path(target="reference")
    path_nphi   = alignment_OBE_nphi.get_warping_path(target="reference")
    path_log_rd = alignment_OBE_log_rd.get_warping_path(target="reference")
    
    dep = dat["DEPT"][path_rhob]
    mis = dat["new_RHOB"].values
    f=interpolate.interp1d(dep,mis,fill_value=(mis[0], mis[-1]), bounds_error=False,kind=1)
    final_rhob=f(dat["DEPT"])
    final_rhob[0]=final_rhob[1]
    final_rhob[-1]=final_rhob[-2]
    
    dep = dat["DEPT"][path_nphi]
    mis = dat["new_NPHI"].values
    f=interpolate.interp1d(dep,mis,fill_value=(mis[0], mis[-1]), bounds_error=False,kind=1)
    final_nphi=f(dat["DEPT"])
    final_nphi[0]=final_nphi[1]
    final_nphi[-1]=final_nphi[-2]
    
    dep = dat["DEPT"][path_log_rd]
    mis = dat["new_log_RD"].values
    f=interpolate.interp1d(dep,mis,fill_value=(mis[0], mis[-1]), bounds_error=False,kind=1)
    final_log_rd=f(dat["DEPT"])
    final_log_rd[0]=final_log_rd[1]
    final_log_rd[-1]=final_log_rd[-2]

    
    if display=='on':
        # Display
        # Show 1
        fig, ax = plt.subplots(4,1,figsize=(25,15),layout="constrained")
        ax[0].set_title("Gamma ray")
        ax[1].set_title("Bulk density")
        ax[2].set_title("Neutron")
        ax[3].set_title("Log10(Resistivity)")
        ax[0].plot(dat["DEPT"], dat["GR"])
        ax[1].plot(dat["DEPT"], dat["RHOB"],'green')
        ax[2].plot(dat["DEPT"], dat["NPHI"],'r')
        ax[3].plot(dat["DEPT"], dat["log_RD"],'orange')
        ax[1].plot(dat["DEPT"], dat["new_RHOB"],'k',alpha=0.7)
        ax[2].plot(dat["DEPT"], dat["new_NPHI"],'k',alpha=0.7)
        ax[3].plot(dat["DEPT"], dat["new_log_RD"],'k',alpha=0.7)
        plt.savefig(f"img/{idx}_01.png")
        plt.close()
        #plt.show()
        
        # Show 2 
        fig, ax = plt.subplots(4,1,figsize=(25,15),layout="constrained")
        ax[0].set_title("Gamma ray")
        ax[1].set_title("Bulk density")
        ax[2].set_title("Neutron")
        ax[3].set_title("Log10(Resistivity)")
        ax[0].plot(dept, gr)
        ax[1].plot(dept, y[idx][:,0],'green')
        ax[2].plot(dept, y[idx][:,1],'r')
        ax[3].plot(dept, y[idx][:,2],'orange')
        ax[1].plot(dept, pred_rhob,'k',alpha=0.7)
        ax[2].plot(dept, pred_nphi,'k',alpha=0.7)
        ax[3].plot(dept, pred_log_rd,'k',alpha=0.7)
        plt.savefig(f"img/{idx}_02.png")
        plt.close()
        #plt.show()
        
        # Show 3
        fig, ax = plt.subplots(4,1,figsize=(25,15),layout="constrained")
        ax[0].set_title("Gamma ray")
        ax[1].set_title("Bulk density")
        ax[2].set_title("Neutron")
        ax[3].set_title("Log10(Resistivity)")
        ax[0].plot(dat["DEPT"], dat["GR"])
        ax[1].plot(dat["DEPT"], dat["RHOB"],'green')
        ax[2].plot(dat["DEPT"], dat["NPHI"],'r')
        ax[3].plot(dat["DEPT"], dat["log_RD"],'orange')
        ax[1].plot(dat["DEPT"][path_rhob],  dat["new_RHOB"],'k',alpha=0.7)
        ax[2].plot(dat["DEPT"][path_nphi],  dat["new_NPHI"],'k',alpha=0.7)
        ax[3].plot(dat["DEPT"][path_log_rd],dat["new_log_RD"],'k',alpha=0.7)
        plt.savefig(f"img/{idx}_03.png")
        plt.close()
        #plt.show()

        # Show 4
        fig, ax = plt.subplots(4,1,figsize=(25,15),layout="constrained")
        ax[0].set_title("Gamma ray")
        ax[1].set_title("Bulk density")
        ax[2].set_title("Neutron")
        ax[3].set_title("Log10(Resistivity)")
        ax[0].plot(dat["DEPT"], dat["GR"])
        ax[1].plot(dat["DEPT"], dat["RHOB"],'green')
        ax[2].plot(dat["DEPT"], dat["NPHI"],'r')
        ax[3].plot(dat["DEPT"], dat["log_RD"],'orange')
        ax[1].plot(dat["DEPT"], final_rhob,'k',alpha=0.7)
        ax[2].plot(dat["DEPT"], final_nphi,'k',alpha=0.7)
        ax[3].plot(dat["DEPT"], final_log_rd,'k',alpha=0.7)
        plt.savefig(f'img/{idx}_04.png')
        plt.close()
        #plt.show()
 
    print("#"*30)
    print(f"sample index : {idx}")
    # Calculation NMSE
    rhob_err   = np.mean((dat["new_RHOB"].values-dat["RHOB"].values)**2) / (np.std(dat["RHOB"].values)**2)
    nphi_err   = np.mean((dat["new_NPHI"].values-dat["NPHI"].values)**2) / (np.std(dat["NPHI"].values)**2)
    log_rd_err = np.mean((dat["new_log_RD"].values-dat["log_RD"].values)**2) / (np.std(dat["log_RD"].values)**2)
    print("\nValue Error - Error Check ")
    print(f"RHOB : {rhob_err}\nNPHI : {nphi_err}\nlog_RD : {log_rd_err}")
    print(f"Total Error : {(rhob_err+nphi_err+log_rd_err)/3}")
    ##################################################################################################
    rhob_err   = np.mean((final_rhob-dat["RHOB"].values)**2) / (np.std(dat["RHOB"].values)**2)
    nphi_err   = np.mean((final_nphi-dat["NPHI"].values)**2) / (np.std(dat["NPHI"].values)**2)
    log_rd_err = np.mean((final_log_rd-dat["log_RD"].values)**2) / (np.std(dat["log_RD"].values)**2)
    print("\nValue Error - Prediction performence")
    print(f"RHOB : {rhob_err}\nNPHI : {nphi_err}\nlog_RD : {log_rd_err}")
    print(f"Total Error : {(rhob_err+nphi_err+log_rd_err)/3}")

    print("")

    # Calculation MAD
    rhob_dept_err   = np.mean(np.abs(dat["DEPT"].values-dat["new_RHOB_depth"].values))
    nphi_dept_err   = np.mean(np.abs(dat["DEPT"].values-dat["new_NPHI_depth"].values))
    log_rd_dept_err = np.mean(np.abs(dat["DEPT"].values-dat["new_log_RD_depth"].values))
    print("\nDepth Error - Error Check")
    print(f"RHOB : {rhob_dept_err}\nNPHI : {nphi_dept_err}\nlog_RD : {log_rd_dept_err}")
    print(f"Total Error : {(rhob_dept_err+nphi_dept_err+log_rd_dept_err)/3}")
    ##################################################################################################
    rhob_dept_err   = np.mean(np.abs(dat["DEPT"].values-dat["DEPT"].values[path_rhob]))
    nphi_dept_err   = np.mean(np.abs(dat["DEPT"].values-dat["DEPT"].values[path_nphi]))
    log_rd_dept_err = np.mean(np.abs(dat["DEPT"].values-dat["DEPT"].values[path_log_rd]))
    print("\nDepth Error - Prediction performence")
    print(f"RHOB : {rhob_dept_err}\nNPHI : {nphi_dept_err}\nlog_RD : {log_rd_dept_err}")
    print(f"Total Error : {(rhob_dept_err+nphi_dept_err+log_rd_dept_err)/3}")
