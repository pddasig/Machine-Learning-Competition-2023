import glob
import pandas as pd
import numpy as np 
from tqdm import tqdm
from tensorflow import keras
from tensorflow.keras.layers import *
from keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras.optimizers import *
import tensorflow as tf

# Hyperparameter or What we need to consider
scaler = "standard" # standard, max, minmax
optimizer = Adam(learning_rate=1e-4) # 0.0001  
metrics = "mse"
epoch = 1000
batch_size = 32    
patience = 5       # For Early Stop
convolve = False   # True #False
conv_window = 30

# custom loss
def custom_loss(y_true, y_pred):
    return tf.math.reduce_mean(tf.square(y_true - y_pred))
    #return tf.math.reduce_mean(tf.abs(y_true - y_pred))

total_data_describe = "/h1/bhoon/SPWLA_proj/total.csv"
train_files_path="/h1/bhoon/SPWLA_proj/new_train/*.csv"
valid_files_path="/h1/bhoon/SPWLA_proj/new_valid/*.csv"

save_model_path ="/h1/bhoon/SPWLA_proj/best_value_model.h5"

total_df_c = pd.read_csv(total_data_describe)
DEPT_mu,    GR_mu,    RHOB_mu,    NPHI_mu,    LOG_RD_mu    = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['mean']
DEPT_sigma, GR_sigma, RHOB_sigma, NPHI_sigma, LOG_RD_sigma = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['std']
DEPT_max,   GR_max,   RHOB_max,   NPHI_max,   LOG_RD_max   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['max']
DEPT_min,   GR_min,   RHOB_min,   NPHI_min,   LOG_RD_min   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['min']

#input_columns = ["DEPT","GR",'new_RHOB','new_NPHI','new_log_RD']
input_columns = ["GR",'new_RHOB','new_NPHI','new_log_RD']
output_columns = ["RHOB","NPHI","log_RD"] #,"new_RHOB_depth","new_NPHI_depth","new_log_RD_depth"]

train_flist = sorted(glob.glob(train_files_path))
valid_flist = sorted(glob.glob(valid_files_path))

train_X = []
train_y = []
for fname in tqdm(train_flist):

    df0=pd.read_csv(fname)
    scaled_df0=df0.copy() 

    # total_df_c = pd.read_csv(fname)
    # DEPT_mu,    GR_mu,    RHOB_mu,    NPHI_mu,    LOG_RD_mu    = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['mean']
    # DEPT_sigma, GR_sigma, RHOB_sigma, NPHI_sigma, LOG_RD_sigma = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['std']
    # DEPT_max,   GR_max,   RHOB_max,   NPHI_max,   LOG_RD_max   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['max']
    # DEPT_min,   GR_min,   RHOB_min,   NPHI_min,   LOG_RD_min   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['min']

    # Depth >> minmax or max scaling

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

    train_X.append(np.array(scaled_df0[input_columns]))
    train_y.append(np.array(scaled_df0[output_columns]))
    
train_X=np.array(train_X)
train_y=np.array(train_y)

valid_X = []
valid_y = []
for fname in tqdm(valid_flist):
    df0=pd.read_csv(fname)
    scaled_df0=df0.copy() 
    
    # total_df_c = pd.read_csv(fname)
    
    # DEPT_mu,    GR_mu,    RHOB_mu,    NPHI_mu,    LOG_RD_mu    = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['mean']
    # DEPT_sigma, GR_sigma, RHOB_sigma, NPHI_sigma, LOG_RD_sigma = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['std']
    # DEPT_max,   GR_max,   RHOB_max,   NPHI_max,   LOG_RD_max   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['max']
    # DEPT_min,   GR_min,   RHOB_min,   NPHI_min,   LOG_RD_min   = total_df_c[["DEPT","GR","RHOB","NPHI","log_RD"]].describe().loc['min']

    # Depth >> minmax or max scaling

    #scaled_df0["DEPT"]-=DEPT_mu
    #scaled_df0["DEPT"]/=DEPT_sigma

    #scaled_df0["DEPT"]/=DEPT_max

    #scaled_df0["DEPT"]-=DEPT_min
    #scaled_df0["DEPT"]/=(DEPT_max-DEPT_min)
    #scaled_df0["new_RHOB_depth"]-=DEPT_min
    #scaled_df0["new_RHOB_depth"]/=(DEPT_max-DEPT_min)
    #scaled_df0["new_NPHI_depth"]-=DEPT_min
    #scaled_df0["new_NPHI_depth"]/=(DEPT_max-DEPT_min)
    #scaled_df0["new_log_RD_depth"]-=DEPT_min
    #scaled_df0["new_log_RD_depth"]/=(DEPT_max-DEPT_min)

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

    valid_X.append(np.array(scaled_df0[input_columns]))
    valid_y.append(np.array(scaled_df0[output_columns]))
    
valid_X=np.array(valid_X)
valid_y=np.array(valid_y)

# Build Model Architecture
inputs = keras.Input(shape=(None,1,4))

# Model - Functional
x = inputs
c1 = Conv2D(64, (200,1),padding='same')(x)
c1 = BatchNormalization()(c1)
c1 = ReLU()(c1)
c1 = Dropout(0.2)(c1)

t1 = concatenate([c1,x])

c2 = Conv2D(64,(200,1), padding='same')(t1)
c2 = BatchNormalization()(c2)
c2 = ReLU()(c2)
c2 = Dropout(0.2)(c2)

t2 = concatenate([c2,x])

c3 = Conv2D(64,(200,1), padding='same')(t2)
c3 = BatchNormalization()(c3)
c3 = ReLU()(c3)
c3 = Dropout(0.2)(c3)

t3 = concatenate([c3,x])

c4 = Conv2D(64, (200,1),padding='same')(t3)
c4 = BatchNormalization()(c4)
c4 = ReLU()(c4)
c4 = Dropout(0.2)(c4)

t4 = concatenate([c4,x])

c5 = Conv2D(64,(200,1), padding='same')(t4)
c5 = BatchNormalization()(c5)
c5 = ReLU()(c5)
c5 = Dropout(0.2)(c5)

t5 = concatenate([c5,x])
outputs = Conv2D(3,(200,1), activation='linear', padding='same')(t5)

model = keras.Model(inputs=inputs, outputs=outputs)
model.summary()

# Early Stopping & Model Checkpoint
earlystopper = EarlyStopping(patience=patience, restore_best_weights=True)
model_checkpoint = ModelCheckpoint(save_model_path,monitor='val_loss',save_best_only=True)

# Model compile 
model.compile(optimizer=optimizer,loss=custom_loss,metrics=metrics)

train_X = np.expand_dims(train_X,-2)
train_y = np.expand_dims(train_y,-2)
valid_X = np.expand_dims(valid_X,-2)
valid_y =  np.expand_dims(valid_y,-2)

print(f"train_X.shape : {train_X.shape}")
print(f"train_y.shape : {train_y.shape}")
print(f"valid_X.shape : {valid_X.shape}")
print(f"valid_y.shape : {valid_y.shape}")

# Model Train
result = model.fit(train_X, train_y,
                   epochs=epoch,
                   batch_size=batch_size,
                   validation_data=(valid_X, valid_y),
                   callbacks=[earlystopper, model_checkpoint])
