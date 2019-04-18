from numpy import *
from matplotlib.pylab import *
import random
import sys
import io
import os
import glob
import h5py


os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"]="0"

from tensorflow.python.client import device_lib
print(device_lib.list_local_devices())


from keras.callbacks import ModelCheckpoint
from keras.models import Model, load_model, Sequential
from keras.layers import Dense, Activation, Dropout, Input, Masking, TimeDistributed, LSTM, Conv1D, Concatenate
from keras.layers import GRU, Bidirectional, BatchNormalization, Reshape
from keras.optimizers import Adam
from keras.regularizers import L1L2
from keras import backend as K
from keras.layers import MaxPooling2D, AveragePooling2D,Conv2D
from keras.layers import Lambda
from keras.backend import transpose
from keras import regularizers


def get_random_time_segment(segment_frames,total_frames=12000):
    '''
    Gets a random time segment of duration segment_frames in a file
    with number of frames: total_frames
    '''
    
    segment_start = randint(0, high = total_frames-
                                   segment_frames)
    segment_end = segment_start + segment_frames
    
    return (segment_start, segment_end)
    
    
def getdataextract(xt,nframes=1000):
    duration, roisize = xt.shape
    start,end = get_random_time_segment(nframes,total_frames=duration)
    xtt = xt[start:end,:].transpose()
    s1 = xtt.mean()
    s0 = xtt.std()/s1
    s2 = xtt.std()
    s3 = mean((xtt-s1)**3)/s2**3
    s4 = mean((xtt-s1)**4)/s2**4
    xtt0 = xtt.mean(axis=0)
    xtt = xtt-xtt.min(axis=0)
    xtt = xtt.mean(axis=0)
    xtt0 = xtt0/xtt0.mean()
    xtt = xtt/xtt.mean()
    s1n = xtt.mean()
    s2n = xtt.std()
    s3n = mean((xtt-s1n)**3)/s2n**3
    s4n = mean((xtt-s1n)**4)/s2n**4
    s0t = array([s2,s3,s4,(s3**2+1)/s4,s2n,s3n,s4n,(s3**2+1)/s4])
    xtt = concatenate((s0t,xtt,xtt0))
    return(xtt)
nframes = 1000 

# ~ print(X.shape, Y.shape)
    
pars = load("Training/roival.npy")

lpar = len(pars)

nframes = 1000 
pars = load("Training/roival.npy")

lpar = len(pars)


Xdev = []
Ydev = []
for i in range(lpar-200,lpar):
    xt = load("Training/roi"+str(i).zfill(4)+".npy")
    par = pars[i]
    for j in range(5):
        Xdev.append(getdataextract(xt))
        #Ydev.append(concatenate((yt[start:end],par[-1:])))
        Ydev.append(par)
Xdev = array(Xdev)
Ydev = array(Ydev)

Xdev = reshape(Xdev,(Xdev.shape[0],Xdev.shape[1],Xdev.shape[2],1))
Ydev = reshape(Ydev,(Ydev.shape[0],1,1,1))


def training(Xdev,Ydev,nf,nconv=20,epochs = 1000):
    nframes = 1000
    def model(input_shape):
        '''
        Function used to create the model's graph in Keras
        
        Argument:
        -- input_shape. Shape of the model's input data (Keras conventions?!)
        
        Returns:
        -- model. Keras model instance
        '''
        tf_session = K.get_session()
        
        X_input = Input(shape = input_shape)
        
        # Layers

        X = X_input
        
        n0par = 8
        n1par = 1008
        Xa = Lambda(lambda x: x[:,n1par:,:], output_shape=(n1par-n0par,1))(X)
        Xb = Lambda(lambda x: x[:,n0par:n1par,:], output_shape=(n1par-n0par,1))(X)
        X  = Lambda(lambda x: x[:,:n0par,:], output_shape=(n0par,1))(X)
        X = Reshape((1,n0par))(X)
        Xa = Concatenate(axis=-1)([Xa,Xb])
        
        Xashape = array(input_shape)
        Xashape[-1] -= 1 
        n1, n2, s = (20,10,4)

        Xa = Dropout(0.2)(Xa)
        Xa = Conv1D(100,200,strides=8,padding="same")(Xa)
        #Xc = AveragePooling1D(40,strides=5,padding="valid")(Xa)
        #Xa = MaxPooling1D(40,strides=5,padding="valid")(Xa)
        Xa = Activation("relu")(Xa)  
        #Xa = Concatenate()([Xa,Xb,Xc])
        Xa = BatchNormalization()(Xa)   

        Xa = MaxPooling1D(100,strides=3,padding="valid")(Xa)
        #Xa = Activation("relu")(Xa)  
        #Xa = Concatenate()([Xa,Xb,Xc])
        Xa = BatchNormalization()(Xa)      
        
        Xa = Dropout(0.2)(Xa)
        Xa = Conv1D(200,3,strides=1,padding="same")(Xa)
        Xa = Activation("relu")(Xa)
        Xa = BatchNormalization()(Xa)    
        
        Xa = AveragePooling1D(3,strides=1,padding="same")(Xa)
        #Xa = Activation("relu")(Xa)  
        #Xa = Concatenate()([Xa,Xb,Xc])
        Xa = BatchNormalization()(Xa)      
        
        Xa = Dropout(0.2)(Xa)
        Xa = Conv1D(400,3,strides=1,padding="same")(Xa)
        Xa = Activation("relu")(Xa)
        Xa = BatchNormalization()(Xa)   
        
        Xa = AveragePooling1D(9,strides=1,padding="valid")(Xa)
        #Xa = Activation("relu")(Xa)  
        #Xa = Concatenate()([Xa,Xb,Xc])
        Xa = BatchNormalization()(Xa)   
        
        #Xa = Dropout(0.2)(Xa)
        #Xa = Conv1D(400,5,strides=1,padding="same")(Xa)
        #Xa = Activation("relu")(Xa)
        #Xa = BatchNormalization()(Xa)  
     
        #Xa = Dropout(0.2)(Xa)
        #Xa = AveragePooling1D(4,strides=1,padding="valid")(Xa)

        #Xa = Dropout(0.2)(Xa)
        #Xa = Conv1D(120,10,strides=3,padding="valid")(Xa)
        #Xa = Activation("relu")(Xa)
        #Xa = BatchNormalization()(Xa)  
        
        #Xa = Reshape((1,16*120))(Xa)
        X = Concatenate(axis=2)([X,Xa])
        #X = BatchNormalization()(Xa) 
        
        X = Dropout(0.2)(X)
        X = Dense(200,activation="sigmoid")(X)
        X = BatchNormalization()(X)
        
        X = Dropout(0.2)(X)
        X = Dense(100,activation="relu")(X)
        X = BatchNormalization()(X)

        X = Dropout(0.2)(X)
        X = Dense(40,activation="relu")(X)
        X = BatchNormalization()(X)
        
        X = Dropout(0.2)(X)
        X = Dense(20,activation="relu")(X)
        X = BatchNormalization()(X)
         
        
        X = Dropout(0.2)(X)
        X = Dense(4,activation="sigmoid")(X)
        X = BatchNormalization()(X)
        
        X = Dense(1,activation="sigmoid")(X)

        #X = Xa
        # Defining the model
        
        model = Model(inputs = X_input, outputs = X)
        
        return model
        
    model = model(input_shape = (8+nframes*2,1))\
    opt = Adam(lr=0.001, beta_1=0.95, beta_2=0.999, decay=0.001)

    model.compile(loss='binary_crossentropy', optimizer=opt, metrics=["accuracy"])

    pars = load("Training/roival.npy")

    lpar = len(pars)

    X = []
    Y = []
    for i in permutation(lpar-200):
        xt = load("Training/roi"+str(i).zfill(4)+".npy")
        par = pars[i]
        nj = 3
        if par == 1: nj = 8
        for j in range(nj):
            X.append(getdataextract(xt))
            Y.append(par)
    X = array(X)
    Y = array(Y)
    X = reshape(X,(X.shape[0],X.shape[1],X.shape[2],1))
    Y = reshape(Y,(Y.shape[0],1,1,1))

    history1 = model.fit(X, Y, batch_size = 500, epochs = epochs)
         
    lpar = len(pars)

    nframes = 1000 

    pars = load("/export/home1/users/bssn/serna/Downloads/Classifier/fullimage/moviefull/roival.npy")
    pars = pars[:,2]
    lpar = len(pars)

    X = []
    Y = []
    for i in permutation(lpar-200):
        xt = load("/export/home1/users/bssn/serna/Downloads/Classifier/fullimage/moviefull/roi_F01A"+str(i).zfill(4)+".npy")
        par = pars[i]
        nj = 3
        if par == 1: nj = 3
        for j in range(nj):
            X.append(getdataextract(xt))
            #Y.append(concatenate((yt[start:end],par[-1:])))
            Y.append(par)
    X = array(X)
    Y = array(Y)
    X = reshape(X,(X.shape[0],X.shape[1],X.shape[2],1))
    Y = reshape(Y,(Y.shape[0],1,1,1))
    history2 = model.fit(X, Y, batch_size = 500, epochs = epochs)
            
    loss, acc = model.evaluate(Xdev, Ydev)

    return(loss,acc,model)
    

nfs = [2,4,6,8,10,20,50]
nconvs = [5,10,20,40,80]

pcond = [ [nf,nconv] for nf in nfs for nconv in nconvs]

res = []
for i, np in enumerate(pcond):
    nf,nconv = np
    loss,acc,_  = training(Xdev,Ydev,nf=nf,nconv=nconv,epochs = 500)
    print(" Training with",nf,nconv,"gives",acc)
    res.append([nf,nconv,loss,acc])
    savetxt("temp.dat",array(res))

savetxt("tempF.dat",array(res))


nfo = 4
nconvo = 10
loss,acc,model = training(Xdev,Ydev,nf=nfo,nconv=nconvo,epochs = 1000)

import time
date = time.localtime()
dates = str(date.tm_year)+str(date.tm_mon)+str(date.tm_mday)+'_'+str(date.tm_hour)

model.save("classifier"+dates+".h5")
