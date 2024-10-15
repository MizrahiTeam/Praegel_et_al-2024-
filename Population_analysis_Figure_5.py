import numpy as np
import os
import scipy.io as spio
from multiprocessing import Pool
from functools import partial
import scipy
from sklearn.cluster import KMeans
from tqdm import tqdm_notebook as tqdm
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.stats
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from numpy import dot
from numpy.linalg import norm
import pickle
from sklearn.model_selection import train_test_split
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import accuracy_score

# from https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    def _check_keys(d):
        '''
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        '''
        for key in d:
            if isinstance(d[key], spio.matlab.mio5_params.mat_struct):
                d[key] = _todict(d[key])
            elif isinstance(d[key], np.ndarray):
                d[key] = _tolist(d[key])
        return d

    def _todict(matobj):
        '''
        A recursive function which constructs from matobjects nested dictionaries
        '''
        d = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, spio.matlab.mio5_params.mat_struct):
                d[strg] = _todict(elem)
            elif isinstance(elem, np.ndarray):
                d[strg] = _tolist(elem)
            else:
                d[strg] = elem
        return d

    def _tolist(ndarray):
        '''
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        '''
        elem_list = []
        for sub_elem in ndarray:
            if isinstance(sub_elem, spio.matlab.mio5_params.mat_struct):
                elem_list.append(_todict(sub_elem))
            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(_tolist(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list
    data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def bin_raw(spike_times,aligned_time,sampled_bins, window_size):
    spike_times=np.array(spike_times)
    sp=spike_times-aligned_time
    sp=sp[np.logical_and(sp>-window_size,sp<window_size)]
    return np.histogram(sp,sampled_bins)[0]

def bin_raw_multievents(aligned_times,sampled_bins,window_size,spike_times):
    return_list=[]
    for aligned_time in aligned_times:
        return_list.append(bin_raw(spike_times,aligned_time,sampled_bins,window_size))
    return np.array(return_list)

def bin_raw(spike_times,aligned_time,sampled_bins, window_size):
    spike_times=np.array(spike_times)
    sp=spike_times-aligned_time
    sp=sp[np.logical_and(sp>-window_size,sp<window_size)]
    return np.histogram(sp,sampled_bins)[0]

def bin_raw_multievents(aligned_times,sampled_bins,window_size,spike_times):
    return_list=[]
    for aligned_time in aligned_times:
        return_list.append(bin_raw(spike_times,aligned_time,sampled_bins,window_size))
    return np.array(return_list)

class session:
    def __init__(self,spike_list,event_list):
        self.spike_list=spike_list
        self.event_list=event_list
    def count_event_times(self,event):
        return len(self.event_list[event])
    def binning(self,sampled_bins,event,sigma=0.005,window_size=10,cpu=32, subsample=None, subsample_event=None):
        if subsample_event is None:
            tem_func=partial(Gaussian_smooth_multievents,self.event_list[event],sampled_bins,sigma,window_size)
        else:
            event=self.event_list[event]
            tempt=np.array(event)
            np.random.shuffle(tempt)
            tem_func=partial(Gaussian_smooth_multievents,tempt[:subsample_event],sampled_bins,sigma,window_size)
        spike_list=np.array(self.spike_list)
        if subsample is not None:
            np.random.shuffle(spike_list)
        if subsample is None:
            subsample=len(self.spike_list)
        if cpu==1:  
            return_list=[]
            for neuron in spike_list[:subsample]:
                return_list.append(tem_func(neuron))
        else:
            with Pool(cpu) as p:
                return_list=p.map(tem_func,spike_list[:subsample])
        return np.array(return_list)
    def binning_raw(self,sampled_bins,event,window_size=10,cpu=32, subsample=None, subsample_event=None):
        if subsample_event is None:
            tem_func=partial(bin_raw_multievents,self.event_list[event],sampled_bins,window_size)
        else:
            event=self.event_list[event]
            tempt=np.array(event)
            np.random.shuffle(tempt)
            tem_func=partial(bin_raw_multievents,tempt[:subsample_event],sampled_bins,window_size)
        spike_list=np.array(self.spike_list)
        if subsample is None:
            subsample=len(self.spike_list)
        else:
            indices=np.arange(len(spike_list))
            np.random.shuffle(indices)
            spike_list=spike_list[indices]
        if cpu==1:  
            return_list=[]
            for neuron in spike_list[:subsample]:
                return_list.append(tem_func(neuron))
        else:
            with Pool(cpu) as p:
                return_list=p.map(tem_func,spike_list[:subsample])
        return np.array(return_list)

def creating_dataset(s,event_list,sampled_bins,window_size,cpu,subsample,subsample_event):
    x_list=[]
    y_list=[]
    for i,event in enumerate(event_list):
        x=s.binning_raw(sampled_bins,event,window_size,cpu,None,subsample_event)
        x=np.swapaxes(x,0,1)
        x_list.append(x)
        y_list.append(np.ones(len(x))*i)
    x_list=np.concatenate(x_list)
    indices=np.arange(x_list.shape[1])
    np.random.shuffle(indices)
    x_list=x_list[:,indices]
    if subsample is None:
        return x_list,np.concatenate(y_list)
    else:
        return x_list[:,:subsample],np.concatenate(y_list)

def cos_sim(a,b):
    return dot(a, b)/(norm(a)*norm(b))
    
def lda_distance(args,event_list,bin_min,bin_max,n, window_size=50,cpu=1,subsample=20,subsample_event=20):
    
    s=args[0]
    seed=args[1]
    np.random.seed(seed)
    lda=LinearDiscriminantAnalysis()
    results=[]
    sampled_bin=np.linspace(bin_min,bin_max,n+1)
    x,y=creating_dataset(s,event_list,sampled_bin,window_size,cpu,subsample,subsample_event)
    x=x.reshape(x.shape[0],-1)
    x=(x-x.mean(axis=0))/(x.std(axis=0)+1e-10)
    
    new_x,new_y=x[np.logical_or(y==0,y==1)],y[np.logical_or(y==0,y==1)]
    lda.fit(new_x,new_y)
    x_trans=lda.transform(new_x)
    x_trans=x_trans/np.linalg.norm(lda.scalings_)
    results.append([lda.means_, np.array((x_trans[new_y==0].mean(), x_trans[new_y==1].mean())), np.array((x_trans[new_y==0].var(), x_trans[new_y==1].var()))])
    vec1=lda.means_[0]-lda.means_[1]
    
    new_x,new_y=x[np.logical_or(y==2,y==3)],y[np.logical_or(y==2,y==3)]
    lda.fit(new_x,new_y)
    x_trans=lda.transform(new_x)
    x_trans=x_trans/np.linalg.norm(lda.scalings_)
    results.append([lda.means_, np.array((x_trans[new_y==2].mean(), x_trans[new_y==3].mean())), np.array((x_trans[new_y==2].var(), x_trans[new_y==3].var()))])
    vec2=lda.means_[0]-lda.means_[1]
    
    results.append(cos_sim(vec1,vec2))
    
    return results

def create_balanced_splits(lda, X, y, test_size_per_class, seed=42, cv=5):
    np.random.seed(seed)
    scores = []
    tempt=[]
    tempt2=[]
    for _ in range(cv):
        # Initialize containers for the split
        X_train_list, X_val_list, y_train_list, y_val_list = [], [], [], []

        for cls in np.unique(y):
            # Separate the dataset by class
            X_cls = X[y == cls]
            y_cls = y[y == cls]

            # Calculate the number of validation samples for this class
            val_size = min(len(X_cls), test_size_per_class)

            # Split the class-specific data
            X_cls_train, X_cls_val, y_cls_train, y_cls_val = train_test_split(X_cls, y_cls, test_size=val_size,random_state=cv*seed+_)

            # Append class-specific splits to the lists
            X_train_list.append(X_cls_train)
            X_val_list.append(X_cls_val)
            y_train_list.append(y_cls_train)
            y_val_list.append(y_cls_val)

        # Combine class-specific splits into training and validation sets
        X_train = np.concatenate(X_train_list)
        y_train = np.concatenate(y_train_list)
        X_val = np.concatenate(X_val_list)
        y_val = np.concatenate(y_val_list)


        # Train and evaluate the model
        lda.fit(X_train, y_train)
        y_pred = lda.predict(X_val)
        score = accuracy_score(y_val, y_pred)
        scores.append(score)

    return scores#,tempt,tempt2

def lda_acc(args,event_list,bin_min,bin_max,n, window_size=50,cpu=1,subsample=20,subsample_event=20,cv=5,test_size_per_class=10):
    s=args[0]
    seed=args[1]
    np.random.seed(seed)
    #lda=LinearDiscriminantAnalysis(solver='eigen',shrinkage='auto',priors=[0.5,0.5])
    lda=LinearDiscriminantAnalysis(priors=[0.5,0.5])
    results=[]
    sampled_bin=np.linspace(bin_min,bin_max,n+1)
    x,y=creating_dataset(s,event_list,sampled_bin,window_size,cpu,subsample,subsample_event)
    x=x.reshape(x.shape[0],-1)
    x=(x-x.mean(axis=0))/(x.std(axis=0)+1e-10)
    
    new_x,new_y=x[np.logical_or(y==0,y==1)],y[np.logical_or(y==0,y==1)]
    results.append(create_balanced_splits(lda, new_x, new_y, test_size_per_class, seed=seed, cv=cv))
    
    new_x,new_y=x[np.logical_or(y==2,y==3)],y[np.logical_or(y==2,y==3)]
    results.append(create_balanced_splits(lda, new_x, new_y, test_size_per_class, seed=seed, cv=cv))
    return results

def gs(A):
    A=A.copy()    
    (n, m) = A.shape    
    for i in range(n):        
        q = A[i] # i-th column of A        
        for j in range(i):
            q = q - np.dot(A[j], A[i]) * A[j]        
        if np.array_equal(q, np.zeros(q.shape)):
            raise np.linalg.LinAlgError("The column vectors are not linearly independent")        
        # normalize q
        q = q / np.sqrt(np.dot(q, q))        
        # write the vector back in the matrix
        A[i] = q    
    return A

def lda_variance(args,event_list,bin_min,bin_max,n, window_size=50,cpu=1,subsample=20,subsample_event=20):
    s=args[0]
    seed=args[1]
    np.random.seed(seed)
    lda=LinearDiscriminantAnalysis()
    results=[]
    sampled_bin=np.linspace(bin_min,bin_max,n+1)
    x,y=creating_dataset(s,event_list,sampled_bin,window_size,cpu,subsample,subsample_event)
    x=x.reshape(x.shape[0],-1)
    x=(x-x.mean(axis=0))/(x.std(axis=0)+1e-10)
    results.append([x[y==ii].var(axis=0).mean() for ii in range(4)])
    lda.fit(x,y)
    x_trans=lda.transform(x)
    results.append([x_trans[y==ii].var(axis=0).mean() for ii in range(4)])
    
    return results

def lda_variance_vis(args,event_list,bin_min,bin_max,n, window_size=50,cpu=1,subsample=20,subsample_event=20):
    
    s=args[0]
    seed=args[1]
    np.random.seed(seed)
    lda=LinearDiscriminantAnalysis()
    results=[]
    sampled_bin=np.linspace(bin_min,bin_max,n+1)
    x,y=creating_dataset(s,event_list,sampled_bin,window_size,cpu,subsample,subsample_event)
    x=x.reshape(x.shape[0],-1)
    x=(x-x.mean(axis=0))/(x.std(axis=0)+1e-10)
    
    new_x,new_y=x[np.logical_or(y==0,y==1)],y[np.logical_or(y==0,y==1)]
    lda.fit(new_x,new_y)
    x_trans=lda.transform(x)
    x_trans=x_trans/np.linalg.norm(lda.scalings_)
    vec1=lda.scalings_.ravel()/np.linalg.norm(lda.scalings_)
    
    new_x,new_y=x[np.logical_or(y==2,y==3)],y[np.logical_or(y==2,y==3)]
    lda.fit(new_x,new_y)
    x_trans=lda.transform(x)
    x_trans=x_trans/np.linalg.norm(lda.scalings_)
    vec2=lda.scalings_.ravel()/np.linalg.norm(lda.scalings_)

    
    vec_mid=(vec2+vec1)/2
    vec_mid=vec_mid/np.linalg.norm(vec_mid)
    
    scalings=gs(np.array([vec_mid,vec1])).T
    x_trans=(x-lda.xbar_)@scalings
    
    results.append([(x_trans[y==i]).mean(axis=0) for i in range(4)])
    results.append([np.cov(x_trans[y==i].T) for i in range(4)])
    results.append(np.array([vec1,vec2])@scalings)
   
    
    return results