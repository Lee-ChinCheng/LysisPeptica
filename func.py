import numpy as np
import sys
#from tensorflow.keras.metrics import *
import tensorflow as tf






def pc6_encode(seq, pram):  #return 6depth

    pc6_zs_d={
        'A': [0.62, -1.239, -0.084, -0.014, 0.7, -0.441], 'C': [0.29, -0.768, -1.05, -0.539, -0.932, -1.115], 
        'D': [-0.9, -0.895, 1.738, -1.839, -1.276, -0.918], 'E': [-0.74, -0.29, 1.477, -1.585, 0.056, -0.447], 
        'F': [1.19, 1.181, -1.161, -0.308, -1.49, 0.026], 'G': [0.48, -1.995, 0.251, -0.031, 0.7, 2.202], 
        'H': [-0.4, 0.178, 0.771, 0.885, -1.533, -0.716], 'I': [1.38, 0.576, -1.161, -0.003, 0.786, -0.219], 
        'K': [-1.5, 0.755, 1.106, 2.1, 0.013, -0.279], 'L': [1.06, 0.576, -1.273, -0.025, 0.786, 0.243], 
        'M': [0.64, 0.593, -0.976, -0.161, 0.442, -0.51], 'N': [-0.78, -0.381, 1.217, -0.347, -0.674, -0.469], 
        'P': [0.12, -0.843, -0.121, 0.156, -0.803, 3.132], 'Q': [-0.85, 0.224, 0.808, -0.212, -0.03, 0.205], 
        'R': [-2.531, 0.892, 0.808, 2.676, -0.03, 0.119], 'S': [-0.18, -1.189, 0.325, -0.195, 0.142, -0.481], 
        'T': [-0.05, -0.584, 0.102, -0.24, -0.374, -0.5], 'V': [1.08, -0.029, -0.901, -0.036, 0.614, 0.325], 
        'W': [0.81, 2.006, -1.087, -0.076, 2.805, 0.032], 'Y': [0.26, 1.231, -0.79, -0.206, 0.099, -0.189]}
    pc6_norm_d={
        'A': [0.806, 0.189, 0.395, 0.404, 0.515, 0.159], 'C': [0.721, 0.307, 0.074, 0.288, 0.139, 0.0], 
        'D': [0.417, 0.275, 1.0, 0.0, 0.059, 0.046], 'E': [0.458, 0.426, 0.913, 0.056, 0.366, 0.157], 
        'F': [0.951, 0.794, 0.037, 0.339, 0.01, 0.269], 'G': [0.77, 0.0, 0.506, 0.4, 0.515, 0.781], 
        'H': [0.545, 0.543, 0.679, 0.603, 0.0, 0.094], 'I': [1.0, 0.643, 0.037, 0.407, 0.535, 0.211], 
        'K': [0.264, 0.687, 0.79, 0.872, 0.356, 0.197], 'L': [0.918, 0.643, 0.0, 0.402, 0.535, 0.32], 
        'M': [0.811, 0.647, 0.099, 0.372, 0.455, 0.142], 'N': [0.448, 0.403, 0.827, 0.33, 0.198, 0.152], 
        'P': [0.678, 0.288, 0.383, 0.442, 0.168, 1.0], 'Q': [0.43, 0.555, 0.691, 0.36, 0.346, 0.311], 
        'R': [0.0, 0.722, 0.691, 1.0, 0.346, 0.291], 'S': [0.601, 0.201, 0.531, 0.364, 0.386, 0.149], 
        'T': [0.634, 0.353, 0.457, 0.354, 0.267, 0.145], 'V': [0.923, 0.491, 0.124, 0.399, 0.495, 0.339], 
        'W': [0.854, 1.0, 0.062, 0.39, 1.0, 0.27], 'Y': [0.714, 0.806, 0.16, 0.362, 0.376, 0.218]}
    
    if pram=='pc6zs':   
        pc6_d=pc6_zs_d  
    else: #pram=='pc6norm'               
        pc6_d=pc6_norm_d 

    encoded_seq = []
    for base in seq:
        if base not in pc6_d: print(base,'no amino acid one hot')
        encoded_seq.append(pc6_d.get(base, [0, 0, 0, 0, 0, 0]))  
    return np.array(encoded_seq)#.T



def ugml_to_uM(seq:str, conc: float) -> float:
    """
    input seq & conc(ug/ml) -> output  conc(uM)
    """
    # Average amino acid residue masses (in Da)
    aa_weights = {
        'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
        'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
        'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
        'M': 131.19, 'F': 147.18, 'P': 97.12,  'S': 87.08,
        'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.14
    }
    # get molecular weight
    mw = sum(aa_weights[aa] for aa in seq) + 18.015  # add H2O 

    conc_uM = (conc * 1000) / mw
    return round(conc_uM, 4)



def pc6_8d_encode( inpli, ugmlli, pram, tar_len ):

    #conc_norm_d = {'uM' : {'Min_v': 0.08, 'Max_v': 128 }, 'ugml' : {'Min_v': 0.2,  'Max_v': 250}} #250 is pr90
    #conc_zs_d   = {'uM' : {'mean': 53 , 'std': 41 },      'ugml' : {'mean': 120,  'std': 87}}
    Min_ugml, Max_ugml = 0.2, 250
    Min_uM, Max_uM     = 0.08, 128
    ugml_mean, ugml_std = 120, 87
    uM_mean, uM_std     = 53, 41
    
    op_li=[]
    for idx,seq in enumerate(inpli):
        ugml = ugmlli[idx]
        uM = ugml_to_uM(seq, ugml)
        en1_seq = pc6_encode(seq, pram) #pc6zs

        # get normed(zs or min-max)  2 conc(ugml & uM)
        if pram=='pc6zs':
            if    ugml == ugml_mean: nm_ugml=0
            else:   nm_ugml = round( ((ugml - ugml_mean) / ugml_std) , 4 )

            if    uM == uM_mean: nm_uM=0
            else:   nm_uM = round( ((uM - uM_mean) / uM_std) , 4 )

        else: #pram=='pc6norm':
            if   ugml <= Min_ugml:  nm_ugml = Min_ugml
            elif ugml >= Max_ugml:  nm_ugml = Max_ugml 
            else:   nm_ugml = round( ((ugml - Min_ugml) / (Max_ugml - Min_ugml)) , 4)  

            if   uM <= Min_uM:  nm_uM = Min_uM
            elif uM >= Max_uM:  nm_uM = Max_uM
            else:   nm_uM = round( ((uM - Min_uM) / (Max_uM - Min_uM)) , 4) 
   
        # Create a column of normed_conc   
        col1 = np.full((en1_seq.shape[0], 1), nm_ugml)
        col2 = np.full((en1_seq.shape[0], 1), nm_uM)
        en2_seq = np.hstack((en1_seq, col1, col2))           
                                  
        # Calculate how many rows to add
        if len(seq) < tar_len:
            rows_to_add = tar_len - en2_seq.shape[0]
            pad_arr = np.pad(en2_seq, ((0, rows_to_add), (0, 0)), mode='constant', constant_values=0)                                                 
        elif len(seq) == tar_len:
            pad_arr = en2_seq
        else:
            sys.exit("seq over target(max) length")   

        op_li.append(pad_arr)   

    op_li = np.array(op_li)
    return op_li   

 

def add_conc_on_pepbert_array(inpli, seqli, ugmlli):
    #Min_ugml, Max_ugml = 0.2, 250 
    Min_uM, Max_uM = 0.02, 128
    pbert_c_min, pbert_c_max, Shift = -2.3, 2.3, 2.3 #-3,3,3

    op_li=[]
    for idx,arr in enumerate(inpli):
        ugml = ugmlli[idx]
        uM = ugml_to_uM( seqli[idx], ugml )
        
        #if ugml <= Min_ugml:    nm_ugml = pbert_c_min 
        #elif ugml >= Max_ugml:  nm_ugml = pbert_c_max 
        #else:
        #    times=2*Shift
        #    nm_ugml = round( (  times*(ugml - Min_ugml) / (Max_ugml - Min_ugml) - Shift  ) , 3)
        if uM <= Min_uM:    nm_uM = pbert_c_min 
        elif uM >= Max_uM:  nm_uM = pbert_c_max 
        else:
            times=2*Shift
            nm_uM = round( (  times*(uM - Min_uM) / (Max_uM - Min_uM) - Shift  ) , 3)

   
        col = np.full((arr.shape[0], 1), nm_uM)
        # Concatenate along axis=1 (horizontally)
        arr161d = np.hstack((arr, col))
        op_li.append(arr161d)
        
    op_li=np.array(op_li)
    return op_li


class CustomModel(tf.keras.Model):

    def __init__(self, **kwargs):

        super(CustomModel, self).__init__(**kwargs)



class GlobalMinPooling1D(tf.keras.layers.Layer):
    def __init__(self, **kwargs):
        super(GlobalMinPooling1D, self).__init__(**kwargs)

    def call(self, inputs):
        # reduce along the time/sequence axis
        return tf.reduce_min(inputs, axis=1)

    def compute_output_shape(self, input_shape):
        # input_shape = (batch_size, steps, features)
        return (input_shape[0], input_shape[2])

    def get_config(self):
        config = super(GlobalMinPooling1D, self).get_config()
        return config

