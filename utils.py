from Bio import SeqIO
import os
import sys
import numpy as np

stderr = sys.stderr
sys.stderr = open(os.devnull, 'w')
from keras.layers import LSTM, Convolution1D, Dropout, BatchNormalization, TimeDistributed, Bidirectional, Input, merge, \
    Dense
from keras.regularizers import l2
from keras.models import Model
import keras.backend as K

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
sys.stderr = stderr


# Adapted from https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta
# Checks whether file is in fasta format
def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

# Reads sequence from Psiblast PSSM file
def get_pssm_sequence(fn):
    c = 0
    seq_list = []
    try:
        with open(fn) as f:
            for line in f:
                if c > 2:
                    try:
                        aa = line.split()
                        seq_list.append(aa[1])
                    except IndexError:
                        break
                c += 1
        f.close()
    except FileNotFoundError:
        pass
    return ''.join(seq_list)


def sigmoid(x):
    return 1 / (1 + np.exp(-x))

# Encodes amino acid sequence in one-hot format
def enc_seq_onehot(seq, pad_length=None, pad_left=0):
    aa1 = list("ACDEFGHIKLMNPQRSTVWY")
    aa_indices = {aa1[k]: k for k in range(0, len(aa1))}
    enc_seq = []
    for aa in seq:
        enc_aa = np.zeros(20)
        enc_aa[aa_indices[aa]] = 1
        enc_seq.append(enc_aa)
    matrix = np.asarray(enc_seq)
    if pad_length:
        pad_matrix = np.zeros((pad_length, 20))
        pad_matrix[pad_left:matrix.shape[0] + pad_left, 0:matrix.shape[1]] = matrix
        return pad_matrix
    return matrix

# Encodes PSSM
def enc_pssm(pssm_file, pad_length=None, pad_left=0):
    pssm_matrix = sigmoid(np.genfromtxt(pssm_file, skip_header=3, skip_footer=5, usecols=(i for i in range(2, 22))))
    if pad_length:
        pad_matrix = np.zeros((pad_length, 20))
        pad_matrix[pad_left:pssm_matrix.shape[0] + pad_left, 0:pssm_matrix.shape[1]] = pssm_matrix
        return pad_matrix
    return pssm_matrix

# Decodes predictions (takes into the account padding of sequence)
def decode(pred, enc_sec):
    return pred[np.any(enc_sec, axis=-1), :]

# Pipred model architecture
def PiPred_Model():
    inp1 = Input(shape=(700, 40), dtype='float32', name='inp')
    a = Convolution1D(64, 3, activation='tanh', padding='same', kernel_regularizer=l2(0.0001))(inp1)
    a_b = BatchNormalization()(a)
    b = Convolution1D(64, 5, activation='tanh', padding='same', kernel_regularizer=l2(0.0001))(inp1)
    b_b = BatchNormalization()(b)
    e = Convolution1D(64, 7, activation='tanh', padding='same', kernel_regularizer=l2(0.0001))(inp1)
    e_b = BatchNormalization()(e)
    x = merge([a_b, b_b, e_b], mode='concat', concat_axis=-1)
    t = TimeDistributed(Dense(200, activation='relu', kernel_regularizer=l2(0.0001)))(x)
    k = Bidirectional(LSTM(200, return_sequences=True, activation='tanh', recurrent_activation='sigmoid', dropout=0.5,
                           recurrent_dropout=0.5))(t)
    k1 = Bidirectional(LSTM(200, return_sequences=True, activation='tanh', recurrent_activation='sigmoid', dropout=0.5,
                            recurrent_dropout=0.5))(k)
    f = TimeDistributed(Dense(200, activation='relu', kernel_regularizer=l2(0.0001)))(k1)
    out = TimeDistributed(Dense(4, activation='softmax', name='out'))(f)
    model = Model(inputs=inp1, outputs=out)
    model.compile(optimizer='adam', loss='categorical_crossentropy')
    return model

# Exit function
def exit():
    print("Run failed!")
    K.clear_session()
    sys.exit(1)
