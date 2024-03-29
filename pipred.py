import argparse
from Bio import SeqIO
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import sys
import random
import numpy as np
from utils import enc_seq_onehot, enc_pssm, is_fasta, get_pssm_sequence, PiPred_Model, decode,exit
import tensorflow as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
import tensorflow.keras.backend as K



my_loc = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description='PiPred')
parser.add_argument('-i',
                    help='Input file with sequence in fasta format.',
                    required=True,
                    metavar='FILE')
parser.add_argument('-out_path',
                    help='Output directory',
                    default='.',
                    metavar='DIR')
parser.add_argument('-pssm_path',
                    metavar='DIR',
                    default='.',
                    help='Directory with PSSM files.')
args = parser.parse_args()


# Verify whether weights files are present

for i in range(1, 11):
    if not os.path.isfile('%s/weights/final_%s.h5' % (my_loc, i)):
        print("Weight files for the PiPred model are not available in weights directory.")
        exit()


# INPUT VERIFICATION #

print("Veryfing input...")
# Check if input file exists
if not os.path.isfile(args.i):
    print('ERROR: Input file does not exist!')
    exit()
# Check if input is valid fasta file
if not is_fasta(args.i):
    print("ERROR: Malformed fasta file. Please check input!")
    exit()
if not os.path.isdir(args.out_path):
    print("ERROR: Output directory does not exist!")
    exit()
# Parse fasta file
input_data = list(SeqIO.parse(args.i, "fasta"))
sequences = [str(data.seq) for data in input_data]
entries = [''.join(e for e in str(data.id) if (e.isalnum() or e == '_')) for data in input_data]
if not len(entries) == len(set(entries)):
    print("ERROR: Sequence identifiers in the fasta file are not unique!")
    exit()
# Check sequence length and presence of non standard residues
aa1 = "ACDEFGHIKLMNPQRSTVWYX"
for entry, seq in zip(entries, sequences):
    if not (len(seq) >= 25 and len(seq) <= 700):
        print('ERROR: Not accepted sequence length (ID %s - %s). Only sequences between 30 and 700 residues are accepted!' % (
        entry, len(seq)))
        exit()
    for aa in seq:
        if aa not in aa1:
            print("ERROR: Sequence (ID %s) contains non-standard residue (%s)." % (entry, aa))
            exit()

# PSSM SPECIFIC INPUT VERIFICATION #
pssm_files = []
# Check if directory exists
if not os.path.isdir(args.out_path):
    print("ERROR: Directory with PSSM files does not exist!")
    exit()
for entry, seq in zip(entries, sequences):
    pssm_fn = '%s/%s.pssm' % (args.pssm_path, entry)
    if not os.path.isfile(pssm_fn):
        print("ERROR: PSSM file for entry %s does not exist!" % entry)
        exit()
    if not get_pssm_sequence(pssm_fn) == seq:
        print("ERROR: Sequence in PSSM file does not match fasta sequence for entry %s!" % entry)
        exit()
    try:
        parsed_pssm = np.genfromtxt(pssm_fn, skip_header=3, skip_footer=5, usecols=(i for i in range(2, 22)))
        if parsed_pssm.shape[0] == len(seq) - 2:
            parsed_pssm = np.genfromtxt(pssm_fn, skip_header=3, skip_footer=3, usecols=(i for i in range(2, 22)))
        elif parsed_pssm.shape[0] != len(seq):
            raise ValueError
    except ValueError:
        print("ERROR: Malformed PSSM file for entry %s!" % entry)
        exit()
    pssm_files.append(pssm_fn)

print("Encoding sequences...")
# Encode sequence into vector format
enc_sequences = []
for seq, pssm_fn in zip(sequences, pssm_files):
    pad_left = random.randint(0, 700 - len(seq))
    enc_sequences.append(np.concatenate((enc_seq_onehot(seq, pad_length=700),
                                         enc_pssm(pssm_fn, pad_length=700)), axis=1))
    
# Create model 
model = PiPred_Model()


enc_sequences = np.asarray(enc_sequences)


ensemble_results = {}
print("Predicting...")
for i in range(1, 11):
    model.load_weights('%s/weights/final_%s.h5' % (my_loc, i))
    predictions = model.predict(enc_sequences)
    decoded_predictions = [decode(pred, encoded_seq) for pred, encoded_seq in
                     zip(predictions, enc_sequences)]
    for decoded_prediction, entry in zip(decoded_predictions, entries):
        if i == 1:
            ensemble_results[entry] = []
        k = 0
        for aa_col in decoded_prediction.T:
            if i == 1:
                ensemble_results[entry].append(aa_col)
            else:
                ensemble_results[entry][k] = np.vstack((ensemble_results[entry][k], aa_col))
            k += 1
K.clear_session()

# Dump the results
for entry, seq in zip(entries, sequences):
    f = open('%s/%s.out' % (args.out_path, entry), 'w')
    f.write('aa prH   prI   prE   prC\n')
    final_results = np.asarray([np.average(ss_col, axis=0) for ss_col in ensemble_results[entry]])
    for aa, probs in zip(seq, final_results.T):
        f.write("%s " % aa)
        for prob in probs:
             f.write(" %.3f" % prob)
        f.write("\n")
    f.close()
print("Done!")
