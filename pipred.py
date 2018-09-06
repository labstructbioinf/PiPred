import argparse
from Bio import SeqIO
import os
import sys
import random
import numpy as np
from utils import enc_seq_onehot, enc_pssm, is_fasta, get_pssm_sequence, PiPred_Model, decode,exit, seq_split, pssm_split
import keras.backend as K
import h5py
# cx_freeze specific
if getattr(sys, 'frozen', False):
    my_loc = os.path.dirname(os.path.abspath(sys.executable))
else:
    my_loc = os.path.dirname(os.path.realpath(__file__))

#my_loc = os.path.dirname(os.path.abspath(__file__))
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

for i in range(1, 4):
    if not os.path.isfile('%s/weights/weights_pipred_%s.h5' % (my_loc, i)):
        print("Weight files for the PiPred model are not available.")
        print("Download weights from http://lbs.cent.uw.edu.pl/")
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
aa1 = "ACDEFGHIKLMNPQRSTVWY"
for entry, seq in zip(entries, sequences):
    if len(seq) < 25:
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
    except ValueError:
        print("ERROR: Malformed PSSM file for entry %s!" % entry)
        exit()
    
    if not parsed_pssm.shape[0] == len(seq) and parsed_pssm.shape[1] == 20:
        if parsed_pssm.shape[0] == len(seq)-2:
            parsed_pssm = np.genfromtxt(pssm_fn, skip_header=3, skip_footer=3, usecols=(i for i in range(2, 22)))
        else:
            print("ERROR: Malformed PSSM file for entry %s!" % entry)
            exit()
    pssm_files.append(pssm_fn)

print("Encoding sequences...")
# Encode sequence into vector format
enc_sequences = []
common = 50
splits_len = []

for seq, pssm_fn in zip(sequences, pssm_files):
    splitted_seq = seq_split(seq,common)
    seq_splitted,n_split, pssm_splitted = splitted_seq[0], splitted_seq[1],pssm_split(pssm_fn, common)
    splits_len.append(n_split)		
    for seq_s, s_pssm in zip(seq_splitted, pssm_splitted):
        enc_sequences.append(np.concatenate((enc_seq_onehot(seq_s, pad_length=700),
					enc_pssm(s_pssm, pad_length=700)), axis=1))


# Create model (supress warnings).
stderr = sys.stderr
sys.stderr = open(os.devnull, 'w')
model = PiPred_Model()
sys.stderr = stderr
enc_sequences = np.asarray(enc_sequences)


ensemble_results = {}
print("Predicting...")
for i in range(1, 4):
    model.load_weights('%s/weights/weights_pipred_%s.h5' % (my_loc, i))
    predictions = model.predict(enc_sequences)[:,:,2:3]
    decoded_predictions = [decode(pred, encoded_seq) for pred, encoded_seq in
                     zip(predictions, enc_sequences)]
    cnt = 0
    for n, n_split in enumerate(splits_len):
            ns = decoded_predictions[cnt].shape[0]
           
            start=0			
            if n_split == 1:
                decoded_prediction = np.zeros((ns,1))
                decoded_prediction+=decoded_predictions[cnt]
            else:
                decoded_prediction = np.array([])
                for j in range(n_split-1):
                    div=min(common,len(decoded_predictions[cnt+j+1]))
                    decoded_prediction = np.concatenate([decoded_prediction,decoded_predictions[cnt+j][:ns-div], (decoded_predictions[cnt+j][ns-div:]+decoded_predictions[cnt+j+1][:div])/2])
                    start = div
                if len(decoded_predictions[n_split-1])>common:	           
                    decoded_prediction = np.concatenate([decoded_prediction,decoded_predictions[n_split-1][common:]])
            entry = entries[n]
            if i == 1:
                ensemble_results[entry] = decoded_prediction
            else:
                ensemble_results[entry] = np.vstack((ensemble_results[entry], decoded_prediction))
            cnt+=n_split
K.clear_session()
# Dump the results
f = h5py.File(args.out_path+'/out', 'w')
for entry, seq in zip(entries, sequences):
    f.create_dataset(data=np.average(ensemble_results[entry], axis=0), name=entry)
f.close()
print("Done!")
