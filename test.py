import subprocess
import sys
import numpy as np

EXIT = 0


def run(command):
    result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result.communicate()
    return result.returncode


# Test single input
cmd = "python3.5 pipred.py -i test/test_1/1twu_A.fasta -pssm_path test/test_1/ -out_path test/test_1/"
code = run(cmd)
if code != 0:
    print("Single input test failed!")
    EXIT = 1
results = np.loadtxt('test/test_1/1twu_A.out', skiprows=1, usecols=(1, 2, 3 ,4))
results_bk = np.loadtxt('test/test_1/1twu_A.out.bk', skiprows=1, usecols=(1, 2, 3 ,4))
if not np.array_equal(results, results_bk):
    print("Single input test failed!")
    EXIT = 1

# Test multiple inputs in one file
cmd = "python3.5 pipred.py -i test/test_2/in.fasta -pssm_path test/test_2/ -out_path test/test_2/"
code = run(cmd)
if code != 0:
    print("Multiple inputs test failed!")
    EXIT = 1
results_0 = np.loadtxt('test/test_2/1twu_A.out', skiprows=1, usecols=(1, 2, 3 ,4))
results_0_bk = np.loadtxt('test/test_2/1twu_A.out.bk', skiprows=1, usecols=(1, 2, 3 ,4))
results_1 = np.loadtxt('test/test_2/1mty_C.out', skiprows=1, usecols=(1, 2, 3 ,4))
results_1_bk = np.loadtxt('test/test_2/1mty_C.out.bk', skiprows=1, usecols=(1, 2, 3 ,4))
if not np.array_equal(results_0, results_0_bk):
    print("Multiple inputs test failed!")
    EXIT = 1
if not np.array_equal(results_1, results_1_bk):
    print("Multiple inputs test failed!")
    EXIT = 1

sys.exit(EXIT)

