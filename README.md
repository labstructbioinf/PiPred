# **PiPred - deep learning based prediction of pi-helices in protein sequences** #

The webserver is available at the MPI Bioinformatics Toolkit (Quick2D tool) - https://toolkit.tuebingen.mpg.de/tools/quick2d

The original paper can be accessed here - https://www.nature.com/articles/s41598-019-43189-4

## **Installation** ##
First clone this repository:
```bash
$ git clone https://github.com/labstructbioinf/PiPred.git
$ cd PiPred
```

Required packages to run PiPred are listed in the **`requirements.txt`** file.
We suggest running PiPred in the virtual environment:
If you don't have virtualenv installed do so:
```bash
$ pip3 install virtualenv
```
Create virtual environment and install required packages:
```bash
$ cd virtual_envs_location
$ virtualenv pipred_env
$ source pipred_env/bin/activate
$ cd PIPRED_LOCATION
$ pip3 install -r requirements.txt
```
Test the installation:
```bash
$ ./run_example.sh
```
This should produce output **`example/1mty_C.out`** identical to **`example/1mty_C.out.bk`**.

## **Usage** ##
```bash
python3.5 pipred.py [-h] -i FILE [-out_path DIR] [-pssm_path DIR]
```
| Option    | Description |
|:----------:|-------------|
| **`-i`** | Input file in FASTA format. Can contain multiple entries. |
| **`-pssm_path`** | Directory with psiblast PSSM files. For each entry in the input fasta file there must be a PSSM file. |
| **`-out_path`** | Directory where the predictions are saved. For each entry one file will be saved. |

PSSM filenames should be based on the identifiers in the fasta file (only alphanumeric characters and '_'). For example if a fasta sequence is as follows:
```
>1mty_C
ERRRGLTDPEMAAV...
```
PSSM file should be named **`1mty_C.pssm`**.

You can generate PSSM files with the following command (requires NR90 database):
```bash
psiblast -query 1mty_C.fasta -db NR90_LOCATION -evalue 0.001 -num_iterations 3 -out_ascii_pssm 1mty_C.pssm
```
In order to generate PSSM file from multiple sequence alignment (MSA) you can use this command:
```bash
psiblast -subject sequence.fasta -in_msa alignment.fasta -out_ascii_pssm output.pssm
```

## **Supplementary files** ##

* Updated CB6133 and CB513 benchmark datasets - [link](https://lbs.cent.uw.edu.pl/pipred)
* Results of PFAM scan - [link](https://lbs.cent.uw.edu.pl/pipred)


