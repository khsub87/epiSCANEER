# epiSCNAEER
---

This repository contains the core computational components of **epiSCANEER**.
**epiSCANEER** identifies and recommends residue combinations that improve enzymatic activity by exploiting epistatic relationships inferred from coevolutionary networks.

**epiSCANEER** is based on **SCANEER**, incorporating parts of the original code. \
SCANEER github link: https://github.com/SBIlab/SCANEER

## Requirements
---
- python (v 3.9.5)
- Java (v 1.8.0)
- Dependencies are listed in the requirements.txt file.

##  Usage
---
1. Clone the repository
```
git clone [repository-url]
cd [repository-name]
```
2. Prepare input files
+ Place input files in the ```input/``` directory.
+ Each input file should be a multiple sequence alignment (MSA) of the target enzyme.
+ Input files must be in CLUSTAL format (*.aln).
You may use the example multiple sequence alignments provided in the ```input/``` directory.

3. Configure input/ouput paths(optional)
+ To change the input directory, modify input_path in ```run_SCANEER.py```.
+ To change the output directory, modify output_path in ```run_SCANEER.py```.

4. Run SCANEER
Excute the following command:
```
python run_epiSCANEER.py
```

5. Output files
+ All results will be saved in the ```output/``` directory.
+ The output directory is automatically created if it does not exist.
+ The contents of each output file are as follows:
     + ```*.aln_cn``` - A CLUSTAL format file containing a multiple sequence alignment to calculate covarying strength.
     + ```*.coe_out_mcbasc``` - A text file containing the calculated covarying strengths of all combination of residue pairs.
     + ```*.cn``` - A text file containing the number of co-evolutionary relationships of residues.
     + ```*.coenet``` - A text file containing a residue-residue co-evolutionary network of the enzyme
     + ```*.txt``` - A text file containing final epiSCI scores of mutations.
