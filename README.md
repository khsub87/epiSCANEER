# epiSCNAEER
---

This repository contains the core computational components of **epiSCANEER**.
**epiSCANEER** identifies and recommends residue combinations that improve enzymatic activity by exploiting epistatic relationships inferred from coevolutionary networks.

**epiSCANEER** is based on **SCANEER**, incorporating parts of the original code. \
SCANEER github link: https://github.com/SBIlab/SCANEER

## Requirements
---
- python (v 2.7.13)
- Biopython (v 1.72)
- Numpy (v 1.15.4)
- Java (v 1.7.0)

##Usage##
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
     
