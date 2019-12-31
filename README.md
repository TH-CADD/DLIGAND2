# DLIGAND2
DLIGAND2 is a knowledge-based method to predict protein-ligand binding affinity based on a distance-scaled, finite, ideal-gas reference (DFIRE) state.

<div style="align: center">
<img src="https://github.com/TH-CADD/DLIGAND2/blob/master/doc/picture.png" width="70%" div align=center />
</div>

# How to run 
When you have downloaded the DLIGAND2 software package, the following files can be obtained:
```
$ tree 
.
|-- bin
|   |-- amino.mol2
|   |-- dfire.2
|   |-- dligand2.gnu
|   |-- dligand2.intel
|   `-- Potential_parameters.txt
|-- doc
|   `--pitcure.png
|-- example
|   |-- 966c_ligand.mol2
|   |-- 966c_protein.pdb
|   `-- run_example.sh
|-- LICENSE
`-- README.md

2 directories, 9 files
```
Two version of executable dligand2 are provided, the envrionment of DATAPATH must be set to find the files of amino.mol2 and dfire.2. Then:
```
$ cd example
$ ../bin/dligand2.intel 
Usage: ../bin/dligand2.intel [-v] [-s dfire.2] [-etype 1,2] [-P protein] [-L ligand]
$ ./dligand2.intel -etype 1 -P ../966c/966c_protein.pdb -L ../966c/966c_ligand.mol2
 -12.4143
$ ./dligand2.intel -etype 2 -P ../966c/966c_protein.pdb -L ../966c/966c_ligand.mol2
 -13.908
$ ./dligand2.intel  -P ../966c/966c_protein.pdb -L ../966c/966c_ligand.mol2
 -12.4143
```
The argument of [-etype 1,2] gives two versions of potentials,represents the atom type in protein subdivided into 167 atom types and 13 atom types.

Enjoy!

# Citing
Chen, P., Ke, Y., Lu, Y. et al. DLIGAND2: an improved knowledge-based energy function for protein–ligand interactions using the distance-scaled, finite, ideal-gas reference state. J Cheminform 11, 52 (2019) doi:10.1186/s13321-019-0373-4
