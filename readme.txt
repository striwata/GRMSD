G-RMSD: Program code to determine molecular similarity in 3D


@ Folders:
-------------------------------------------

algo   : program codes for calculating RMSD
query  : query structures in csv format. All the molecules in the same file are composed of the same number of atoms.
target : dataset structures in csv format. All the molecules in the same file are composed of the same number of atoms.
sdf2csv   : program code for converting sdf to csv files.



@ Input and output of calculate_RMSDs.m
--------------------------------------------

input:
 csv files for query and target (see below for the format)
 parameters (optional)

output:
 q_t_result.csv
 (q: query file name, t: target file name)



Input file (csv)
------------------
query and target files include molecular information in csv.

The number of atoms in a molecule: n
The first line: n,n,n,n
From the second line to n+1-th line: x-coordinate, y-coordinate, z-coordinate, atomic number


Output file (csv)
------------------
The results are output to q_t_result.csv
The value of (i,j) is the RMSD value between i-th target and j-th query.
The value "Inf" is written if it was error.

@ Compile C++ programs
--------------------------------------------

Those two C++ codes need to compile:

algo/sub/Hungarian2.cpp
algo/sub/Hungarian3.cpp


@ Parameters
--------------------------------------------
options : a structure array with attributes of e_label,permit_mirror,ignore_atom
The default values are used if not specified.

use_label     : Consider atom types (true) or not (false) in the RMSD calculation. Default is true.
               
permit_mirror : Permit chirality (true) or not (false) in the RMSD calculation. Default is false.

ignore_atom   : Atomic numbers to omit in the RMSD calculation. Default is 1 (:H atoms are omitted). Set [] (null) to consider H atoms. Use comma to specify multiple atom types (ex: [1,8] to omit H and O atoms.

clus_mode     : If true, RMSD is calculated when the query and target data has the same file name. Default is false.

reduction     : If true, saving the steps in the RMSD calculation. Default is true. 

iter_num      : How many initial values are specified if reduction is true. Default is 4.

