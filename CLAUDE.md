1) here is molscat 2020 package coded in old-style fortran 
2) please ignoret all .o files produced during compilation
3) in the sub-directory examples/  one can find simplest examplest for running molscat. 
4) in the root directory you cna find example run_molscat_sro2_example.sh which was used to submit jobs into SLURM system of present cluster. Please adopt it to the case of atom-diatom case. In the source_code there are numerous molscat_*   files compiled on cluster (using intel and blas). Do not try to compile codes unless I am in the interactive mode 
5) read potential_generator_plan.md in which I summarized my work on running molscat test for potential sampling.
6) I also vibe coded potential_generator.py 
7)  Initial test can be done on set of random potentials which are already in the file potentials_library.csv with parameters defined in  potential_generator_plan.md 

8)molscat runs should be performed for log-uniform grid between 0.1 and 200 cm-1 

