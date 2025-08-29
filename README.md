This folder contains files to execute a project in which a short structure will be mutated, folded and then scored.
Initially, it was supposed to be done by FoldX, but it simply refused to run on my Windows or Linux WSL machines.

So I delved into publications and selected alternative software to achieve the same goal. The requirement was the calculation of ddG and the ability to upload multiple mutation lists either to server or installastion.
Based on these requirements, I selected:
- MAESTRO https://pbwww.services.came.sbg.ac.at/maestro/web
- mCSM https://biosig.lab.uq.edu.au/mcsm/ (later upgraded by SDM in DUET)
- PoPMuSiC (problem with licence atm so not included)

Please follow the licence agreement for each software (the data is generated via the web servers)

Literature used for researching tools:
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04238-w#Sec8

The protein software used for visualisation and cleaning is PyMOL (https://pymol.org/2/).

Two files are required to run:
- pdb file placed in Mutation_modelling_complete/pdb/ folder
- region file placed in Mutation_modelling_complete/region_select/ folder (a txt. file where residues to be mutated are listed and separated by spaces in the same line - mine was generated in pyMol around residue 35 within 6A)

At the moment, the code is designed so  that all files will be overwritten for each analysis - no option to give alternative names

It produces 3 types of outputs
- graphs showing each residue and corresponding results from both software
- a graph showing the residues for which both software agreed to stabilising/destabilising effect and when one of the values was above 0.75, they were shown and listed (the standard error for most of softwares die ddG is above 1.1)
- a mutation input file, which is going to be fed to ddG of Rosetta software to further validate the hits

  Comment: Funny enough mCSM and Maestro are in strong disagreement about the majority of the impact of each mutation. So choosing AA, which will produce the same type of change was a very good selector. It gives me some increased hope that those mutations will actually produce some effect.
