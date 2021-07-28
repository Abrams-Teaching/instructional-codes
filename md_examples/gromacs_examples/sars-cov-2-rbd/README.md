## The SARS-CoV-2 S Receptor Binding Domain

1. Download 7c8j from PDB using VMD
2. `[atomselect top "chain B"] writepdb my_7c8j.pdb` (this file is provided here)
3. Follow Justin Lemkul's tutorial (steps below)

```bash
# generate the topology; use OPLS-AA force field (sel. 15)
$ gmx pdb2gmx -f my_7c8j.pdb -o my_7c8j_processed.pdb -water spce
# topol.top is also created
# enlarge box
$ gmx editconf -f my_7c8j_processed.pdb -o my_7c8j_newbox.gro -c -d 1.0 -bt cubic
# solvate
$ gmx solvate -cp my_7c8j_newbox.gro -cs spc216.gro -o my_7c8j_solv.gro -p topol.top
# copy JK's ions.mdp; add ions
$ gmx grompp -f ions.mdp -c my_7c8j_solv.gro -p topol.top -o ions.tpr
$ gmx genion -s ions.tpr -o my_7c8j_solv.gro -p topol.top -pname NA -nname CL -neutral
# selected group 14 SOL
# copy JK's minim.mdp
# minimize potential energy
$ gmx grompp -f minim.mdp -c my_7c8j_solv.gro -p topol.top -o em.tpr
$ gmx mdrun -v -deffnm em
# have a look at the potential energy
$ gmx energy -f em.edr -o potential.xvg
# plot using xvg_plot.py forked from Joao Rodrigues
# copy JK's nvt.mdp; run NVT MD for 50,000 steps -- this will take a few hours...
$ gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
$ gmx mdrun -v -deffnm nvt
```
