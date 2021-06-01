# GROMACS

![Gromacs](http://www.gromacs.org) is a popular molecular-dynamics package.  Most Linux distrubtions offer Gromacs as precompiled packages that can be installed; e.g., for Ubuntu (even inside WSL):

```bash
$ sudo apt install gromacs
```

In most high-performance computing settings, Gromacs is compile from source code in order to link in hardware-specific libraries for things like communication and base arithmetic.  For this introductory survey though, we can just run the precompiled version on a laptop.  (If you have macOS, you might have to compile Gromacs from source.)

Gromacs has been actively maintained since 1991 and several new releases appear each year with performance improvements and feature enhancements.  As it stands now, Gromacs is highly optimized to simulate systems of large biomacromolecules in biologically relevant settings, like water solvent or in lipid bilayer membranes, but the community of users continuously demostrates how to use it in a much broader variety of settings.  It is simply a good general-purpose MD engine.

The ![definitive Gromacs tutorials](mdtutorials.com) are maintained by Justin Lemkul and are a great place to start.  For this course, we'll consider very simple systems we can build ourselves.

## A Box of Water

For example, let's just consider a box of water.  First, let's use the predefined single-point-charge (SPC) water topology that comes with gromacs already to fill a small 3x3x3 nm<sup>3</sup> box:

```bash
$ gmx insert-molecules -ci /usr/share/gromacs/top/spc216.gro -nmol 100 -box 3 3 3 -o water_box.gro
```

This gets me 432 water molecules in this box (your results may vary since the boxes are packed randomly).  Now, we need to generate the topology.  I'm going to select the CHARMM force field even though it doesn't matter for a system of pure water.  Next, I'll stipulate that I'm using the TIP3P water model.  BOth of these selections are made interactively while running the `pdb2gmx` command:

```bash
$ gmx pdb2gmx -f water_box.gro -o new_water_box.gro -p topol.top
```
With that, we can now generate a `tpr` file to run a minimization.  The required parameter file is `minim.mdp`, provided here.  The `tpr` is built using the `grompp`, followed by `mdrun` to run the minimization:
```
$ gmx grompp -f minim.mdp -c new_water_box.gro -p topol.top -o min.tpr
$ gmx mdrun -v -deffnm min
```
With the system minimized, we can feed the resulting coordinates, in `min.gro`, along with `nvt.mdp` parameter file, to `gmx grompp` to build `nvt.tpr` for an NVT MD simulation:
```bash
$ gmx grompp -f nvt.mdp -c min.gro -p topol.top -o nvt.tpr
$ gmx mdrun -v -deffnm nvt
```

Once that MD simulation finishes, you can try an NPT simulation using the `npt.mdp` parameter file:
```bash
$ gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
$ gmx mdrun -v -deffnm npt
```
Energy-like output appears in `npt.edr` but in a binary format.  To extract ASCII data, use
```bash
$ gmx energy -f npt.edr
```
