# Run multiple simulations in parallel on a multiprocessor machine
# Cameron Abrams cfa22@drexel.edu

import os
from functools import partial
from multiprocessing.dummy import Pool
from subprocess import call
import argparse as ap
import numpy as np
import matplotlib.pyplot as plt

# Determine number of CPUs available
nPe_detected=os.cpu_count()

parser=ap.ArgumentParser()
parser.add_argument('-namd2',metavar='<str>',type=str,default='/home/cfa/namd/NAMD_2.14_Source/Linux-x86_64-g++/namd2',help='NAMD2 executable')
parser.add_argument('-base-config',metavar='<str>',type=str,default='base.namd',help='base config file')
parser.add_argument('-np',metavar='<int>',type=int,default=4,help='Number of processors to use in parallel')
parser.add_argument('-n',metavar='<int>',type=int,default=16,help='Number of windows')
parser.add_argument('-T',metavar='<float>',type=float,default=298.0,help='Temperature (K)')
parser.add_argument('-nsteps',metavar='<int>',type=int,default=5000000,help='Number of MD steps')
parser.add_argument('-k',metavar='<float>',type=float,default=100.0,help='Window potential spring constant')
parser.add_argument('-zlim',metavar=('zmin<float>','zmax<float>'),type=float,default=[2.5,4.5],nargs=2,help='limits on histogram domain')
parser.add_argument('-info',metavar='<str>',type=str,default='info.in',help='WHAM input file to create')
parser.add_argument('-dry-run',action='store_true')
args=parser.parse_args()
# Determine number requested (if any)
nPe=args.np
if nPe>nPe_detected:
    print('Warning: requested {:d} processors but only {:d} are found.'.format(nPe,nPe_detected))
    nPe=nPe_detected

# windows 
win_edges = np.linspace(args.zlim[0],args.zlim[1],args.n+1)
wz=win_edges[:-1]+0.5*(args.zlim[1]-args.zlim[0])/args.n
# Default values for program name and command-line arguments
prg=args.namd2
options={}
base_arguments=['+p1']
cfgs={'%T%':str(args.T),'%K%':str(args.k),'%NSTEPS%':str(args.nsteps)}
with open(args.base_config,"r") as f:
    base_config = f.read()

# build list of fully-resolved command names
with open(args.info,"w") as f:
    f.write('# i  k  z_0\n')
commands=[]
for i,w in enumerate(wz):
    with open(args.info,"a") as f:
        f.write('{:d} {:.4f} {:.5f}\n'.format(i,args.k,w))
    arguments=base_arguments[:]
    new_config=base_config[:]
    cfgs['%W%']=str(i)
    cfgs['%Z0%']=str(w)
    for k,v in cfgs.items():
        subst=new_config.replace(k,v)
        new_config=subst
    with open('{:d}.namd'.format(i),"w") as f:
        f.write(new_config)
    arguments.append('{:d}.namd'.format(i))
    log='log{:d}.out'.format(i)
    commands.append(prg+' '+' '.join([' '.join(['-'+k,str(v)]) for k,v in options.items()])+' '+' '.join(arguments)+' > {:s}'.format(log))

if not args.dry_run:
    # deploy simulations onto a pool of processor elements
    pool = Pool(nPe)
    for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
        if returncode != 0:
            print("%d command failed: %d" % (i, returncode))
