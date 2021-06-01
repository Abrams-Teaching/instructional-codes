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
parser.add_argument('-np',metavar='<int>',type=int,default=4,help='number of processors to use in parallel')
parser.add_argument('-n',metavar='<int>',type=int,default=16,help='number of windows')
parser.add_argument('-T',metavar='<float>',type=float,default=1.0,help='Temperature')
parser.add_argument('-nsteps',metavar='<int>',type=int,default=2000000,help='Number of BD steps')
parser.add_argument('-hist-n',metavar='<int>',type=int,default=1000,help='Number of z-histogram bins')
parser.add_argument('-k',metavar='<float>',type=float,default=15.0,help='Window potential spring constant')
parser.add_argument('-zlim',metavar=('zmin<float>','zmax<float>'),type=float,default=[-8,8],nargs=2,help='limits on histogram domain')
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
prg='./bd-w'
options={'ns':args.nsteps,'k-win':args.k,'hist-n':args.hist_n,'T':args.T}

# build list of fully-resolved command names
commands=[]
for i,w in enumerate(wz):
    options['x-win']=w
    options['which-win']=i
    options['plot-w']='W{:d}.dat'.format(i)
    if i==0:
        options['plot-f']='V.dat'
    else:
        if 'plot-f' in options:
            del options['plot-f']
    #print(options)
    log='log{:d}.out'.format(i)
    commands.append(prg+' '+' '.join([' '.join(['-'+k,str(v)]) for k,v in options.items()])+' > {:s}'.format(log))

# deploy simulations onto a pool of processor elements
pool = Pool(nPe)
for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
    if returncode != 0:
        print("%d command failed: %d" % (i, returncode))
