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
parser.add_argument('-np',type=int,default=4,help='number of processors to use in parallel')
parser.add_argument('-o',type=str,default='plot.png',help='plot file name')
parser.add_argument('-n',type=int,default=16,help='number of windows')
parser.add_argument('-T',type=float,default=1.0,help='Temperature')
parser.add_argument('-k',type=float,default=10.0,help='Window potential spring constant')
parser.add_argument('-zlim',type=float,default=[-8,8],nargs='+',help='limits on histogram domain')
args=parser.parse_args()
# Determine number requested (if any)
nPe=args.np
# Don't allow number of CPU's to be exceeded
if nPe>nPe_detected:
    print('Warning: requested {:d} processors but only {:d} are found.'.format(nPe,nPe_detected))
    nPe=nPe_detected

# windows 
win_edges = np.linspace(args.zlim[0],args.zlim[1],args.n+1)
wz=win_edges[:-1]+0.5*(args.zlim[1]-args.zlim[0])/args.n

# Default values for program name and command-line arguments
prg='bd-w'
options={'ns':50000000,'k-win':args.k,'hist-n':1000,'T':args.T}

# build list of fully-resolved command names
commands=[]
for i,w in enumerate(wz):
    options['x-win']=w
    options['which-win']=i
    options['plot-w']='win-pot{:d}.dat'.format(i)
    if i==0:
        options['plot-f']='f-pot.dat'
    log='log{:d}.out'.format(i)
    commands.append(r'./'+prg+' '+' '.join([' '.join(['-'+k,str(v)]) for k,v in options.items()])+' > {:s}'.format(log))

# deploy simulations onto a pool of processor elements
pool = Pool(nPe)
for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
    if returncode != 0:
        print("%d command failed: %d" % (i, returncode))
