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
parser.add_argument('-n',type=int,default=4,help='number of processors to use in parallel')
parser.add_argument('-plotonly',action='store_true',help='only make plot')
parser.add_argument('-o',type=str,default='plot.png',help='plot file name')
args=parser.parse_args()
# Determine number requested (if any)
nPe=args.n
# Don't allow number of CPU's to be exceeded
if nPe>nPe_detected:
    print('Warning: requested {:d} processors but only {:d} are found.'.format(nPe,nPe_detected))
    nPe=nPe_detected

# Specify the temperature-pressure parameter space
k = 16.0
wxs = np.linspace(-8,8,16)

# Default values for program name and command-line arguments
prg='bd-w'
options={'ns':50000000,'k-win':k,'hist-n':1000}

# build list of fully-resolved command names
if not args.plotonly:
    commands=[]
    for i,wx in enumerate(wxs):
        options['x-win']=wx
        options['which-win']=i
        options['plot-w']='win-pot{:d}.dat'.format(i)
        options['plot-f']='f-pot{:d}.dat'.format(i)
        log='log{:d}.out'.format(i)
        commands.append(r'./'+prg+' '+' '.join([' '.join(['-'+k,str(v)]) for k,v in options.items()])+' > {:s}'.format(log))

    # deploy simulations onto a pool of processor elements
    pool = Pool(nPe)
    for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
        if returncode != 0:
            print("%d command failed: %d" % (i, returncode))
