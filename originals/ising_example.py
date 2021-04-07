'''
runs the sample ising_mc C program and plots results
Cameron Abrams cfa22@drexel.edu
'''
import os
from random import randint
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse as ap

parser=ap.ArgumentParser()
parser.add_argument("-x",default='./ising_mc',help='executable name')
parser.add_argument("-L",type=int,default=40,help='magnet side length')
parser.add_argument("-nc",type=int,default=50000,help='number of cycles')
parser.add_argument("-fs",type=int,default=100,help='sampling interval (cycles)')
parser.add_argument("-nt",type=int,default=6,help='number of trials per T-set')
parser.add_argument("-T",default='1,2,3,4,5,6,7',help='comma-separated temperatures')

args=parser.parse_args()

outfile='T{:d}-t{:d}.dat'
cmd='{:s} -L {:d} -T {:.1f} -s {:d} -nc {:d} -fs {:d}'
Ts=[float(x) for x in args.T.split(',')]
trials = args.nt
e=[]
s1=[]
for i in range(trials):
    e.append([])
    s1.append([])
    for j,T in enumerate(Ts):
        seed = randint(99999,999999)
        this_cmd=cmd.format(args.x,args.L,T,seed,args.nc,args.fs)+' > '+outfile.format(j,i)
        print('Generating',outfile.format(j,i))
        os.system(this_cmd)
        with open(outfile.format(j,i)) as f:
            for line in f:
                tokens=line.split()
                if tokens[0] == '#' and tokens[2] == 'average':
                    if tokens[3] == 'magnetization':
                        s1[i].append(float(tokens[-1]))
                    if tokens[3] == 'energy':
                        e[i].append(float(tokens[-1]))

fig,ax=plt.subplots(2,1)
ax[0].set_ylim([-1,1])
ax[1].set_ylim([-2,0])
ax[1].set_xlabel("T (J/k$_B$)")
ax[0].set_ylabel(r"$\left<s_i\right>$")
ax[1].set_ylabel(r"$\left<e\right>$")
cmap=cm.get_cmap('plasma_r')
for j,T in enumerate(Ts):
    ax[0].scatter([T]*trials,[s1[i][j] for i in range(trials)],color=cmap(j/len(Ts)))
    ax[1].scatter([T]*trials,[e[i][j] for i in range(trials)],color=cmap(j/len(Ts)))
plt.savefig('ising_ex.png')


