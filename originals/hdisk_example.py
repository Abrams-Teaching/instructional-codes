'''
runs the sample hdisk C program and plots results
'''
import os
from random import randint
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

outfile='dR{:d}-r{:d}.dat'
cmd='./hdisk -N 200 -rho {:.1f} -dr {:.1f} -seed {:d} -nc 2000'

rhos=[0.2,0.4,0.6]
drs=np.concatenate((np.linspace(0,2.5,22),np.linspace(3,10,16)))

ar=[]
for i in range(len(rhos)):
    ar.append([])
    for j,dr in enumerate(drs):
        seed = randint(99999,999999)
        this_cmd=cmd.format(rhos[i],dr,seed)+' > '+outfile.format(j,i)
        os.system(this_cmd)
        with open(outfile.format(j,i)) as f:
            for line in f:
                tokens=line.split()
                if tokens[0]=='Acceptance':
                    ar[i].append(float(tokens[-1]))
fig,ax=plt.subplots(1,1)
ax.set_ylim([0,1])
ax.set_xlim([0,10])
ax.set_ylabel("Acceptance Ratio")
ax.set_xlabel(r"$\Delta R$")
cmap=cm.get_cmap('plasma_r')
for i,r in enumerate(rhos):
    ax.scatter(drs,ar[i],color=cmap((i+1)/(len(rhos)+2)),label=r'$\rho$ = {:.1f}'.format(r))
ax.legend()
plt.savefig('hdisk_ex.png')