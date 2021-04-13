'''
runs the sample hdisk C program and plots results
'''
import os
from random import randint
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

outfile='dR{:d}-dW{:d}.dat'
cmd='./hddb -N 200 -rho 0.5 -dr {:.1f} -seed {:d} -nc 2000 -dw {:.1f}'

dws=[0.2,0.5,0.8]
drs=np.concatenate((np.linspace(0,2.5,22),np.linspace(3,10,16)))

ar=[]
for i in range(len(dws)):
    ar.append([])
    for j,dr in enumerate(drs):
        seed = randint(99999,999999)
        this_cmd=cmd.format(dr,seed,dws[i])+' > '+outfile.format(j,i)
        os.system(this_cmd)
        with open(outfile.format(j,i)) as f:
            for line in f:
                tokens=line.split()
                if tokens[0]=='Displacement' and tokens[1]=='acceptance':
                    ar[i].append(float(tokens[-1]))
fig,ax=plt.subplots(1,1)
ax.set_ylim([0,1])
ax.set_xlim([0,10])
ax.set_ylabel("Acceptance Ratio")
ax.set_xlabel(r"$\Delta R$")
cmap=cm.get_cmap('plasma_r')
for i,dw in enumerate(dws):
    ax.scatter(drs,ar[i],color=cmap((i+1)/(len(dws)+2)),label=r'disp. wt. = {:.1f}'.format(dw))
ax.legend()
plt.savefig('hdiskdumbbell_ex.png')