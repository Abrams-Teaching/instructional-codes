'''
runs the sample mclj C program and plots results
'''
import os
from random import randint
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

outfile='rho{:d}-T{:d}.dat'
cmd='./mclj -N 512 -rho {:.1f} -s {:d} -nc 500 -ne 500 -rc 2.5 -T {:.1f} -dr 0.5 -sh'

Ts=[0.9,2.0]
rhos=np.linspace(0.1,0.9,9)

P=[]
for i in range(len(Ts)):
    P.append([])
    for j,r in enumerate(rhos):
        seed = randint(99999,999999)
        this_cmd=cmd.format(r,seed,Ts[i])+' > '+outfile.format(j,i)
        os.system(this_cmd)
        print(this_cmd)
        with open(outfile.format(j,i)) as f:
            for line in f:
                tokens=line.split()
                if tokens[0]=='Total' and tokens[1]=='Pressure:':
                    P[i].append(float(tokens[-1]))
fig,ax=plt.subplots(1,1)
ax.set_ylim([-1,10])
ax.set_xlim([0,1])
ax.set_ylabel("Pressure ($\epsilon/\sigma^3$)")
ax.set_xlabel("Density ($\sigma^{-3}$)")
cmap=cm.get_cmap('plasma_r')
for i,T in enumerate(Ts):
    ax.scatter(rhos,P[i],color=cmap((i+1)/(len(Ts)+2)),label=r'$T$ = {:.1f}'.format(T))
ax.legend()
plt.savefig('mclj_ex_2.png')