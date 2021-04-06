'''
runs the sample ising_mc C program and plots results
'''
import os
from random import randint
import matplotlib.pyplot as plt
import matplotlib.cm as cm

outfile='T{:d}-t{:d}.dat'
cmd='./ising_mc -L 40 -T {:.1f} -s {:d} -nc 50000 -fs 100'
Ts=[7,6,5,4,3,2,1]
trials = 6
e=[]
s1=[]
for i in range(trials):
    e.append([])
    s1.append([])
    for j,T in enumerate(Ts):
        seed = randint(99999,999999)
        this_cmd=cmd.format(T,seed)+' > '+outfile.format(j,i)
        print(outfile.format(j,i))
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


