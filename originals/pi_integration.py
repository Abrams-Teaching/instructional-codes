from random import random
from math import pi, sqrt
import sys
import matplotlib.pyplot as plt

def f(x):
    return sqrt(1-x*x)

N = int(sys.argv[1])
M = int(sys.argv[2])
x=[]

I = [0.0]*M
base = 2
initbase = 2
Y=[]
for i in range(1,N+1):
    for j in range(M):
        I[j] += f(random())
    if i%initbase == 0:
        x.append(i)
        Y.append(sum([4*I[j]/(i+1) for j in range(M)])/M)
        initbase*=base

fig,ax=plt.subplots(2,1)
ax[0].plot(x,[pi]*len(x),'k--')

ax[0].scatter(x,Y,color='red')
ax[0].set_xscale('log')
ax[0].set_ylabel('$4I$')
ax[1].plot(x,[(v-pi)**2 for v in Y],color='red')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_ylabel('$(4I-\pi)^2$')
ax[1].set_xlabel('number of random draws')
plt.savefig('myplot.png')