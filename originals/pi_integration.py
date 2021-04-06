'''
Compute pi using Monte Carlo integration with uniform variants on the 1-D domain [0,1]

Cameron F Abrams cfa22@drexel.edu

Run using
/usr/bin/python3 /home/cfa/cheT580-202035/instructional-codes/originals/pi_integration.py 1000000 100
where first argument is the number of samples and the second is the number of trials

generates a logarithmically scaled set of two plots
'''
from random import random
from math import pi, sqrt
import sys
import matplotlib.pyplot as plt
import scipy.integrate as si
import numpy as np

def f(x):
    return np.sqrt(1-x*x)

N = int(sys.argv[1]) # number of trials per cycle
M = int(sys.argv[2]) # number of cycles
x = []

I = [0.0]*M
srI = []
base = 2
initbase = 2
Y=[]
for i in range(1,N+1):
    for j in range(M):
        I[j] += f(random())
    if i%initbase == 0:
        x.append(i*M)
        Y.append(sum([4*I[j]/(i+1) for j in range(M)])/M)
        srx=np.linspace(0,1,i*M)
        sry=f(srx)
        srI.append(4*si.simpson(sry,srx))
        initbase*=base

fig,ax=plt.subplots(2,1)
ax[0].plot(x,[pi]*len(x),'k--',label=None)
ax[0].scatter(x,Y,color='red',label='MC')
ax[0].scatter(x,srI,color='blue',label='Simpson\'s')
ax[0].set_xscale('log')
ax[0].set_ylabel('$4I$')
ax[0].legend()
ax[1].plot(x,[(v-pi)**2 for v in Y],color='red',label='MC')
ax[1].plot(x,[(v-pi)**2 for v in srI],color='blue',label='Simpson\'s')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_ylabel('$(4I-\pi)^2$')
ax[1].set_xlabel('number of points on [0,1]')
plt.savefig('myplot.png')