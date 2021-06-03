import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
import matplotlib.cm as cm

cmap=cm.get_cmap('inferno')

parser=ap.ArgumentParser()
parser.add_argument('-nrep',type=int,default=0)
parser.add_argument('-o',type=str,default='plot.png')
parser.add_argument('-figsize',type=float,nargs=2,default=[9,6])
args=parser.parse_args()

fig,ax=plt.subplots(1,1,figsize=args.figsize)
for i in range(args.nrep):
    with open('p{:d}.dat'.format(i)) as f:
        ln=f.readline()
        tok=ln.split()
        T=float(tok[3])
    xl,xr,h=np.loadtxt('p{:d}.dat'.format(i),unpack=True)
    ax.plot(0.5*(xl+xr),h,color=cmap(i/args.nrep),label='T = {:.3f}'.format(T))
ax.set_xlabel('$x$')
ax.set_ylabel('$H_i(x)$')
ax.legend()
plt.savefig(args.o,bbox_inches='tight')

t,x,a=np.loadtxt("rep0.log",unpack=True,)
fig,ax=plt.subplots(1,1,figsize=args.figsize)
ax.plot(t,x)
plt.savefig('trace.png')