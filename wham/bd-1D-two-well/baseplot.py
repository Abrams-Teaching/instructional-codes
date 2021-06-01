import matplotlib.pyplot as plt
import numpy as np
import argparse as ap

parser=ap.ArgumentParser()
parser.add_argument('-T',type=float,default=1.0)
parser.add_argument('-f',type=str,default='p99.dat')
parser.add_argument('-V',type=str,default='V.dat')
parser.add_argument('-ylim',type=float,nargs=2,default=[0,25])
parser.add_argument('-o',type=str,default='baseplot.png')
args=parser.parse_args()
l,r,p=np.loadtxt(args.f,unpack=True)
z,v,dvdz=np.loadtxt(args.V,unpack=True)

fig,ax=plt.subplots(1,1,figsize=(6,5))
F=-args.T*np.log(p)
ax.set_ylim(args.ylim)
ax.plot((l+r)/2,F-F.min(),label='$F(x)=-kT\ln p_0$')
ax.plot(z,v-v.min(),label='V(x)',alpha=0.5)
ax.legend()
ax.set_label('$x$')
ax.set_ylabel('$F(x)$, $V(x)$')
plt.savefig(args.o,bbox_inches='tight')
