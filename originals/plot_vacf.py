import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
import scipy.integrate as si

parser=ap.ArgumentParser()
parser.add_argument('-i',type=str,default=[],action='append',help='input file name')
parser.add_argument('-o',type=str,default='msd.png',help='output file name')
parser.add_argument('-z',type=int,default=10,help='zoom factor')

args=parser.parse_args()

fig,ax=plt.subplots(1,2,figsize=(12,5))
ax[0].set_xlabel("t")
ax[0].set_ylabel("VACF(t)")
ax[1].set_xlabel("t")
ax[1].set_ylabel("VACF(t)")

for f in args.i:
    t,vacf=np.loadtxt(f,unpack=True)
    D=1./3.*si.simpson(vacf,x=t)
    print(f,D)
    ax[0].plot(t,vacf)
    ax[1].plot(t[:len(t)//args.z],vacf[:len(vacf)//args.z],label=f)
ax[1].legend()
plt.savefig(args.o,bbox_inches='tight')