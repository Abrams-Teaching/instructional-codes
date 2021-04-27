import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
from scipy.optimize import curve_fit

def msd_D (t, D):
    return 6*D*t

parser=ap.ArgumentParser()
parser.add_argument('-i',type=str,default=[],action='append',help='input file name')
parser.add_argument('-o',type=str,default='msd.png',help='output file name')
parser.add_argument('-lowt',type=float,default=0.0,help='lower limit of time used for fitting')

args=parser.parse_args()

fig,ax=plt.subplots(1,2,figsize=(12,5))
ax[0].set_xlabel("t")
ax[0].set_ylabel("MSD(t)")
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[1].set_xlabel("1/t")
ax[1].set_ylabel("MSD(t)/(6t)")
ax[1].set_xscale('log')
ax[1].set_yscale('log')
for f in args.i:
    t,msd=np.loadtxt(f,unpack=True)
    ax[0].plot(t,msd,label='MD')
    ax[1].plot(np.reciprocal(t[1:]),msd[1:]/6*np.reciprocal(t[1:]),label='MD')
    if args.lowt>0.0:
        for i in range(len(t)):
            if t[i]>args.lowt:
                break
        par, cov = curve_fit(msd_D,t[i:],msd[i:],p0=[0.01])
        ax[0].plot(t,msd_D(t,*par),label='fit 6Dt; D={:.2e}'.format(par[0]),color='magenta')
        ax[1].plot(np.reciprocal(t[1:]),[par[0]]*len(t[1:]),label='fit; D={:.2e}'.format(par[0]),color='magenta')
ax[0].legend()
ax[1].legend()
plt.savefig(args.o,bbox_inches='tight')