import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
from scipy.optimize import curve_fit
from scipy.special import gamma

def my_gamma(x,beta,N):
    return beta**(N/2)/gamma(N/2)*x**(N/2-1)*np.exp(-beta*x)

parser=ap.ArgumentParser()
parser.add_argument('-i',type=str,default=[],action='append',help='input file name(s)')
parser.add_argument('-l',type=str,default=[],action='append',help='labels')
parser.add_argument('-o',type=str,default='edf.png',help='output file name')


args=parser.parse_args()

fig,ax=plt.subplots(1,1,figsize=(6,5))
ax.set_xlabel(r"$e$ ($\epsilon$)")
ax.set_ylabel(r"$f(e)$ ($\epsilon^{-1}$)")
ax.set_yscale('log')
for f,l in zip(args.i,args.l):
    e,d=np.loadtxt(f,unpack=True)
    #par,cov=curve_fit(my_gamma,e,d,p0=[1.0,200.0])
    if len(l)==0:
        label=f
    else:
        label=l
    ax.scatter(e,d,label=label)
    #ax.plot(e,my_gamma(e,*par),label='{:s}-fit({:.2f})'.format(label,par[1]))
ax.legend()
plt.savefig(args.o,bbox_inches='tight')
