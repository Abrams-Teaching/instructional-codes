import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
import scipy.integrate as si

def lj(r,eps,sig):
    r6=r*r*r*r*r*r
    r6i=pow(sig,6)*np.reciprocal(r6)
    return 4*eps*(r6i*r6i-r6i)

parser=ap.ArgumentParser()
parser.add_argument('-i',type=str,default=[],action='append',help='input file name')
parser.add_argument('-o',type=str,default='rdf.png',help='output file name')
parser.add_argument('-rho',type=float,default=0.0,help='density (for integrations)')
parser.add_argument('-R',type=float,default=0.0,help='outer radius (for integrations)')
parser.add_argument('-eps',type=float,default=0.0,help='LJ epsilon (for integrations)')
parser.add_argument('-sig',type=float,default=0.0,help='LJ sigma (for integrations)')

args=parser.parse_args()

fig,ax=plt.subplots(1,1,figsize=(6,5))
ax.set_xlabel("r ($\sigma$)")
ax.set_ylabel("g(r)")
for f in args.i:
    r,g=np.loadtxt(f,unpack=True)
    if args.R>0.0:
        for i in range(len(r)):
            if r[i]>args.R:
                break
        ig=g[:i]*r[:i]*r[:i]
        n=4*np.pi*args.rho*si.simpson(ig,x=r[:i])
        print("n={:.3f}".format(n))
        if args.eps>0.0:
            ig=ig[1:]*lj(r[1:i],args.eps,args.sig)
            u=2*np.pi*args.rho*si.simpson(ig,x=r[1:i])
            print("u={:.3f}".format(u))
    ax.plot(r,g)
plt.savefig(args.o,bbox_inches='tight')
