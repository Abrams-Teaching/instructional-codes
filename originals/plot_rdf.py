import matplotlib.pyplot as plt
import numpy as np
import argparse as ap

parser=ap.ArgumentParser()
parser.add_argument('-i',type=str,default=['rdf.dat'],action='append',help='input file name')
parser.add_argument('-o',type=str,default='rdf.png',help='output file name')
args=parser.parse_args()

fig,ax=plt.subplots(1,1,figsize=(6,5))
ax.set_xlabel("r ($\sigma$)")
ax.set_ylabel("g(r)")
for f in args.i:
    r,g=np.loadtxt(f,unpack=True)
    ax.plot(r,g)
plt.savefig(args.o,bbox_inches='tight')
