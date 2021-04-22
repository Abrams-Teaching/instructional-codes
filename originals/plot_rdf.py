import matplotlib.pyplot as plt
import numpy as np
import argparse as ap

parser=ap.ArgumentParser()
parser.add_argument('-i',type=str,default='rdf.dat',help='input file name')
parser.add_argument('-o',type=str,default='rdf.png',help='output file name')
args=parser.parse_args()

r,g=np.loadtxt(args.i,unpack=True)
fig,ax=plt.subplots(1,1,figsize=(4,3))
ax.set_xlabel("r ($\sigma$)")
ax.set_ylabel("g(r)")
ax.plot(r,g,'b-')
plt.savefig(args.o,bbox_inches='tight')
