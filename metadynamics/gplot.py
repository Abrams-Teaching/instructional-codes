import matplotlib.pyplot as plt
import argparse as ap
import numpy as np
import matplotlib.cm as cm
import warnings

parser=ap.ArgumentParser()
parser.add_argument('-i',metavar='<str>',default=[],nargs='+',type=str,help='pmf files')
parser.add_argument('-l',metavar='<str>',default=[],nargs='+',type=str,help='labels')
parser.add_argument('-o',metavar='<str>',default='plot.png',type=str,help='output image file')
parser.add_argument('-ylim',default=[0,10],nargs=2,type=float,help='plot ylimits')
parser.add_argument('-xlim',default=[2.5,4.5],nargs=2,type=float,help='plot xlimits')
parser.add_argument('-figsize',type=float,nargs=2,default=[7,6])
parser.add_argument('-fontsize',type=float,default=14)
args=parser.parse_args()
if len(args.l)==0:
    args.l=args.i 

plt,ax=plt.subplots(1,1,figsize=args.figsize)
ax.set_ylim(args.ylim)
ax.set_xlim(args.xlim)
ax.set_xlabel('$z$ (Angstrom)',fontsize=args.fontsize)
ax.set_ylabel('$F(z)$ (kcal/mol)',fontsize=args.fontsize)
if len(args.i)>0:
    cmap=cm.get_cmap('inferno')
    pmf={}
    for i,(f,l) in enumerate(zip(args.i,args.l)):
        z,pmf[i]=np.loadtxt(f,unpack=True)
        ax.plot(z,pmf[i],color=cmap(i/len(args.i)),label=l)
ax.legend()

plt.savefig(args.o,bbox_inches='tight')
