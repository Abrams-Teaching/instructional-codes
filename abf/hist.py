# read the series history files that abf generates in namd
# cameron f abrams cfa22@drexel.edu
#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import argparse as ap

def read_hist(fn):
    frm=0
    data={}
    with open(fn,'r') as f:
        lines=f.readlines()
        i=0
        while i < len(lines):
            #print(lines[i])
            i+=1
            hl=lines[i].split()
            #print(hl)
            minx=float(hl[1])
            dx=float(hl[2])
            nx=int(hl[3])
            i+=2
            tx=[]
            ty=[]
            for j in range(i,i+nx):
                tln=lines[j].split()
                #print(j,tln)
                tx.append(float(tln[0]))
                ty.append(float(tln[1]))
            i+=nx+1
            data[frm]=[np.array(tx),np.array(ty)]
            frm+=1
    return data

if __name__=='__main__':
    parser=ap.ArgumentParser()
    parser.add_argument('-f',type=str,default='',help='input')
    parser.add_argument('-o',type=str,default='plot.png',help='output')
    parser.add_argument('-cmap-name',type=str,default='viridis',help='colormap name')
    parser.add_argument('-figsize',type=float,nargs=2,default=[5,4])
    parser.add_argument('-xlabel',type=str,default='x')
    parser.add_argument('-ylabel',type=str,default='y')
    parser.add_argument('-ylim',type=float,nargs=2,default=[0,10])
    parser.add_argument('-xlim',type=float,nargs=2,default=[0,1])
    parser.add_argument('-lw',type=float,default=0.5)
    args=parser.parse_args()

    mydat=read_hist(args.f)
    print('Read {:d} datasets from {:s}.'.format(len(mydat),args.f))

    fig,ax=plt.subplots(1,1,figsize=args.figsize)
    cmap=cm.get_cmap(args.cmap_name)
    for frm,ds in mydat.items():
        ax.plot(ds[0],ds[1],color=cmap(frm/len(mydat)),linewidth=args.lw)
    
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(args.ylabel)
    ax.set_xlim(args.xlim)
    ax.set_ylim(args.ylim)
    plt.savefig(args.o,bbox_inches='tight')

