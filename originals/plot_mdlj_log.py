# Plots columns generated by mdlj.c
# Cameron F Abrams cfa22@drexel.edu

import numpy as np
import matplotlib.pyplot as plt
import argparse as ap
import statsmodels.api as sm
lowess = sm.nonparametric.lowess

def get_data ( logs, cols, col_labels ):
    alldat={}
    pdat={}
    plot_labels={}
    for log in args.logs:
        alldat[log]=np.loadtxt(log)
        x=alldat[log][:,1]
        y=[]
        this_labels=[]
        for c in cols:
            y.append(alldat[log][:,c].copy())
            if len(col_labels)>0:
                this_labels.append(col_labels[c])
            else:
                this_labels.append('')
        pdat[log]=[x,y]
        plot_labels[log]=this_labels
    return pdat,plot_labels

def get_NAMD_data ( logs, cols, col_labels ):
    alldat={}
    pdat={}
    plot_labels={}
    for log in args.logs:
        linedat=[]
        with open(log) as f:
            for l in f:
                tok=l.strip().split()
                if len(tok)>0 and tok[0]=='ENERGY:':
                    linedat.append([int(tok[1])] + [float(_) for _ in tok[2:]])
        alldat[log]=np.array(linedat).T
        print(alldat[log])
#        print(alldat[log][0,:])
        exit()
        x=alldat[log][:,0]
        y=[]
        this_labels=[]
        for c in cols:
            y.append(alldat[log][:,c].copy())
            if len(col_labels)>0:
                this_labels.append(col_labels[c])
            else:
                this_labels.append('')
        pdat[log]=[x,y]
        plot_labels[log]=this_labels
    return pdat,plot_labels

def block ( A ):
    Ab=[]
    for i in range(len(A)//2):
        Ab.append(0.5*(A[2*i]+A[2*i-1]))
    return np.array(Ab)

def flyberg ( y, minblocks=4):
    c1n=[]
    c1m=[]
    c1s=[]
    j=0
    while (len(y)>minblocks):
        mn=y.mean()
        s2=y.var()
        av=np.sqrt(s2/(len(y)-1))
        sav=av/np.sqrt(2*(len(y)-1))
        print('{:d} {:.5f} {:.5f} {:.5f}'.format(j,mn,av,sav))
        c1n.append(j)
        c1m.append(mn)
        c1s.append(av)
        y=block(y)
        j+=1
    return c1n,c1m,c1s

parser=ap.ArgumentParser()
parser.add_argument("-logs",default=[],action='append',type=str,help="log file")
parser.add_argument("-fmt",default='',type=str,help="format * or NAMD")
parser.add_argument("-labels",default=[],action='append',type=str,help="labels")
parser.add_argument("-d",default='2',type=str,help="comma-separated-list of data columns to plot vs time")
parser.add_argument("-ylabel",default='energy ($\epsilon$)',type=str,help="y-axis label")
parser.add_argument("-o",default='out.png',type=str,help="output file name")
parser.add_argument("-N",default=1,type=int,help="number of particles")
parser.add_argument("-show-column-labels",default=False,action='store_true',help="show column labels")
parser.add_argument("-xlog",default=False,action='store_true',help="make x axis log scale")
parser.add_argument("-do-flyvberg",default=False,action='store_true',help="Do Flyvberg analysis")
parser.add_argument("-divide-by-N",default=False,action='store_true',help="what do you think?")
parser.add_argument("-fluc-in-leg",default=False,action='store_true',help="what do you think?")
parser.add_argument("-every",default=1,type=int,help="only plot every # points")
args=parser.parse_args()

cols=list(map(int,args.d.split(',')))

if len(args.logs)==0:
    args.logs.append("log")
col_labels=[]
with open(args.logs[0]) as f:
    for l in f:
        tok=l.split()
        #print(len(tok),tok)
        if len(tok)>0 and (tok[0]=="#LABELS" or tok[0]=='ETITLE:'):
            col_labels=tok[1:]
            break
if args.show_column_labels:
    for i,c in enumerate(col_labels):
        print(i,c)
    exit()

if args.fmt=='NAMD':
    pdat,plot_labels = get_NAMD_data(args.logs,cols,col_labels)
else:
    pdat,plot_labels = get_data(args.logs,cols,col_labels)

fig,ax=plt.subplots(1,1,figsize=(6,4))
if args.xlog:
    ax.set_xscale('log')
if len(args.labels)==len(args.logs):
    for l,label in zip(args.logs,args.labels):
        plot_labels[l]=[label]
for log in args.logs:
    x,yy=pdat[log]
    ll=plot_labels[log]
    for y,label in zip(yy,ll):
        relfluc=y[len(y)//2:].var()/(y[len(y)//2:].mean()**2)*args.N
#        print('{:.3f}'.format(relfluc))
        if args.divide_by_N:
            y/=args.N
        raw_alpha=1
        if args.every>1:
            lowessfrac=args.every/len(y)
            smres=lowess(y,x,frac=lowessfrac)
            plot_x=smres[:,0]
            plot_y=smres[:,1]
            raw_alpha=0.5
        if args.fluc_in_leg:
            label+=', {:.3f}'.format(relfluc)
        if args.every>1:
            ax.plot(plot_x,plot_y,alpha=1,label=r'{:s}'.format(label))
            ax.plot(x,y,alpha=raw_alpha)
        else:
            ax.plot(x,y,alpha=raw_alpha,label=r'{:s}'.format(label))
        if args.do_flyvberg:
            n,m,s=flyberg(y)

if len(col_labels)>0:
    ax.set_xlabel(col_labels[1])
ax.set_ylabel(args.ylabel)
ax.legend()
plt.savefig(args.o,bbox_inches='tight')

