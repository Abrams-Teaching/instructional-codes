'''
   Blocked averages computation a la Flyvbjerg and Petersen.
   [J. Chem. Phys., 91:461-466,1989]

   Reads in a Y vs. X data, and recursively block adjacent 
   elements, computing variance.

   Expects input in the form of progress lines from mdlj.c
   <timestep> <time> <PE> <KE> <TE> <drift> <T> <P>

   Blocks any one of columns 2 (<PE>), 3, 4, 5, 6, 7

   Outputs column-oriented data as
   <generation> <mean> <variance>

   Cameron F. Abrams

   Written for the course CHE T580, Modern Molecular Simulations
   Spring 20-21

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2021
'''
import numpy as np
import argparse as ap
import matplotlib.pyplot as plt

def block ( A ):
    Ab=[]
    for i in range(len(A)//2):
        Ab.append(0.5*(A[2*i]+A[2*i-1]))
    return np.array(Ab)

if __name__ == '__main__':
    parser=ap.ArgumentParser()
    parser.add_argument("-l",type=str,default="log",help="mdlj log file")
    parser.add_argument("-d",default=2,type=int,help="data column to analyze")
    parser.add_argument("-o",default='out.png',type=str,help="output file name")
    parser.add_argument("-N",default=1,type=int,help="number of particles in simulation")
    parser.add_argument("-show-column-labels",default=True,type=bool,help="show column labels")
    args=parser.parse_args()

    data=np.loadtxt(args.l)
    y=data[:,args.d]
    print(len(y),'# data points')
    c1n=[]
    c1s=[]
    j=0
    while (len(y)>4):
        mn=y.mean()/args.N
        s2=y.var()/(args.N*args.N)
        av=np.sqrt(s2/(len(y)-1))
        sav=av/np.sqrt(2*(len(y)-1))
        print('{:d} {:.5f} {:.5f} {:.5f}'.format(j,mn,av,sav))
        c1n.append(j)
        c1s.append(av)
        y=block(y)
        j+=1

    data=np.loadtxt(args.l)
    y=data[:149000,args.d]
    print(len(y),'# data points')
    c2n=[]
    c2s=[]
    j=0
    while (len(y)>4):
        mn=y.mean()/args.N
        s2=y.var()/(args.N*args.N)
        av=np.sqrt(s2/(len(y)-1))
        sav=av/np.sqrt(2*(len(y)-1))
        print('{:d} {:.5f} {:.5f} {:.5f}'.format(j,mn,av,sav))
        c2n.append(j)
        c2s.append(av)
        y=block(y)
        j+=1

    fig,ax=plt.subplots(1,1,figsize=(5,4))
    ax.set_xlabel('M')
    ax.set_ylabel('$\sigma$')
    ax.scatter(c2n,c2s,marker='^',label='150,000')
    ax.scatter(c1n,c1s,marker='s',label='600,000')
    plt.legend()
    plt.savefig(args.o,bbox_inches='tight')
