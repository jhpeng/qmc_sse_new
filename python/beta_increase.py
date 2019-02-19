import subprocess as sp
import numpy as np
import time
import sys

L=int(sys.argv[1])
BI=float(sys.argv[2])
BF=float(sys.argv[3])
BV=float(sys.argv[4])
N=int(sys.argv[5])
T=int(sys.argv[6])
P=float(sys.argv[7])
NS=int(sys.argv[8])
NT=int(sys.argv[9])
SS=int(sys.argv[10])

lx=L
ly=L
nsweep=N
thermal=T
beta_i=BI
beta_f=BF
beta_v=BV
J=1.823
dJ=0.5
seed=SS
n_thread=NT

loop = range(1,NS+1)

pipe=[]
count=n_thread
finish=0
i=0
while i<len(loop):
    if count>0:
        seed_t=seed+i
        args=[]
        args.append('./exe')
        args.append('-x')
        args.append(str(lx))
        args.append('-y')
        args.append(str(ly))
        args.append('-m')
        args.append('2')
        args.append('-j')
        args.append(str(J))
        args.append('-i')
        args.append(str(beta_i))
        args.append('-f')
        args.append(str(beta_f))
        args.append('-v')
        args.append(str(beta_v))
        args.append('-p')
        args.append(str(P))
        args.append('-t')
        args.append(str(thermal))
        args.append('-n')
        args.append(str(nsweep))
        args.append('-s')
        args.append(str(seed_t))
        p = sp.Popen(args)
        pipe.append(p)
        count=count-1
        i=i+1

    else:
        finish=0
        for j in range(i):
            if(pipe[j].poll()==0):
                finish=finish+1

        count = finish+n_thread-i
        time.sleep(0.5)
        
        
for i in range(len(loop)):
    pipe[i].wait()
