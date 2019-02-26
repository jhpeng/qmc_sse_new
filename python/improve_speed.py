import subprocess as sp
import numpy as np
import time
import sys

L=32

lx=L
ly=L
nblock=50
J=1.990
dJ=0.5
seed=21873
n_thread=3

it = 30
beta=100.0

pipe=[]
count=n_thread
finish=0
i=0
while i<it:
    if count>0:
        seed_t=seed+i
        args=[]
        args.append('./exe')
        args.append('-x')
        args.append(str(lx))
        args.append('-y')
        args.append(str(ly))
        args.append('-m')
        args.append('5')
        args.append('-j')
        args.append(str(J))
        args.append('-b')
        args.append(str(beta))
        args.append('-k')
        args.append(str(nblock))
        args.append('-s')
        args.append(str(seed_t))
        p = sp.Popen(args)
        pipe.append(p)
        count=count-1
        i=i+1
        beta = beta-2

    else:
        finish=0
        for j in range(i):
            if(pipe[j].poll()==0):
                finish=finish+1

        count = finish+n_thread-i
        time.sleep(0.5)
        
        
for i in range(len(loop)):
    pipe[i].wait()
