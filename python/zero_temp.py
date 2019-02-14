import subprocess as sp
import numpy as np
import time
import sys

L=int(sys.argv[1])
N=int(sys.argv[2])
T=int(sys.argv[3])
B=2*L
P=float(sys.argv[4])
NS=int(sys.argv[5])
NT=int(sys.argv[6])

lx=L
ly=L
nsweep=N
thermal=T
beta=B
J=1.823
dJ=0.5
P=P
seed=21873
n_thread=NT

loop = range(1,NS+1)

pipe=[]
count=n_thread
finish=0
i=0
while i<len(loop):
    if count>0:
        seed_t=seed+i
        args='./exe -n '+str(nsweep)+' -t '+str(thermal)+' -x '+str(lx)+' -y '+str(ly)+' -m 4 -s '+str(loop[i])+' -j '+str(J)+' -d '+str(dJ)+' -p '+str(P)+' -b '+str(beta)
        p = sp.Popen(args,shell=True)
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
