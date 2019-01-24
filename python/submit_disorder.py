import subprocess as sp
import numpy as np
import time

lx=64
ly=64
nsweep=8000
thermal=5000
beta_i=1.0
beta_f=3.2
interv=0.1
J=1.990
dJ=0.5
P=0.52
seed=21873
n_thread=4


loop = range(1,4+1)

pipe=[]
count=n_thread
finish=0
i=0
while i<len(loop):
    if count>0:
        seed_t=seed+i
        args='./exe -n '+str(nsweep)+' -t '+str(thermal)+' -x '+str(lx)+' -y '+str(ly)+' -m 3 -s '+str(loop[i])+' -j '+str(J)+' -i '+str(beta_i)+' -f '+str(beta_f)+' -v '+str(interv)+' -d '+str(dJ)+' -p '+str(P)
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
