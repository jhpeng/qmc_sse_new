import subprocess as sp
import numpy as np
import time

lx=16
ly=16
nsweep=4000
thermal=2000
beta_i=1.0
beta_f=32.0
interv=0.5
J=1.9
dJ=0.5
p=0.8
seed=21873
n_thread=3


loop = range(0,10)

pipe=[]
count=n_thread
finish=0
i=0
while i<len(loop):
    if count>0:
        seed_t=seed+i
        args='./exe -n '+str(nsweep)+' -t '+str(thermal)+' -x '+str(lx)+' -y '+str(ly)+' -m 2 -s '+str(loop[i])+' -j '+str(J)+' -i '+str(beta_i)+' -f '+str(beta_f)+' -v '+str(interv)+' -d '+str(dJ)+' -p '+str(p)
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
