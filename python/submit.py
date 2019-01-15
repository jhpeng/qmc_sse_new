import subprocess as sp
import time

lx=4
ly=4
beta=4
J=1.9
seed=21873
n_thread=3

loop = [1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2]

pipe=[]
count=n_thread
finish=0
i=0
while i<len(loop):
    if count>0:
        print(i)
        seed_t=seed+i
        args='/Users/jhpeng/Works/qmc/qmc_sse_new/exe -x '+str(lx)+' -y '+str(ly)+' -m 1 -b '+str(beta)+' -j '+str(loop[i])+' -s '+str(seed_t)
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
