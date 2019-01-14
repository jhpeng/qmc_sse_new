import subprocess as sp

lx=4
ly=4
beta=5
J=1.9
seed=21873
pipe=[]

loop = [1.8,1.9,2.0,2.1]

for J in loop:
    seed=seed+1
    args='/Users/jhpeng/Works/qmc/qmc_sse_new/exe -x '+str(lx)+' -y '+str(ly)+' -m 1 -b '+str(beta)+' -j '+str(J)+' -s '+str(seed)
    p = sp.Popen(args,shell=True)
    pipe.append(p)

for i in range(len(loop)):
    pipe[i].wait()
