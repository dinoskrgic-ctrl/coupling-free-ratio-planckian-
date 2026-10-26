import sys, numpy as np, time
sys.path.insert(0,''+__import__('os').path.dirname(__import__('os').path.abspath(__file__))+'')
from real import *
mode=sys.argv[1]
D=1.0; out=''+__import__('os').path.dirname(__import__('os').path.abspath(__file__))+f'/out/'
Ts=[2.8e-3,2e-3,1.4e-3,1e-3,7e-4,5e-4,3.5e-4,2.5e-4,1.77e-4,1.25e-4]
g=RGrid(L=4.0,N=2**21); A=None; res=[]
for T in Ts:
    t0=time.time()
    if mode=='uf0.1':
        A,Sf,it,err=solve_f(g,T,0.1,A0=A); Ac,Sc,itc,errc=solve_c(g,T,D,0.3,A)
    else:
        A,Ac,Sf,Sc,it,err=solve_fc(g,T,D,0.3,0.05,1.0,Af0=A)
    tt=tau_tr(g,T,D,Sc); res.append([T,1/(tt*T)])
    print(T,1/(tt*T),it,err,f'{time.time()-t0:.0f}s',flush=True)
    np.save(out+f'real_alpha_{mode}.npy',np.array(res))
