import sys, numpy as np, time
sys.path.insert(0,''+__import__('os').path.dirname(__import__('os').path.abspath(__file__))+'')
from real import *
Uf=0.05; D=1.0; Us=[0.17,0.3,0.5]
Ts=[2.8e-3,2e-3,1.4e-3,1e-3,7e-4,5e-4,3.5e-4,2.5e-4,1.77e-4,1.25e-4]
g=RGrid(L=4.0,N=2**21)
A=None; res=[]
for T in Ts:
    t0=time.time()
    A,Sf,it,err=solve_f(g,T,Uf,A0=A,tol=1e-9)
    row=[T]
    for U in Us:
        Ac,Sc,itc,errc=solve_c(g,T,D,U,A)
        tt=tau_tr(g,T,D,Sc)
        # lambda from mass renormalization at w ~ 0:  -dReS/dw at w=0 (T-dependent log)
        i0=np.argmin(np.abs(g.w)); dRe=-(Sc.real[i0+5]-Sc.real[i0-5])/(g.w[i0+5]-g.w[i0-5])
        row+= [1/(tt*T), dRe, -Sc.imag[i0]/(np.pi*T)]
        np.save(f''+__import__('os').path.dirname(__import__('os').path.abspath(__file__))+f'/out/Sc_U{U}_T{T}.npy',Sc[np.abs(g.w)<0.3])
    res.append(row)
    print(' '.join(f'{x:.6g}' for x in row), f'f_it={it} err={err:.1e} t={time.time()-t0:.0f}s',flush=True)
    np.save(''+__import__('os').path.dirname(__import__('os').path.abspath(__file__))+f'/out/real_alpha.npy',np.array(res))
np.save(''+__import__('os').path.dirname(__import__('os').path.abspath(__file__))+f'/out/w_small.npy',g.w[np.abs(g.w)<0.3])
