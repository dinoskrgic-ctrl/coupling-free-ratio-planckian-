import sys, numpy as np, json
sys.path.insert(0,''+__import__('os').path.dirname(__import__('os').path.abspath(__file__))+'')
from matsu import solve, free_energy
Uf=float(sys.argv[4]) if len(sys.argv)>4 else 0.05; D=1.0
U=float(sys.argv[1]); dtau=float(sys.argv[2]); r=float(sys.argv[3]) if len(sys.argv)>3 else 0.0
Ts=np.exp(np.linspace(np.log(1.5e-4),np.log(4e-3),14))
out=[]
for T in Ts:
    s=solve(T,D=D,Uf=Uf,U=U,r=r,dtau=dtau,tol=2e-9,maxit=(1500 if r else 800),mix=(0.3 if r else 0.5))
    if r!=0:
        s0=solve(T,D=D,Uf=Uf,U=0.0,r=0.0,dtau=dtau,tol=3e-10,maxit=2000)
        fe=free_energy(s,(s0['Sf'],s0['Gf'],s0['gf']))
    else:
        fe=free_energy(s)
    out.append((T,fe,s['it'],s['err']))
    print(T,out[-1][1],s['it'],s['err'],flush=True)
np.save(f''+__import__('os').path.dirname(__import__('os').path.abspath(__file__))+f'/out/F_U{U}_dt{dtau}_r{r}_Uf{Uf}.npy',np.array(out))
