import sys, numpy as np, json
sys.path.insert(0,''+__import__('os').path.dirname(__import__('os').path.abspath(__file__))+'')
from matsu import solve, free_energy
Uf=0.05; D=1.0
U=float(sys.argv[1]); dtau=float(sys.argv[2]); r=float(sys.argv[3]) if len(sys.argv)>3 else 0.0
Ts=np.exp(np.linspace(np.log(1e-4),np.log(4e-3),24))
out=[]
for T in Ts:
    s=solve(T,D=D,Uf=Uf,U=U,r=0.0,dtau=dtau,tol=3e-10,maxit=2000)
    out.append((T,free_energy(s),s['it'],s['err']))
    print(T,out[-1][1],s['it'],s['err'],flush=True)
np.save(f''+__import__('os').path.dirname(__import__('os').path.abspath(__file__))+f'/out/F_U{U}_dt{dtau}.npy',np.array(out))
