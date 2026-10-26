"""Large-M saddle for the thermal mass of the self-damped, randomly coupled 2D z=2 boson (same model as
bath-realizations/delta_selfconsistent.py), scanned over the bare phi^4 coupling u, the frequency cutoff Ec and T/E_L.
Units: damping gamma=1, boson bandwidth E_L=1, mode density rho=1. Writes ../data/universal_curves.json (used by make_figures_v6.py, drift_test_v6.py, v6_numbers.py).
Key output: r=m/T saturates as u->infinity (the constraint limit), so beta=S(0) has a coupling-independent MINIMUM
that depends on T/E_L only logarithmically."""
import numpy as np, os, warnings, json
HERE=os.path.dirname(os.path.abspath(__file__))
warnings.filterwarnings('ignore')
from scipy.integrate import quad
from scipy.optimize import brentq
import sys
sys.path.insert(0, HERE)
import bath   # same transport functional as the old shapes/cdf_shape.py (checked identical)
rho=1.0; EL=1.0
def Cb(r): return bath.Cb(bath.w_arc(r))
def S0(r): return bath.S0(bath.w_arc(r))
sig=lambda x: x/np.expm1(x)-np.log(-np.expm1(-x))
def s(T,m): return rho/(2*np.pi)*quad(lambda W: sig(W/T)*np.log(((m+EL)**2+W*W)/(m*m+W*W)),1e-14,60*T,limit=400,points=[m,T])[0]
B0=np.pi/3*rho
res={}
for Ec,u in [(1.0,1e3),(1.0,1.0),(1.0,0.3),(0.3,1e3),(3.0,1e3)]:
    def phi2_0(m): return rho/np.pi*quad(lambda e: np.log((m+e+Ec)/(m+e)),0,EL,limit=200,points=[max(m,1e-12)*10])[0]
    def phi2_T(T,m):
        h=lambda W: 2/np.expm1(W/T)*(np.arctan((m+EL)/W)-np.arctan(m/W))
        return rho/np.pi*(quad(h,1e-14,50*T,limit=400,points=[m,T])[0])
    P00=rho/np.pi*quad(lambda e: np.log((e+Ec)/e),0,EL,limit=200)[0]
    lnTs=np.linspace(np.log(1e-6),np.log(5e-2),29); Ts=np.exp(lnTs)
    ms=np.array([brentq(lambda m: m/u-(phi2_0(m)+phi2_T(T,m)-P00),1e-16,300*T) for T in Ts])
    Sv=np.array([s(T,m) for T,m in zip(Ts,ms)])
    CT=np.gradient(Sv,Ts); B=-np.gradient(CT,lnTs); d=B/B0-1
    r=ms/Ts; be=np.array([S0(x) for x in r]); cb=np.array([Cb(x) for x in r])
    dbeta=np.gradient(be,lnTs)
    key=f"Ec{Ec}_u{u:g}"
    res[key]=dict(T=Ts.tolist(),r=r.tolist(),beta=be.tolist(),Cb=cb.tolist(),delta=d.tolist(),dbeta_dlnT=dbeta.tolist())
    print(key)
    for i in range(2,27,3):
        print(f"  T/EL={Ts[i]:.1e} r={r[i]:.3f} beta={be[i]:.2f} dbeta/dln10T={dbeta[i]*np.log(10):+.2f} Cb={cb[i]:.2f} delta={d[i]:+.3f} Q={cb[i]/(1+d[i]):.2f}")
json.dump(res,open(os.path.join(HERE,'..','data','universal_curves.json'),'w'))
