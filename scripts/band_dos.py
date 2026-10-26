"""Band systematic (SM Sec. S6). b = B/(gamma_b*lambda) for LSCO p=0.24 from Michon et al. Nat Commun 14, 3033 (2023), arXiv:2205.04030.
Michon match the ln T slope of the optical mass (m*/m = 1 + 2g ln, so lambda = 2g) to the C/T slope of Girod et al.
ASSUMING C/T = gamma_b * m*_qp/m (i.e. b = 1). That fixes (m/m_e)K = 583 meV. Their tight-binding gives 582 meV.
If instead B = b*gamma_b*lambda, the data require (m/m_e)K_band = 583/b meV, so b = 583 meV / (m/m_e)K_TB.
(m/m_e)K is independent of the overall hopping t (K ~ t, m ~ 1/t), so it tests only band SHAPE: t'/t, t''/t, filling.
Here: recompute (m/m_e)K for several published LSCO-type bands at n = 0.76 to get the spread of b."""
import numpy as np
a = 3.78; hb2_me = 7.6200  # hbar^2/m_e in eV A^2
E0 = hb2_me / a**2          # hbar^2/(m_e a^2) in eV
def band(kx, ky, t, tp, tpp):
    e = -2*t*(np.cos(kx)+np.cos(ky)) - 4*tp*np.cos(kx)*np.cos(ky) - 2*tpp*(np.cos(2*kx)+np.cos(2*ky))
    vx = 2*t*np.sin(kx) + 4*tp*np.sin(kx)*np.cos(ky) + 4*tpp*np.sin(2*kx)
    return e, vx
def props(tp_t, tpp_t, n_target=0.76, t=1.0, N=1200, eta=0.01):
    k = (np.arange(N)+0.5)/N*2*np.pi - np.pi
    KX, KY = np.meshgrid(k, k)
    e, vx = band(KX, KY, t, tp_t*t, tpp_t*t)
    es = np.sort(e.ravel()); mu = es[int(n_target/2*es.size)]   # n per Cu, both spins
    w = np.exp(-((e-mu)/eta)**2)/(np.sqrt(np.pi)*eta)              # broadened delta
    Ntot = 2*w.mean()            # states / (energy * cell), both spins
    K = 2*(vx**2*w).mean()       # (hbar/e)^2 Phi_2D, energy units
    m_me = np.pi*Ntot*E0         # m = pi hbar^2 N_tot / a^2
    return Ntot, K, m_me
print("band (t'/t, t''/t)               N(0)/(eV a^2) at t=0.3   K (meV) at t=0.3   (m/m_e)K (meV)   b = 583/(m/m_e)K")
cases = [("Michon TB (-0.17, 0)", -0.17, 0.0),
         ("(-0.12, 0)", -0.12, 0.0), ("(-0.20, 0)", -0.20, 0.0), ("(-0.25, 0)", -0.25, 0.0),
         ("Yoshida-type (-0.15, +0.075)", -0.15, 0.075),
         ("Nd-LSCO ADMR band (-0.1364, +0.0682)", -0.1364, 0.0682)]
for name, tp, tpp in cases:
    Ntot, K, m_me = props(tp, tpp, t=0.3)
    prod = m_me*K*1e3
    print(f"{name:38s} {Ntot:6.3f}              {K*1e3:6.0f}             {prod:6.0f}           {583/prod:5.2f}")

# --- extra cases: small |t'| (what conformal b = 1.5 would need), and the Nd-LSCO band with its kz term ---
for name, tp in [("(0, 0)", 0.0), ("(-0.05, 0)", -0.05), ("(-0.08, 0)", -0.08)]:
    Ntot, K, m_me = props(tp, 0.0, t=0.3)
    print(f"{name:38s} {Ntot:6.3f}              {K*1e3:6.0f}             {m_me*K*1e3:6.0f}           {583/(m_me*K*1e3):5.2f}")
def nd_lsco_kz(N=500, Nz=16, eta=0.01):
    t, tp, tpp, tz = 1.0, -0.1364, 0.0682, 0.0651   # Grissonnanche 2021 band, units of t
    k = (np.arange(N)+0.5)/N*2*np.pi - np.pi
    KX, KY = np.meshgrid(k, k); Ns, Ks, es_all = [], [], []
    for kz in (np.arange(Nz)+0.5)/Nz*2*np.pi - np.pi:
        e, vx = band(KX, KY, t, tp, tpp)
        cz = 2*tz*np.cos(kz/2)
        e = e + cz*np.cos(KX/2)*np.cos(KY/2)*(np.cos(KX)-np.cos(KY))**2
        vx = vx + cz*(-0.5*np.sin(KX/2)*np.cos(KY/2)*(np.cos(KX)-np.cos(KY))**2
                      - 2*np.cos(KX/2)*np.cos(KY/2)*(np.cos(KX)-np.cos(KY))*np.sin(KX))
        es_all.append((e, vx))
    E = np.array([x[0] for x in es_all]); V = np.array([x[1] for x in es_all])
    mu = np.sort(E.ravel())[int(0.76/2*E.size)]
    w = np.exp(-((E-mu)/eta)**2)/(np.sqrt(np.pi)*eta)
    Ntot = 2*w.mean(); K = 2*(V**2*w).mean(); return Ntot, K, mu
Ntot, K, mu = nd_lsco_kz()
# rescale to t = 0.3 eV for comparison (product is t-independent)
prod = np.pi*(Ntot/0.3)*E0*K*0.3*1e3
print(f"{'Nd-LSCO band with kz (tz=0.0651t)':38s} {Ntot/0.3:6.3f}              {K*0.3e3:6.0f}             {prod:6.0f}           {583/prod:5.2f}   (mu/t at n=0.76: {mu:.4f}; Grissonnanche -0.8243)")
