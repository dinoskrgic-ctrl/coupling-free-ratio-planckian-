import sys, json, os
import numpy as np
# [repro] I/O only: the original required `fit_michon2019_B.py PDF curves.json out.json` and always opened the PDF
# for the calibration check. Here every argument is optional: PDF '-' (or missing file) skips the calibration
# check; curves default to ../data/michon2019_curves.json; output defaults to ../outputs/michon2019_B_fits.json.
_H=os.path.dirname(os.path.abspath(__file__))
_a=sys.argv[1:]+[None]*3
_pdf=_a[0] if _a[0] not in (None,'-') else None
_cur=_a[1] or os.path.join(_H,'..','data','michon2019_curves.json')
_out=_a[2] or os.path.join(_H,'..','outputs','michon2019_B_fits.json')
cur=json.load(open(_cur))
if _pdf and os.path.exists(_pdf):
    import pymupdf as f
    d=f.open(_pdf)
else:
    d=None; print('calibration check skipped (no arXiv:1804.08502 PDF given)')
# calibration check: reference horizontal lines (Eu p=0.16 green = 4.2; Nd p=0.12 blue = 3.6)
for pg,col,y0,y25,lab in ([(9,(0.03,0.38,0.1),512.0,334.5,'Eu p=0.16 (expect 4.2)'),(11,(0.0,0.0,0.83),519.0,338.0,'Nd p=0.12 (expect 3.6)')] if d is not None else []):
    p=d[pg]; W,H=p.rect.width,p.rect.height
    for dr in p.get_drawings():
        r=dr['rect']
        if r.x0>W*0.5 and H*0.38<r.y0<H*0.70 and dr.get('color') and tuple(round(x,2) for x in dr['color'])==col:
            print(lab, round((y0-r.y0)*25/(y0-y25),2), round((y0-r.y1)*25/(y0-y25),2))
def fit(T,C,lo,hi):
    m=(T>=lo)&(T<=hi); x=np.log(T[m]); A=np.vstack([-x,np.ones_like(x)]).T
    (B,c),*_=np.linalg.lstsq(A,C[m],rcond=None); return B, np.exp(c/B)
def at(T,C,t,w=0.04):
    m=np.abs(np.log(T/t))<w; return C[m].mean()
res={}
for k in ['Eu_p24','Nd_p24','Eu_p21or22','Nd_p21or22']:
    T=np.array(cur[k]['T']); C=np.array(cur[k]['CT'])
    row={}
    for lo,hi in [(0.55,10),(0.55,2),(2,10),(1.8,10)]:
        if T.min()<=lo*1.05:
            B,T0=fit(T,C,lo,hi); row[f'B_{lo}-{hi}K']=round(B,3)
    row['C/T@0.55']=round(at(T,C,0.55),2) if T.min()<0.56 else None
    row['C/T@2']=round(at(T,C,2.0),2); row['C/T@10']=round(at(T,C,9.6),2)
    res[k]=row; print(k,row)
json.dump(res,open(_out,"w"),indent=1)  # [repro] was open(sys.argv[3],"w")
