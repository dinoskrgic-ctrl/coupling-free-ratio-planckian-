import sys, json, math
import numpy as np
import pymupdf as f
d=f.open(sys.argv[1])
cal={9:dict(y0=512.0,y25=334.5,x1=391.0,x10=526.5,name='Eu'),
     11:dict(y0=519.0,y25=338.0,x1=387.0,x10=526.0,name='Nd')}
cols={'p24':(0.87,0.03,0.02),'p21or22':(1.0,0.55,0.11)}
out={}
for pg,c in cal.items():
    p=d[pg]; W,H=p.rect.width,p.rect.height
    # tick calibration from black strokes: print short horizontal ticks near left axis
    for key,col in cols.items():
        pts=[]
        for dr in p.get_drawings():
            r=dr['rect']
            if not(r.x0>W*0.5 and H*0.38<r.y0<H*0.70): continue
            if not dr.get('color') or tuple(round(x,2) for x in dr['color'])!=col: continue
            for it in dr['items']:
                for q in it[1:]:
                    if hasattr(q,'x'): pts.append((q.x,q.y))
        pts=np.array(sorted(set(pts)))
        T=10**((pts[:,0]-c['x1'])/(c['x10']-c['x1']))
        CT=(c['y0']-pts[:,1])*25/(c['y0']-c['y25'])
        out[c['name']+'_'+key]=dict(T=T.tolist(),CT=CT.tolist())
        print(c['name'],key,len(T),'T range',T.min().round(3),T.max().round(2))
json.dump(out,open(sys.argv[2],'w'))
