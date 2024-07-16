import numpy as np
from sage.numerical.optimize import minimize
from scipy.optimize import line_search

# controlepaneel
featureDim = 2
middenGraad = 15
heenGraad = 1
probleemVar = var('x')
f = function('f')
probleem = [ # logistic DE
                diff(f(x),x) == (f(x) * (1 - f(x)))*4*(3/2),
                f(0) == 1/2
            ]
bereik = (x,-1,1)
of = f(x)

# controledingetjes
middenGraadTeller = middenGraad+1
heenGraadTeller   = heenGraad+1

# opstellen loss
loss = 0
for i in probleem:
    loss += (i.lhs() - i.rhs())^2
print(loss)

# variabelen en beginpunt enz.
vrs = []
bwds = []

# stel heentransformatie op
heen = []
htvrs = []
for i in range(0,featureDim):
    hhevel = 0
    for j in range(0,heenGraadTeller):
        lv = var("H_"+str(i)+"_"+str(j))
        vrs.append(lv)
        if j == 1:
            bwds.append(1)
        else:
            bwds.append(0)
        hhevel += lv * chebyshev_T(j,probleemVar).horner(probleemVar)
    heen = heen + [hhevel]

# stel middenmodel op (M)
mtvrs = []
mm = 1
for i in range(0,featureDim):
    mh = 0
    tv = heen[i]
    for j in range(0,middenGraadTeller):
        lv = var("M_"+str(i)+"_"+str(j))
        mh += lv * chebyshev_T(j,tv)
        vrs.append(lv)
        if j == 0:
            bwds.append(1)
        else:
            bwds.append(0)
    mm *= mh

model = mm

def maakOplFn(mdl,vrs,pt):
    sbs = dict(zip(vrs,pt))
    mh = model.subs(sbs)
    return mh

def vindFt(loss,fn,modelFi,bds):
    l = loss.substitute_function(fn == modelFi)
    return integrate(l,bds)

def diffComp(vr,mdl,ls,pt,fn,bds):
    sbs = dict(zip(vrs,pt))
    ep = sbs[vr]
    sbs[vr] = vr
    mm = mdl.subs(sbs)
    lfni = ls.substitute_function(fn == mm)
    ls = integrate(lfni,bds)
    afgl = diff(ls,vr)
    return afgl.subs({vr:ep})


def lossGrd(vrs,mdl,ls,pt,fn,bds):
    grd = []
    li = len(vrs)
    teller = -1
    for vr in vrs:
        teller += 1
        print("grd: " + str(teller) + "/" + str(li),end='\r')
        hevel = diffComp(vr,mdl,ls,pt,fn,bds)
        grd.append(hevel)
    return grd

# print( lossGrd(vrs,model,loss,bwds,f(x),bereik) )

def fnc(pt):
    lpt = list(pt)
    modelFi = maakOplFn(model,vrs,lpt)
    uit = float(vindFt(loss,of,modelFi,bereik))
    print("loss: " + str(uit))
    return uit

def grd(pt):
    lpt = list(pt)
    gl = lossGrd(vrs,model,loss,lpt,f(x),bereik)
    uit = np.array([float(i) for i in gl])
    print(uit)
    return 0.001*uit

bp = np.array([float(i) for i in bwds])

# modelParameters  = minimize(fnc,bp,gradient=grd,verbose=True)
modelParameters = minimize(fnc,bp,verbose=True,maxiter=100000)
# modelParameters = line_search(fnc,grd,bp,0.00001*grd(bp))
lmp = list(modelParameters)
fml = maakOplFn(model,vrs,lmp)
plot(fml,(x,-1,1))
