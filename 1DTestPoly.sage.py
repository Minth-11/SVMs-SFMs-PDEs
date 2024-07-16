

# This file was *autogenerated* from the file 1DTestPoly.sage
from sage.all_cmdline import *   # import sage library

_sage_const_4 = Integer(4); _sage_const_3 = Integer(3); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_2 = Integer(2); _sage_const_1000 = Integer(1000); _sage_const_2000 = Integer(2000)
import numpy as np
from sage.numerical.optimize import minimize
from scipy.optimize import line_search

# opstellen model
supports = _sage_const_4 
d = _sage_const_3 
bwds = []
vrs = []
mm = _sage_const_0 
for i in range(_sage_const_0 ,supports):
    wgt = var("alpha_"+str(i))
    vrs.append(wgt)
    # bwds.append((-1)^round(random()) * random())
    bwds.append(_sage_const_0 )
    vrt = var("x_"+str(i))
    vrs.append(vrt)
    bwds.append((-_sage_const_1 )**round(random()) * random())
    cl = var("c_"+str(i))
    vrs.append(cl)
    # bwds.append( (-1)^round(random()) * random() )
    bwds.append( _sage_const_0  )
    mm += wgt * (vrt*x  + cl)**d

model = mm
# print()
# print("model:")
# print(model)

# probleem
probleemVar = var('x')
f = model
probleem = [ 
                # diff(f,x) == f,
                # f.subs({x:0}) == 1
                diff(f(x),x) == (f(x) * (_sage_const_1  - f(x)))*_sage_const_4 *(_sage_const_3 /_sage_const_2 ),
                f.subs({x:_sage_const_0 }) == _sage_const_1 /_sage_const_2 
            ]
bereik = (x,-_sage_const_1 ,_sage_const_1 )

# controledingetjes

# opstellen loss
loss = _sage_const_0 
for i in probleem:
    loss += (i.lhs() - i.rhs())**_sage_const_2 
    
# print("loss:")
# print(loss)

def evalLoss(pt):
    dt = dict(zip(vrs,pt))
    lni = loss.subs(dt)
    # uit = integral(lni,bereik)
    # print(uit)
    # uit = monte_carlo_integral(lambda v: lni.subs({x:v}),[-1,1],10000)
    uit = monte_carlo_integral(lni,[-_sage_const_1 ],[_sage_const_1 ],_sage_const_1000 ,algorithm='vegas')
    print(uit,end='\r')
    return uit[_sage_const_0 ]



modelParameters = minimize(loss,bwds,verbose=True,maxiter=_sage_const_2000 )

print()

fnl = model.subs(dict(zip(vrs,modelParameters)))

plot(fnl,(x,-_sage_const_1 ,_sage_const_1 ))

# def maakOplFn(mdl,vrs,pt):
    # sbs = dict(zip(vrs,pt))
    # mh = model.subs(sbs)
    # return mh

# def vindFt(loss,fn,modelFi,bds):
    # l = loss.substitute_function(fn == modelFi)
    # # l = loss.substitute_function(fn, modelFi)
    # return integrate(l,bds)

# def diffComp(vr,mdl,ls,pt,fn,bds):
    # sbs = dict(zip(vrs,pt))
    # ep = sbs[vr]
    # sbs[vr] = vr
    # mm = mdl.subs(sbs)
    # lfni = ls.substitute_function(fn == mm)
    # ls = integrate(lfni,bds)
    # afgl = diff(ls,vr)
    # return afgl.subs({vr:ep})

# def lossGrd(vrs,mdl,ls,pt,fn,bds):
    # grd = []
    # li = len(vrs)
    # teller = -1
    # for vr in vrs:
        # teller += 1
        # print("grd: " + str(teller) + "/" + str(li),end='\r')
        # hevel = diffComp(vr,mdl,ls,pt,fn,bds)
        # grd.append(hevel)
    # return grd

# # print( lossGrd(vrs,model,loss,bwds,f(x),bereik) )

# def fnc(pt):
    # lpt = list(pt)
    # modelFi = maakOplFn(model,vrs,lpt)
    # uit = float(vindFt(loss,of,modelFi,bereik))
    # print("loss: " + str(uit))
    # return uit

# def grd(pt):
    # lpt = list(pt)
    # gl = lossGrd(vrs,model,loss,lpt,f(x),bereik)
    # uit = np.array([float(i) for i in gl])
    # print(uit)
    # return uit

# bp = np.array([float(i) for i in bwds])
# modelParameters = minimize(fnc,bp,gradient=grd,verbose=True,maxiter=100000)
# modelParameters = minimize(fnc,bp,verbose=True,maxiter=100000)
# # modelParameters = minimize(fnc,bp,verbose=True)
# # modelParameters = line_search(fnc,grd,bp,0.00001*grd(bp))
# lmp = list(modelParameters)
# fml = maakOplFn(model,vrs,lmp)
# plot(fml,(x,-1,1))
