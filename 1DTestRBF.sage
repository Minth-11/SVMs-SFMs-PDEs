import numpy as np
from sage.numerical.optimize import minimize
from scipy.optimize import line_search

# opstellen model
supports = 64
bwds = []
vrs = []
mm = 0
for i in range(0,supports):
    wgt = var("alpha_"+str(i))
    vrs.append(wgt)
    bwds.append((-1)^round(random()) * random())
    # bwds.append( 0 )
    vrt = var("x_"+str(i))
    vrs.append(vrt)
    bwds.append((-1)^round(random()) * random())
    sig = var("sigma_"+str(i))
    vrs.append(sig)
    bwds.append( random() * 10 )
    mm += wgt * exp( -1 * (x - vrt)^2 / (2*sig*sig))

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
                diff(f(x),x) == (f(x) * (1 - f(x)))*4*(3/2),
                f.subs({x:0}) == 1/2
            ]
bereik = (x,-2,2)

# controledingetjes

# opstellen loss
loss = 0
for i in probleem:
    loss += (i.lhs() - i.rhs())^2
    
# print("loss:")
# print(loss)
tlr = 0
def evalLoss(pt):
    lni = loss.subs(dict(zip(vrs,pt)))
    # print(lni)
    # uit = integral(lni,bereik)
    # uit = monte_carlo_integral(lambda v: lni.subs({x:v}),[-1,1],10000)
    uit = monte_carlo_integral(lni,[bereik[1]],[bereik[2]],10000,algorithm='vegas')
    global tlr
    tlr += 1
    print(str(tlr)+" "+str(uit),end='\r')
    return uit[0]

modelParameters = minimize(evalLoss,bwds,verbose=True,maxiter=10^6)

print()

fnl = model.subs(dict(zip(vrs,modelParameters)))

plot(fnl,(x,-2,2))

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
