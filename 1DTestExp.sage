# controlepaneel
featureDim = 2
graad = 6

# stel heenmatrix op (H)
vrs = []
hevel = []
for i in range(0,featureDim):
    lv = var("H_"+str(i))
    vrs.append(lv)
    hevel.append(lv)
H = matrix(hevel).transpose()

# stel heenbias op (h)
hevel = []
for i in range(0,featureDim):
    lv = var("h_"+str(i))
    vrs.append(lv)
    hevel.append(lv)
h = matrix(hevel).transpose()

# stel middenmodel op (M)
tvrs = []
mm = 1
for i in range(0,featureDim):
    mh = 0
    tv = var("t_" + str(i))
    tvrs.append(tv)
    for j in range(0,graad):
        lv = var("M_"+str(i)+"_"+str(j))
        mh += lv * chebyshev_T(j,tv)
        vrs.append(lv)
    mm *= mh

    # Geen terugtransformatie - model zorgt hiervoor!
# # stel terugmatrix op (T)
# hevel = []
# for i in range(0,featureDim):
#     lv = var("T_"+str(i))
#     vrs.append(lv)
#     hevel.append(lv)
# T = matrix(hevel)
#
# # stel hterugbias op (t)
# hevel = []
# for i in range(0,featureDim):
#     lv = var("t_"+str(i))
#     vrs.append(lv)
#     hevel.append(lv)
# t = matrix(hevel).transpose()

# tests

var('x')

heen = H*matrix([x]) + h
model = mm
for i, t in enumerate(tvrs):
    model = model.subs({t:heen[i][0]})
