# controlepaneel
featureDim = 2
middenGraad = 5
heenGraad = 3
probleemVar = var('x')

# controledingetjes
middenGraadTeller = middenGraad+1
heenGraadTeller   = heenGraad+1

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

# tests
print(mm)
