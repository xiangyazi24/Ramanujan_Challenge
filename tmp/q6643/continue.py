import serial32 as P

def show(label, value):
    print(label, repr(value), flush=True)

show('BASE_CLEAN', P.run(150000, 0, (), P.clean()))
show('BASE_GARBAGE1', P.run(150000, 1))
show('BASE_GARBAGE2', P.run(150000, 2))
for sw in (('SHARED-LANE',), ('NO-SEAM-DRAIN',), ('NO-MIDPOINT-OWNER',), ('REMOTE-BUSY',), ('HEADLESS-REGION',), ('REMOTE-BUSY','HEADLESS-REGION')):
    X = P.clean()
    if sw == ('REMOTE-BUSY','HEADLESS-REGION'):
        g = P.glob(0,1,7); z = X[g][0]; X[g] = (z[:4] + ((0,1,1,0,0),), X[g][1])
        for l in range(16):
            g = P.glob(0,1,l); z = X[g][0]; X[g] = ((0,1) + z[2:], X[g][1])
    show('CONTROL ' + repr(sw), P.run(100000, 4, sw, X))
