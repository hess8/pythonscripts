from numpy import array

pnts =[array([ 0.        ,  0.45402743, -0.64209154]), 
       array([-0.27803373,  0.45402743, -0.64209154]), 
       array([-0.64874546,  0.45402743,  0.        ]), 
       array([-0.27803373,  0.45402743,  0.64209154])]

# labels = range(len(pnts))

ds = [0.64209154071986396, 0.69970301064027385, 0.64874545642786008, 
        0.69970301064027385]

# labels = [label for (d,label) in sorted(zip(ds,labels))]

print sorted(zip(ds, pnts), key=lambda x: x[0])

# pnts = [pnt for (d,pnt) in sorted(zip(ds,pnts))]

# print pnts