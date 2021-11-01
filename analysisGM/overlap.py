#!/usr/bin/env python

def distance(x1, y1, x2, y2):
    d = ((x1 - x2)**2.0 + (y1 - y2)**2.0)**(0.5)
    return d

steplines = []
with open("INIT","rt") as file:
    for line in file:
        steplines.append(line.rstrip('\n'))

coords=[x.split() for x in steplines]

xcoords = []
ycoords = []
Rcoords = []
types   = []
for i in range(1,len(coords)):
    xcoords.append(float(coords[i][2]))
    ycoords.append(float(coords[i][3]))
    Rcoords.append(float(coords[i][1]))
    types.append(int(coords[i][0]))

for i in range(0,len(xcoords)):
    for j in range(i+1,len(xcoords)):
        if distance(xcoords[i],ycoords[i],xcoords[j],ycoords[j]) < Rcoords[i] + Rcoords[j]:
            print("overlap!!", i+2 , j+2, Rcoords[i], Rcoords[j],  distance(xcoords[i],ycoords[i],xcoords[j],ycoords[j]))
        #print(Rcoords[i],Rcoords[j])
        else:
            pass
            #print("no overlap")
