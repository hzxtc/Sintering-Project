#!/usr/bin/env python

import csv

init = []
with open("INIT","rt") as file:
    for line in file:
        init.append(line)

splitinit=[init[x].split() for x in range(0,len(init))]

GM_data=[[2,'3.0218605','1.09778751'],[3,'3.0040290','0.90056183'],[4,'3.6684145','0.72210192'],[5,'3.7748643','0.54175170'],[6,'4.1423044','0.36589531'],[7,'4.2377875','0.17865857'],[8,'4.6229707','0.00000000']]

for i in range(0,len(splitinit)):
    if splitinit[i][0] == "2":
        splitinit[i][1] = GM_data[0][1]
        splitinit[i][4] = GM_data[0][2]
    elif splitinit[i][0] == "3":
        splitinit[i][1] = GM_data[1][1]
        splitinit[i][4] = GM_data[1][2]
    elif splitinit[i][0] == "4":
        splitinit[i][1] = GM_data[2][1]
        splitinit[i][4] = GM_data[2][2]
    elif splitinit[i][0] == "5":
        splitinit[i][1] = GM_data[3][1]
        splitinit[i][4] = GM_data[3][2]
    elif splitinit[i][0] == "6":
        splitinit[i][1] = GM_data[4][1]
        splitinit[i][4] = GM_data[4][2]
    elif splitinit[i][0] == "7":
        splitinit[i][1] = GM_data[5][1]
        splitinit[i][4] = GM_data[5][2]
    elif splitinit[i][0] == "8":
        splitinit[i][1] = GM_data[6][1]
        splitinit[i][4] = GM_data[6][2]


with open("GM_INIT","w") as f:
    wr=csv.writer(f,delimiter='\t')
    wr.writerows(splitinit)
