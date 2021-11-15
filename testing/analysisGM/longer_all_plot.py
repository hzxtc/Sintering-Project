import numpy as np
import matplotlib.pyplot as plt
import statistics as st

steplines = []
with open("metropolis1","rt") as file:
    for line in file:
        if line.startswith("s"):
            steplines.append(line.rstrip('\n'))
            

splitlines=[]
clustersize1=[]
for i in range(0,len(steplines)):
    splitlines.append(steplines[i].split())
    
for i in range(0,len(splitlines)):
    clustersize1.append(int(splitlines[i][5]))


steplines = []
with open("metropolis2","rt") as file:
    for line in file:
        if line.startswith("s"):
            steplines.append(line.rstrip('\n'))
            

splitlines=[]
clustersize2=[]
for i in range(0,len(steplines)):
    splitlines.append(steplines[i].split())
    
for i in range(0,len(splitlines)):
    clustersize2.append(int(splitlines[i][5]))


steplines = []
with open("metropolis3","rt") as file:
    for line in file:
        if line.startswith("s"):
            steplines.append(line.rstrip('\n'))
            

splitlines=[]
clustersize3=[]
for i in range(0,len(steplines)):
    splitlines.append(steplines[i].split())
    
for i in range(0,len(splitlines)):
    clustersize3.append(int(splitlines[i][5]))


steplines = []
with open("metropolis4","rt") as file:
    for line in file:
        if line.startswith("s"):
            steplines.append(line.rstrip('\n'))
            

splitlines=[]
clustersize4=[]
for i in range(0,len(steplines)):
    splitlines.append(steplines[i].split())
    
for i in range(0,len(splitlines)):
    clustersize4.append(int(splitlines[i][5]))


steplines = []
with open("metropolis5","rt") as file:
    for line in file:
        if line.startswith("s"):
            steplines.append(line.rstrip('\n'))
            
steps=[]
for i in range(0,100001):
    steps.append(int(i))

splitlines=[]
clustersize5=[]
for i in range(0,len(steplines)):
    splitlines.append(steplines[i].split())
    
for i in range(0,len(splitlines)):
    clustersize5.append(int(splitlines[i][5]))

steplines = []
with open("metropolis6","rt") as file:
    for line in file:
        if line.startswith("s"):
            steplines.append(line.rstrip('\n'))


splitlines=[]
clustersize6=[]
for i in range(0,len(steplines)):
    splitlines.append(steplines[i].split())

for i in range(0,len(splitlines)):
    clustersize6.append(int(splitlines[i][5]))


steplines = []
with open("metropolis7","rt") as file:
    for line in file:
        if line.startswith("s"):
            steplines.append(line.rstrip('\n'))


splitlines=[]
clustersize7=[]
for i in range(0,len(steplines)):
    splitlines.append(steplines[i].split())

for i in range(0,len(splitlines)):
    clustersize7.append(int(splitlines[i][5]))


steplines = []
with open("metropolis8","rt") as file:
    for line in file:
        if line.startswith("s"):
            steplines.append(line.rstrip('\n'))


splitlines=[]
clustersize8=[]
for i in range(0,len(steplines)):
    splitlines.append(steplines[i].split())

for i in range(0,len(splitlines)):
    clustersize8.append(int(splitlines[i][5]))

steplines = []
with open("metropolis9","rt") as file:
    for line in file:
        if line.startswith("s"):
            steplines.append(line.rstrip('\n'))


splitlines=[]
clustersize9=[]
for i in range(0,len(steplines)):
    splitlines.append(steplines[i].split())

for i in range(0,len(splitlines)):
    clustersize9.append(int(splitlines[i][5]))

steplines = []
with open("metropolis10","rt") as file:
    for line in file:
        if line.startswith("s"):
            steplines.append(line.rstrip('\n'))


splitlines=[]
clustersize10=[]
for i in range(0,len(steplines)):
    splitlines.append(steplines[i].split())

for i in range(0,len(splitlines)):
    clustersize10.append(int(splitlines[i][5]))



plt.plot(steps,clustersize1)
plt.plot(steps,clustersize2)
plt.plot(steps,clustersize3)#,label="isomer")#,c="g")
plt.plot(steps,clustersize4)#,c="b",label="GM + isomer")
plt.plot(steps,clustersize5)#,c="b")
plt.plot(steps,clustersize6,)#c="b")
plt.plot(steps,clustersize7)#,c="r",label="all isomer")
plt.plot(steps,clustersize8)#,c="r")
plt.plot(steps,clustersize9)#,c="r")
plt.plot(steps,clustersize10)
plt.title(input("plot title: "))
plt.xlabel("step number")
plt.ylabel("number of clusters")
#plt.legend()
plt.show()
