#!/usr/bin/env python

import numpy as np
import param, math
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LinearLocator
from matplotlib.ticker import ScalarFormatter


def distance(x1, y1, x2, y2):
    d = ((x1 - x2)**2.0 + (y1 - y2)**2.0)**(0.5)
    return d

# for autoinit version
# current_LCG is current largest cluster radius
def boundaryOverlapCheck(coords, current_x, current_y, current_radius, current_LCG ):
    overlap = False

    # check whether is the beginning
    if current_LCG <= 0:
        return overlap

    # check whether this LCG is smaller than current radius
    if current_LCG <= current_radius :
        current_LCG = current_radius
    
    isRight = False
    isLeft = False
    isUp = False
    isDown = False

    # determine which kinds overlap is happening here
    overlapPotentialType = []

    # right most case, left most, up case, down case
    if current_x + current_radius + current_LCG >=  param.maxx - current_y/np.sqrt(3): # shifting in x boundary
        isRight = True
        overlapPotentialType.append( "right")
    if current_x - current_radius - current_LCG<= param.minx - current_y/np.sqrt(3): # shifting in x boundary
        isLeft = True
        overlapPotentialType.append("left" )
    if current_y + current_radius + current_LCG>= param.maxy:
        isUp = True
        overlapPotentialType.append("up")
    if current_y - current_radius - current_LCG<= param.miny:
        isDown = True
        overlapPotentialType.append("down")
        
    # no boundary overlap case
    if overlapPotentialType == [] :
        return overlap
    
    #print(overlapPotentialType)
    
    # newClusterList[][0] = x
    # newClusterList[][1] = y
    # newClusterList[][2] = radius
    # newClusterList[][3] = types
    newClusterList = []
    
    for types in overlapPotentialType:
        if types == "right" :
            # delta is the longest length that other potential overlap cluster can exist
            delta = current_x + current_radius - param.maxx/2 + 2 * current_LCG
            for cluster in coords:
                if cluster[0] <=  param.minx - cluster[1]/np.sqrt(3) + delta:
                    newClusterList.append([cluster[0] + param.maxx, cluster[1], cluster[2], cluster[3]])       
        if types == "left" :
            # this delt should be a negative number
            delta = current_x - current_radius + param.minx - 2 * current_LCG
            for cluster in coords:
                if cluster[0] >= param.maxx - cluster[1]/np.sqrt(3) + delta:
                    newClusterList.append([cluster[0] - param.maxx, cluster[1], cluster[2], cluster[3]])
        if types == "up" :
            delta = current_y + current_radius - param.maxy + 2 * current_LCG
            for cluster in coords:
                if cluster[1] <= param.miny + delta:
                    newClusterList.append([cluster[0] - param.maxy/np.sqrt(3), cluster[1] + param.maxy, cluster[2], cluster[3]])
        if types == "down" :
            delta = current_y - current_radius + param.miny - 2 * current_LCG
            for cluster in coords:
                if cluster[1] >= param.maxy +delta:
                    newClusterList.append([cluster[0] + param.maxy/np.sqrt(3), cluster[1] - param.maxy, cluster[2],cluster[3]])
    
    # check for corner cases
    if isRight and isUp :
        for cluster in coords:
                if (cluster[0] <= current_LCG) and (cluster[1] <= current_LCG):
                    newClusterList.append([cluster[0] + param.maxx - param.maxy/np.sqrt(3), cluster[1] + param.maxy, cluster[2], cluster[3]])
    if isRight and isDown :
        for cluster in coords:
                if (cluster[0] <= current_LCG) and (cluster[1] >= param.maxy - current_LCG):
                    newClusterList.append([cluster[0] + param.maxx + param.maxy/np.sqrt(3), cluster[1] - param.maxy, cluster[2], cluster[3]])
    if isLeft and isUp :
        for cluster in coords:
                if (cluster[0] >= param.maxx - current_LCG) and (cluster[1] <= current_LCG):
                    newClusterList.append([cluster[0] - param.maxx - param.maxy/np.sqrt(3), cluster[1] + param.maxy, cluster[2], cluster[3]])
    if isLeft and isDown :
         for cluster in coords:
                if (cluster[0] >= param.maxx - param.maxy/np.sqrt(3) - current_LCG) and (cluster[1] >= param.maxy - current_LCG):
                    newClusterList.append([cluster[0] - param.maxx + param.maxy/np.sqrt(3), cluster[1] - param.maxy, cluster[2], cluster[3]])
    
    # checking overlap
    for testCluster in newClusterList:
        if distance(current_x,current_y, testCluster[0],testCluster[1]) <= current_radius + testCluster[2]:
            overlap = True
            break
                
    return overlap 




numprimcell   = param.numprimecellFactor*param.num_clust + param.num_single_atom
xdups = ydups = (int(numprimcell**0.5) + 2)
PES           = [] # potential energy surface element, x, y, z, E
Clusters      = [] # all possible cluster R and E
N_mesh        = param.N_mesh # total number of PES points
N_meshx       = N_meshy = (N_mesh)**(0.5) # number of mesh points in x and y direction

if (param.num_clust == 0 and param.largest_cluster > 1):
    raise ValueError('Largest cluster must be 1 atom! set num_clust > 0')

with open('PES','r') as f1:
    data = f1.readlines()
    for line in data:   
        PES.append([ str(line.split()[0]), 
                     float(line.split()[1]), 
                     float(line.split()[2]), 
                     float(line.split()[3]) ])
with open('INIT', 'w') as f2:
    f2.write('%3s  %10s  %16s  %15s  %14s\n' % ('Pt', 'R', 'X', 'Y', 'E'))

with open('DATA','r') as f3:
    f3.readline()
    data = f3.readlines()
    for line in data:  
        Clusters.append([ int(line.split()[0]), 
                         float(line.split()[1]), 
                         float(line.split()[2]), 
                         float(line.split()[3]),
                         float(line.split()[4]) ])   

        
        
        
        
# Boltzmann population information, specific to Al2O3
T=750
hartree_to_ev = 27.2114
k_J           = 1.38064852e-23
k_ev          = 8.6173303e-5
k_hartree     = k_ev/hartree_to_ev
kT            = k_hartree*T
beta          = 1.0/kT

Pt2=[1.12859649,1.12859725,1.12928328,1.12928439,1.14117534,1.14123858,1.14179444,1.14181141]
Pt2_index=[0,1,2,3,4,5,6,7]
Pt2_boltz_weight=[math.exp(-Pt2[i]/kT) for i in range(0,len(Pt2))]
Q2=sum(Pt2_boltz_weight)
Pt2_prob=[Pt2_boltz_weight[i]/Q2 for i in range(0,len(Pt2_boltz_weight))]
Pt2_data=[[2.79100264,1.12859649],[2.79114648,1.12859725],[2.73305264,1.12928328],[2.78843079,1.12928439],[2.81726766,1.14117534], 
          [2.81883172,1.14123858], [2.79021837,1.14179444], [2.79231593,1.14181141]]

Pt3=[0.94983639, 0.95054819, 0.95061838, 0.95417737, 0.95418000, 0.95472066, 0.95949123, 0.95968591,0.95969046,
    0.95970682,0.96177741,0.96192135]
Pt3_index=[0,1,2,3,4,5,6,7,8,9,10,11]
Pt3_boltz_weight=[math.exp(-Pt3[i]/kT) for i in range(0,len(Pt3))]
Q3=sum(Pt3_boltz_weight)
Pt3_prob=[Pt3_boltz_weight[i]/Q3 for i in range(0,len(Pt3_boltz_weight))]
Pt3_data=[[2.98361137,0.94983639], [3.10320453,0.95054819], [2.97592016,0.95061838], [2.97440047,0.95417737], [2.97561209,0.95418000], 
          [2.97042544,0.95472066], [2.98699161,0.95949123], [2.98667235,0.95968591],[2.98633592,0.95969046], 
          [2.98702001,0.95970682],[3.13507762,0.96177741],[2.90420026,0.96192135]]

Pt4=[0.76719351,0.76744319,0.76744154,0.76744277,0.76884763,0.76884965]
Pt4_index=[0,1,2,3,4,5]
Pt4_boltz_weight=[math.exp(-Pt4[i]/kT) for i in range(0,len(Pt4))]
Q4=sum(Pt4_boltz_weight)
Pt4_prob=[Pt4_boltz_weight[i]/Q4 for i in range(0,len(Pt4_boltz_weight))]
Pt4_data=[[3.40871735,0.76719351], [3.40532512,0.76744319], [3.40600913,0.76744154], [3.40575678,0.76744277], [3.80705530,0.76884763], 
          [3.80781109,0.76884965]]

Pt5=[0.56915931,0.56916098,0.57035311,0.57117317,0.57185810,0.57278937,0.57352190,0.57449383,0.57599677,0.57894476,0.57898836,
    0.57901575,0.57923960,0.58003551,0.58007359,0.58073439,0.58127293,0.58153200,0.58207156]
Pt5_index=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
Pt5_boltz_weight=[math.exp(-Pt5[i]/kT) for i in range(0,len(Pt5))]
Q5=sum(Pt5_boltz_weight)
Pt5_prob=[Pt5_boltz_weight[i]/Q5 for i in range(0,len(Pt5_boltz_weight))]
Pt5_data=[[4.20564574,0.56915931], [4.20609141,0.56916098], [3.36972158,0.57035311], [4.20831736,0.57117317],[3.37496538,0.57185810], 
          [3.32086807,0.57278937], [3.41465927,0.57352190], [4.21767649,0.57449383], [3.55393743,0.57599677],[4.22288627,0.57894476],
          [4.07599753,0.57898836], [4.19531181,0.57901575], [4.13517505,0.57923960], [3.96464630,0.58003551], [3.95519852,0.58007359],
          [4.12011909,0.58073439],[3.85085560,0.58127293], [3.84258034,0.58153200], [4.19670411,0.58207156]]

Pt6=[0.37688221,0.38046857,0.38930451,0.38931809]
Pt6_index=[0,1,2,3]
Pt6_boltz_weight=[math.exp(-Pt6[i]/kT) for i in range(0,len(Pt6))]
Q6=sum(Pt6_boltz_weight)
Pt6_prob=[Pt6_boltz_weight[i]/Q6 for i in range(0,len(Pt6_boltz_weight))]
Pt6_data=[[4.28919697,0.37688221], [4.14839347,0.38046857], [4.17721450,0.38930451], [4.18055781,0.38931809]]

Pt7=[0.16949534, 0.17216772, 0.17781069, 0.17937308, 0.18004505, 0.18359426, 0.18386071, 0.18419310,0.18382938, 0.18042367, 0.17866330]
Pt7_index=[0,1,2,3,4,5,6,7,8,9,10]
Pt7_boltz_weight=[math.exp(-Pt7[i]/kT) for i in range(0,len(Pt7))]
Q7=sum(Pt7_boltz_weight)
Pt7_prob=[Pt7_boltz_weight[i]/Q7 for i in range(0,len(Pt7_boltz_weight))]
Pt7_data=[[4.27667410,0.16949534], [4.56780096,0.17216772], [4.23757940,0.17781069], [4.39487018,0.17937308], [4.22717566,0.18004505], 
          [4.98654557,0.18359426], [4.56976816,0.18386071], [4.63064528,0.18419310], [4.74374118,0.18382938], [4.50218286,0.18042367], 
          [4.36580910,0.17866330 ]]

Pt8=[0.00000000,0.01055526, 0.01332663]
Pt8_index=[0,1,2]
Pt8_boltz_weight=[math.exp(-Pt8[i]/kT) for i in range(0,len(Pt8))]
Q8=sum(Pt8_boltz_weight)
Pt8_prob=[Pt8_boltz_weight[i]/Q8 for i in range(0,len(Pt8_boltz_weight))]
Pt8_data=[[5.23684571,0.00000000],[4.61251355,0.01055526], [4.62601100,0.01332663]]

# for the case of pt2 to pt8
def boltzmannPopulationForNewCluster(numberOfAtoms) :
    assignedRadius = 0
    assignedEnergy = 0
    if numberOfAtoms==2:
        index = np.random.choice(Pt2_index,p=Pt2_prob)
        assignedRadius = Pt2_data[index][0]
        assignedEnergy = Pt2_data[index][1]
        return assignedRadius, assignedEnergy
    if numberOfAtoms==3:
        index = np.random.choice(Pt3_index,p=Pt3_prob)
        assignedRadius = Pt3_data[index][0]
        assignedEnergy = Pt3_data[index][1]
        return assignedRadius, assignedEnergy
    if numberOfAtoms==4:
        index = np.random.choice(Pt4_index,p=Pt4_prob)
        assignedRadius = Pt4_data[index][0]
        assignedEnergy = Pt4_data[index][1]
        return assignedRadius, assignedEnergy
    if numberOfAtoms==5:
        index = np.random.choice(Pt5_index,p=Pt5_prob)
        assignedRadius = Pt5_data[index][0]
        assignedEnergy = Pt5_data[index][1]
        return assignedRadius, assignedEnergy
    if numberOfAtoms==6:
        index = np.random.choice(Pt6_index,p=Pt6_prob)
        assignedRadius = Pt6_data[index][0]
        assignedEnergy = Pt6_data[index][1]
        return assignedRadius, assignedEnergy
    if numberOfAtoms==7:
        index = np.random.choice(Pt7_index,p=Pt7_prob)
        assignedRadius = Pt7_data[index][0]
        assignedEnergy = Pt7_data[index][1]
        return assignedRadius, assignedEnergy
    if numberOfAtoms==8:
        index = np.random.choice(Pt8_index,p=Pt8_prob)
        assignedRadius = Pt8_data[index][0]
        assignedEnergy = Pt8_data[index][1]
        return assignedRadius, assignedEnergy

    return assignedRadius, assignedEnergy


np.random.seed()  # seed for random number generator
count = 0
counter = 0
coords = []
LCG = 0 # LCG is largest cluster radius

# this is correct only if param.largest_cluster <=9
while (count < param.num_clust):
    
    R, E = boltzmannPopulationForNewCluster(param.largest_cluster)
      
    # old method, not as accurate as using the PES data directly
    #x_normal = param.xstep_max * np.random.randint(0, (N_meshx - 1) * (xdups + 1)) 
    #y_normal = param.ystep_max * np.random.randint(0, (N_meshy - 1) * (ydups + 1)) 
    #y = y_normal * np.sin(np.deg2rad(60))
    #x = x_normal - y_normal * np.cos(np.deg2rad(60))
    
    # new method, using the PES data directly
    temp2 = np.random.randint(0, len(PES))
    yIncreaseFactor =  param.primcell_b * np.random.randint(0, ydups + 1)
    y = PES[temp2][2] + yIncreaseFactor * np.sqrt(3)/2
    x = PES[temp2][1] + param.primcell_a * np.random.randint(0, xdups + 1) -  yIncreaseFactor* 1 /2
    
    overlap = False
    
    for i in range(len(coords)):
        if (distance(x, y, coords[i][0], coords[i][1]) <= coords[i][2] + R ):
            overlap = True # overlaping atoms

    # boundary overlap case check
    if overlap == False:
        overlap = boundaryOverlapCheck(coords, x, y, R, LCG)

    if (overlap == False):
        if LCG <= R:
            LCG = R
        coords.append([x, y, R, param.largest_cluster])
        with open('INIT', 'a') as f2:
            f2.write('%3i %16.8f %16.8f %16.8f %16.8f \n' %  
                    (param.largest_cluster, R, x, y, E))
        count += 1
    if (counter > param.CounterLimit):
        print('total num of generated clusters =', count)
        raise ValueError('small cell is used! use a larger cell!')
    counter += 1

count = 0
counter = 0
while (count < param.num_single_atom):
    temp = np.random.randint(0, len(PES))
    yIncreaseFactor =  param.primcell_b * np.random.randint(0, ydups + 1)
    y = PES[temp][2] + yIncreaseFactor * np.sqrt(3)/2
    x = PES[temp][1] + param.primcell_a * np.random.randint(0, xdups + 1) -  yIncreaseFactor* 1 /2
    E = PES[temp][3]
    overlap = False
    for i in range(len(coords)):
        if (distance(x, y, coords[i][0], coords[i][1]) <= coords[i][2] +  param.Ratom):
            overlap = True # overlaping atoms
            
    # boundary overlap case check
    if overlap == False:
        overlap = boundaryOverlapCheck(coords, x, y, param.Ratom, LCG)

    if (overlap == False):
        if LCG <= Clusters[temp][1]:
            LCG = Clusters[temp][1]
        coords.append([x, y, param.Ratom, 1])
        with open('INIT', 'a') as f2:
            f2.write('%3i %16.8f %16.8f %16.8f %16.8f \n' %
                    (1, param.Ratom, x, y, E))
        count += 1
    if (counter > param.CounterLimit):
        print('total num of generated single atoms =', count)
        raise ValueError('small cell is used! use a larger cell!')
    counter += 1

if (counter <= param.CounterLimit):
    print('total num of generated clusters =', param.num_clust + param.num_single_atom)
    print('set maxx in param.py =', (xdups+1), 'X', 'primcell_a')
    print('set maxy in param.py =', (ydups+1), 'X', 'primcell_b')

    # plot setting

    #minorLocator = AutoMinorLocator()
    #mlx  = MultipleLocator(param.xstep_max)
    #mly  = MultipleLocator(param.ystep_max)

    #title_font = {'fontname':'Times New Roman', 'size':'18', 'color':'black', 'weight':'normal', 'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    #axis_font = {'fontname':'Times New Roman', 'size':'18'}
    #mpl.rc('font',family='Times New Roman')

    #ax = plt.subplot() # Defines ax variable by creating an empty plot

    # Set the tick labels font
    #for label in (ax.get_yticklabels() + ax.get_xticklabels()):
        #label.set_fontname('Times New Roman')
        #label.set_fontsize(16)

    #ax.yaxis.set_minor_locator(mly)
    #ax.xaxis.set_minor_locator(mlx)

    xcoords = []
    ycoords = []
    Rcoords = []
    types   = []
    for i in range(len(coords)):
        xcoords.append(coords[i][0])
        ycoords.append(coords[i][1])
        Rcoords.append(30*int(coords[i][2])**2)
        types.append(int(coords[i][3]))

    #frame1 = plt.gca()
    #frame1.axes.yaxis.set_ticklabels([])
    #plt.yticks([])
    #frame1.axes.xaxis.set_ticklabels([])
    #plt.xticks([])
    #plt.ylim(0, param.maxy)
    #plt.xlim(- param.maxy / np.sqrt(3) , param.maxx)
    
    plt.figure(figsize=(7,7))
    boundaryDots_X = [0, param.maxx, param.maxx - param.maxy / np.sqrt(3), - param.maxy / np.sqrt(3), 0 ]
    boundaryDots_Y = [0, 0, param.maxy, param.maxy, 0]
    plt.plot(boundaryDots_X, boundaryDots_Y,marker="o", markerfacecolor="none", ms=1)
    
    for i,type in enumerate(types):
        x = xcoords[i]
        y = ycoords[i]
        plt.scatter(xcoords, ycoords, color='blue', s=Rcoords, marker='o')
        plt.text(x, y, type, fontsize=(12+type/3), color='yellow', horizontalalignment='center', 
                 verticalalignment='center')
    #plt.axvline(x=(xdups+1)*param.primcell_a, color='black', linestyle='-')
    #plt.axhline(y=(ydups+1)*param.primcell_b, color='black', linestyle='-')
    #plt.grid(which='minor', axis='both', linestyle='--')
    plt.show()
