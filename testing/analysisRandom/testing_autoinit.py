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
    if current_x + current_radius + current_LCG >=  param.maxx:
        isRight = True
        overlapPotentialType.append( "right")
    if current_x - current_radius - current_LCG<= param.minx:
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
    
    # newClusterList[][0] = x
    # newClusterList[][1] = y
    # newClusterList[][2] = radius
    # newClusterList[][3] = types
    newClusterList = []
    
    for types in overlapPotentialType:
        if types == "right" :
            # delta is the longest length that other potential overlap cluster can exist
            delta = current_x + current_radius - param.maxx + 2 * current_LCG
            for cluster in coords:
                if cluster[0] <=  param.minx + delta:
                    newClusterList.append([cluster[0] + param.maxx, cluster[1], cluster[2], cluster[3]])       
        if types == "left" :
            # this delt should be a negative number
            delta = current_x - current_radius + param.minx - 2 * current_LCG
            for cluster in coords:
                if cluster[0] >= param.maxx + delta:
                    newClusterList.append([cluster[0] - param.maxx, cluster[1], cluster[2], cluster[3]])
        if types == "up" :
            delta = current_y + current_radius - param.maxy + 2 * current_LCG
            for cluster in coords:
                if cluster[1] <= param.miny + delta:
                    newClusterList.append([cluster[0], cluster[1] +param.maxy, cluster[2], cluster[3]])
        if types == "down" :
            delta = current_y - current_radius + param.miny - 2 * current_LCG
            for cluster in coords:
                if cluster[1] >= param.maxy +delta:
                    newClusterList.append([cluster[0], cluster[1] - param.maxy, cluster[2],cluster[3]])
    
    # check for corner cases
    if isRight and isUp :
        for cluster in coords:
                if (cluster[0] <= current_LCG) and (cluster[1] <= current_LCG):
                    newClusterList.append([cluster[0] + param.maxx, cluster[1] + param.maxy, cluster[2], cluster[3]])
    if isRight and isDown :
        for cluster in coords:
                if (cluster[0] <= current_LCG) and (cluster[1] >= param.maxy - current_LCG):
                    newClusterList.append([cluster[0] + param.maxx, cluster[1] - param.maxy, cluster[2], cluster[3]])
    if isLeft and isUp :
        for cluster in coords:
                if (cluster[0] >= param.maxx - current_LCG) and (cluster[1] <= current_LCG):
                    newClusterList.append([cluster[0] - param.maxx, cluster[1] + param.maxy, cluster[2], cluster[3]])
    if isLeft and isDown :
         for cluster in coords:
                if (cluster[0] >= param.maxx - current_LCG) and (cluster[1] >= param.maxy - current_LCG):
                    newClusterList.append([cluster[0] - param.maxx, cluster[1] - param.maxy, cluster[2], cluster[3]])
    
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

np.random.seed()  # seed for random number generator
count = 0
counter = 0
coords = []
LCG = 0 # LCG is largest cluster radius

# random 
while (count < param.num_clust):
    if (param.largest_cluster == 2):
        temp = np.random.randint(0, 5)
    elif (param.largest_cluster == 3):
        temp = np.random.randint(5, 8)
    elif (param.largest_cluster == 4):
        temp = np.random.randint(8, 14)
    elif (param.largest_cluster == 5):
        temp = np.random.randint(14, 16)
    elif (param.largest_cluster == 6):
        temp = np.random.randint(16, 21)
    elif (param.largest_cluster == 7):
        temp = np.random.randint(21, 28)
    elif (param.largest_cluster == 8):
        temp = np.random.randint(28, 33)
    elif (param.largest_cluster > 8):
        temp = np.random.randint(24 + param.largest_cluster, 25 + param.largest_cluster)
    x = param.xstep_max * np.random.randint(0, (N_meshx - 1) * (xdups + 1)) 
    y = param.ystep_max * np.random.randint(0, (N_meshy - 1) * (ydups + 1)) 
    E = Clusters[temp][4]
    overlap = False
    for i in range(len(coords)):
        if (distance(x, y, coords[i][0], coords[i][1]) <= coords[i][2] + Clusters[temp][1] ):
            overlap = True # overlaping atoms

    # boundary overlap case check
    if overlap == False:
        overlap = boundaryOverlapCheck(coords, x, y, Clusters[temp][1], LCG)

    if (overlap == False):
        if LCG <= Clusters[temp][1]:
            LCG = Clusters[temp][1] 
        coords.append([x, y, Clusters[temp][1], Clusters[temp][0]])
        with open('INIT', 'a') as f2:
            f2.write('%3i %16.8f %16.8f %16.8f %16.8f \n' %  
                    (Clusters[temp][0], Clusters[temp][1], x, y, E))
        count += 1
    if (counter > param.CounterLimit):
        print('total num of generated clusters =', count)
        raise ValueError('small cell is used! use a larger cell!')
    counter += 1

count = 0
counter = 0
while (count < param.num_single_atom):
    temp = np.random.randint(0, len(PES))
    x = PES[temp][1] * np.random.randint(0, (xdups + 1))
    y = PES[temp][2] * np.random.randint(0, (ydups + 1))
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

    minorLocator = AutoMinorLocator()
    mlx  = MultipleLocator(param.xstep_max)
    mly  = MultipleLocator(param.ystep_max)

    title_font = {'fontname':'Times New Roman', 'size':'18', 'color':'black', 'weight':'normal',
                  'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    axis_font = {'fontname':'Times New Roman', 'size':'18'}
    mpl.rc('font',family='Times New Roman')

    ax = plt.subplot() # Defines ax variable by creating an empty plot

    # Set the tick labels font
    for label in (ax.get_yticklabels() + ax.get_xticklabels()):
        label.set_fontname('Times New Roman')
        label.set_fontsize(16)

    ax.yaxis.set_minor_locator(mly)
    ax.xaxis.set_minor_locator(mlx)

    xcoords = []
    ycoords = []
    Rcoords = []
    types   = []
    for i in range(len(coords)):
        xcoords.append(coords[i][0])
        ycoords.append(coords[i][1])
        Rcoords.append(30*int(coords[i][2])**2)
        types.append(int(coords[i][3]))

    frame1 = plt.gca()
    frame1.axes.yaxis.set_ticklabels([])
    plt.yticks([])
    frame1.axes.xaxis.set_ticklabels([])
    plt.xticks([])
    plt.ylim(0,(ydups+1)*param.primcell_b)
    plt.xlim(0,(xdups+1)*param.primcell_a)
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
