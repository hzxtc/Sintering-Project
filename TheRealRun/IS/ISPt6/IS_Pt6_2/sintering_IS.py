#!/usr/bin/env python

'''
Sintering simulation via Ostwald Ripening of metalic clusters deposited on the surface 
based on NVT ensemble and metropolis moves.

Borna Zandkarimi 2020

Now, the sintering simulation is modifided for the case on Al2O3

Tom Hong 2021
'''

import numpy as np
import timeit, math, copy
import param
from scipy.special import expit          # for handling very small exp
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LinearLocator
from matplotlib.ticker import ScalarFormatter

# FUNCTIONS

def distance(x1, y1, x2, y2):
    r = ((x1 - x2)**2.0 + (y1 - y2)**2.0)**(0.5)
    return r

# for sintering version
# current_LCG is current largest cluster radius
def boundaryOverlapCheck(coords, current_x, current_y, current_radius, current_LCG ):
    overlap = False
    targetPosition = -1

    # check whether is the beginning
    if current_LCG <= 0:
        return overlap, targetPosition, current_LCG

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
        return overlap, targetPosition, current_LCG
    
    #print(overlapPotentialType)
    
    # newClusterList[][0] = x
    # newClusterList[][1] = y
    # newClusterList[][2] = radius
    # newClusterList[][3] = types
    # newClusterList[][4] = position
    newClusterList = []
    
    for types in overlapPotentialType:
        if types == "right" :
            # delta is the longest length that other potential overlap cluster can exist
            delta = current_x + current_radius - param.maxx/2 + 2 * current_LCG
            #counting position
            n = 0
            # for coords [0] = type, [1] = R, [2] = x, [3] = y, [4] = e
            for cluster in coords:
                if cluster[2] <=  param.minx - cluster[1]/np.sqrt(3) + delta:
                    newClusterList.append([cluster[2] + param.maxx, cluster[3], cluster[1], cluster[0],n])  
                n += 1
        if types == "left" :
            # this delt should be a negative number
            delta = current_x - current_radius + param.minx - 2 * current_LCG
            n = 0
            for cluster in coords:
                if cluster[2] >= param.maxx - cluster[1]/np.sqrt(3) + delta:
                    newClusterList.append([cluster[2] - param.maxx, cluster[3], cluster[1], cluster[0],n])
                n += 1
        if types == "up" :
            delta = current_y + current_radius - param.maxy + 2 * current_LCG
            n = 0
            for cluster in coords:
                if cluster[3] <= param.miny + delta:
                    newClusterList.append([cluster[2] - param.maxy/np.sqrt(3), cluster[3] + param.maxy, cluster[1], cluster[0],n])
                n += 1
        if types == "down" :
            delta = current_y - current_radius + param.miny - 2 * current_LCG
            n = 0
            for cluster in coords:
                if cluster[3] >= param.maxy +delta:
                    newClusterList.append([cluster[2] + param.maxy/np.sqrt(3), cluster[3] - param.maxy, cluster[1],cluster[0],n])
                n += 1
    
    # check for corner cases
    if isRight and isUp :
        n = 0
        for cluster in coords:
                if (cluster[2] <= current_LCG) and (cluster[3] <= current_LCG):
                    newClusterList.append([cluster[2] + param.maxx - param.maxy/np.sqrt(3), cluster[3] + param.maxy, cluster[1], cluster[0],n])
                n += 1
    if isRight and isDown :
        n = 0
        for cluster in coords:
                if (cluster[2] <= current_LCG) and (cluster[3] >= param.maxy - current_LCG):
                    newClusterList.append([cluster[2] + param.maxx + param.maxy/np.sqrt(3), cluster[3] - param.maxy, cluster[1], cluster[0],n])
                n += 1
    if isLeft and isUp :
        n = 0
        for cluster in coords:
                if (cluster[2] >= param.maxx - current_LCG) and (cluster[3] <= current_LCG):
                    newClusterList.append([cluster[2] - param.maxx - param.maxy/np.sqrt(3), cluster[3] + param.maxy, cluster[1], cluster[0],n])
                n += 1
    if isLeft and isDown :
        n = 0
        for cluster in coords:
                if (cluster[2] >= param.maxx - param.maxy/np.sqrt(3) - current_LCG) and (cluster[3] >= param.maxy - current_LCG):
                    newClusterList.append([cluster[0] - param.maxx + param.maxy/np.sqrt(3), cluster[1] - param.maxy, cluster[2], cluster[3],n])
                n += 1
    
    # checking overlap
    for testCluster in newClusterList:
        if distance(current_x,current_y, testCluster[0],testCluster[1]) <= current_radius + testCluster[2]:
            #print("labla: " + str(distance(current_x,current_y, testCluster[0],testCluster[1])))
            #print("current radius: " + str(current_radius))
            #print(testCluster[2])
            targetPosition = testCluster[4]
            overlap = True
            break
                
    return overlap, targetPosition, current_LCG

def overlap_check(Clusters, OUTPUT_data, LCG):
    ovlp = False
    idx = []
    idx_pair = []
    for i in range(len(OUTPUT_data) - 1):
        if (ovlp == True):
            break
        else:
            for j in range(i+1, len(OUTPUT_data)):
                
                R1 = OUTPUT_data[i][1] 
                R2 = OUTPUT_data[j][1] 
                d  = distance(OUTPUT_data[i][2],OUTPUT_data[i][3],OUTPUT_data[j][2],OUTPUT_data[j][3])
                if (d < R1 + R2):
                   # print("normal overlap")
                    ovlp = True
                    idx.append(i)
                    idx.append(j)
                    idx_pair.append([i,j])
                    #merge
                    numnew = OUTPUT_data[i][0] + OUTPUT_data[j][0]
                    if numnew > 8 :
                        for k in range(len(Clusters)):
                            if (Clusters[k][0] == numnew):
                                Rnew = Clusters[k][1]
                                Enew = Clusters[k][4]
                                if LCG < Rnew:
                                    LCG = Rnew
                    if numnew <= 8 :
                        Rnew, Enew = boltzmannPopulationForNewCluster(numnew)
                        if LCG < Rnew:
                            LCG = Rnew
                            
                    # determine the new location
                    # this is the old method, the new location would be the one who has the larger number of atoms 
                    #if (OUTPUT_data[i][0] > OUTPUT_data[j][0]):
                        #xnew = OUTPUT_data[i][2]
                        #ynew = OUTPUT_data[i][3]
                    #else:
                        #xnew = OUTPUT_data[j][2]
                        #ynew = OUTPUT_data[j][3]
                        
                    # here is the new method:
                    xnew, ynew = newLocationFinder(OUTPUT_data[i][2], OUTPUT_data[i][3], OUTPUT_data[i][0], OUTPUT_data[j][2], OUTPUT_data[j][3],  OUTPUT_data[j][0])
                    OUTPUT_data.append([numnew,Rnew,xnew,ynew,Enew])
                    break
            # overlap boundary check
            ovlp1 = False
            ovlp1, overlapPosition, current_LCG = boundaryOverlapCheck(OUTPUT_data, OUTPUT_data[i][2],  OUTPUT_data[i][3], OUTPUT_data[i][1], LCG)
            
            #print(str(ovlp1) + " " + str(overlapPosition) + " " + str(current_LCG))
            
            if ovlp1 == True:
                ovlp = True
                LCG = current_LCG
                j1 = overlapPosition
                
                # eliminate the possiblity that overlap itself
                if (i == j1):
                    ovlp = False
                    break

                idx.append(i)
                idx.append(j1)
                idx_pair.append([i,j1])
                #print("boundary overlap")
                #print("LCG: " + str(LCG))
                #print("x1: " + str(OUTPUT_data[i][2]) + " y1: " + str(OUTPUT_data[i][3])+ " r1: " + str(OUTPUT_data[i][1]) + " x2: " + str(OUTPUT_data[j1][2]) + " y2: " + str(OUTPUT_data[j1][3])+ " r2: " + str(OUTPUT_data[j1][1]))
                #merge
                numnew = OUTPUT_data[i][0] + OUTPUT_data[j1][0]
                if numnew > 8 :
                    for k in range(len(Clusters)):    
                        if (Clusters[k][0] == numnew):
                            Rnew = Clusters[k][1]
                            Enew = Clusters[k][4]
                            if LCG < Rnew:
                                LCG = Rnew
                if numnew <= 8 :
                    Rnew, Enew = boltzmannPopulationForNewCluster(numnew)
                    if LCG < Rnew:
                        LCG = Rnew
                # determine the new location
                # this is the old method, the new location would be the one who has the larger number of atoms        
                if (OUTPUT_data[i][0] > OUTPUT_data[j1][0]):
                    xnew = OUTPUT_data[i][2]
                    ynew = OUTPUT_data[i][3]
                else:
                    xnew = OUTPUT_data[j1][2]
                    ynew = OUTPUT_data[j1][3]
                    
                    
                # here is the new method: which is not correct for boundary overlap
                #xnew, ynew = newLocationFinder(OUTPUT_data[i][2], OUTPUT_data[i][3], OUTPUT_data[i][0], OUTPUT_data[j1][2], OUTPUT_data[j1][3],  OUTPUT_data[j1][0])
                OUTPUT_data.append([numnew,Rnew,xnew,ynew,Enew])
                break
            
    #removing duplicates from the list
   # idx = list(dict.fromkeys(idx))
    #idx = list(dict.fromkeys(idx))
    
   # for i in sorted(idx, reverse=True):
       # del (OUTPUT_data[i])

    # new remove methond
    if (ovlp == True):
        if (idx[0]<idx[1]):
            OUTPUT_data.pop(idx[0]) # remove i 
            OUTPUT_data.pop(idx[1]-1) # remove j -1
        if (idx[0]>idx[1]):
            OUTPUT_data.pop(idx[0]) # remove i 
            OUTPUT_data.pop(idx[1]) # remove j  
    return OUTPUT_data, ovlp, idx_pair, LCG


def PES_finder(X_new, Y_new, PES_copy, N_mesh):
    finder = False
    n = 0
    
    #print(str(X_new) + "    "+ str(Y_new))
    
    for ii in range(N_mesh - 1):
        #print("count" + str(n))
        n = n + 1
        # old code
        #dx =  np.abs( param.xstep_max*( math.ceil(X_new/param.xstep_max) % N_meshx ) - PES_copy[ii][1] )
        #dy =  np.abs( param.ystep_max*( math.ceil(Y_new/param.ystep_max) % N_meshy ) - PES_copy[ii][2] )
        dy =  np.abs(  Y_new % (param.primcell_b *  np.sqrt(3)/2)- PES_copy[ii][2] )
        yShiftForX = Y_new - ( Y_new % (param.primcell_b *  np.sqrt(3)/2))
        dx =  np.abs( ((X_new +  yShiftForX /np.sqrt(3))% param.primcell_a) + 0.0000001 - PES_copy[ii][1] ) % param.primcell_a # extra mode for negative x case
        # we need to add 0.0000001 for data correction, sometimes the data will lost some very small number after storage and we need to correct it manually in order to fix the calculation
        #print(( Y_new % (param.primcell_b *  np.sqrt(3)/2)))
        #print(str(n) + ": " + str(dx) + ' ' + str(dy))
        
        if ( dx  < 0.1 and dy  < 0.1 ):    
            # PES found 
            # print("find!")
            E_PES  = PES_copy[ii][3] 
            finder = True
            break    
    if (finder == False):
        raise ValueError('Could not find PES for single atom. Bad move! Check you initial setup')      
    return E_PES

def Cluster_finder(Clusters, OUTPUT_data, irmv, iadd, which='both'):
    prob1     = []
    prob2     = []
    prob1_tmp = []
    prob2_tmp = []
    for j in range(Ncluster_tot):     
        if (Clusters[j][0] == OUTPUT_data[iadd][0] + 1 and (which == 'add' or which == 'both') ) :
            if (prob1_tmp == []):
                E_min_add = Clusters[j][4]
                prob1_tmp.append(1.0)     
            else:
                # choose a larger cluster
                prob1_tmp.append( expit(-beta * (Clusters[j][4] - E_min_add) ))     
        if (Clusters[j][0] == OUTPUT_data[irmv][0] - 1 and (which == 'remove' or which == 'both') ) :
            if (prob2_tmp == []):
                E_min_rmv = Clusters[j][4]
                prob2_tmp.append(1.0)     
            else:
                # choose a smaller cluster
                prob2_tmp.append( expit(-beta * (Clusters[j][4] - E_min_rmv) ))    

    if ( which == 'add' or which == 'both' ): 
        prob1 = [icount / sum(prob1_tmp) for icount in prob1_tmp]
        # weighted random chosen number
        clust_idx_add = np.random.choice(np.arange(len(prob1)), 1, p=prob1, replace=False)     
        Enew_ad   = - np.log(prob1_tmp[clust_idx_add[0]]) / beta + E_min_add

     
    if ( which == 'remove' or which == 'both' ): 
        prob2 = [icount / sum(prob2_tmp) for icount in prob2_tmp]
        # weighted random chosen number
        clust_idx_rmv = np.random.choice(np.arange(len(prob2)), 1, p=prob2, replace=False)     
        Enew_remv = - np.log(prob2_tmp[clust_idx_rmv[0]]) / beta + E_min_rmv

    remove = False
    add    = False
    for k in range(Ncluster_tot):          
        if ( (Clusters[k][0] == OUTPUT_data[indx][0] - 1) and (remove==False) and (which == 'remove' or which == 'both') ):
            Rnew_remv = Clusters[k+clust_idx_rmv[0]][1]
            remove    = True            
        if ( (Clusters[k][0] == OUTPUT_data[i][0] + 1) and (add==False) and (which == 'add' or which == 'both') ):
            Rnew_ad  = Clusters[k+clust_idx_add[0]][1]
            add      = True 
 
     # return new energy and radius
    if (which == 'both'):
        return Enew_ad, Enew_remv, Rnew_ad, Rnew_remv
    elif (which == 'add'):
        return Enew_ad, Rnew_ad
    elif (which == 'remove'):
        return Enew_remv, Rnew_remv
 
 
def X_finder(y_coord,primcell,maxx): #this won't work for a supercell- have to play around with it! maybe do all rel maxx?
    #gives the minimum and maximum x values for a given y coordinate
    frac= y_coord /  (primcell * (np.sqrt(3)/2))
    min_x = frac * primcell * (-0.5)
    max_x = min_x + maxx
    return min_x, max_x

def newLocationFinder(x1, y1, type1, x2, y2, type2):
    # using the Parametrization of line by t to help us figure out the proposed positions
    t = type1 / (type1 + type2)
    
    proposedX = x2 + t * (x1 - x2)
    proposedY = y2 + t * (y1 - y2)
    
    shortestDistance =  distance (allPossibleLocation[0][0], allPossibleLocation[0][1], proposedX, proposedY)
    
    new_X = 0
    new_Y = 0
    
    # finding the nearest grid point that is close to the proposed position
    for location in allPossibleLocation:
        if (distance (location[0], location[1], proposedX, proposedY) < shortestDistance):
            shortestDistance = distance (location[0],location[1], proposedX, proposedY)
            new_X = location[0]
            new_Y = location[1]
            
    return new_X, new_Y



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

# READ INPUT   

PES           = [] # potential energy surface element, x, y, z, E
INIT_data     = [] # initial cluster R, x, y, E
Clusters      = [] # all possible cluster R and E              
beta          = param.beta 
hartree_to_ev = param.hartree_to_ev   
Ncluster = len(open('INIT').readlines(  )) - 1 # Number of clusters on the surface during simulation which can be changed
Ncluster_tot = len(open('DATA').readlines(  )) - 1 # Total number of unique different clusters
Metro_Max     = param.MMAX   # num of Metropolis steps
write_step    = param.wstep  # writing metropolis each wstep
N_mesh        = param.N_mesh # total number of PES points
N_meshx       = N_meshy = (N_mesh)**(0.5) # number of mesh points in x and y direction
choose_dir    = [[1,0],[1/2,np.sqrt(3)/2],[-1/2,np.sqrt(3)/2],[-1/2,-np.sqrt(3)/2],[-1,0],[1/2,np.sqrt(3)/2]]

with open('INIT','r') as f1:
    f1.readline()
    data = f1.readlines()
    for line in data:  
        INIT_data.append([ int(line.split()[0]), 
                           float(line.split()[1]), 
                           float(line.split()[2]), 
                           float(line.split()[3]),
                           float(line.split()[4]) ])

with open('DATA','r') as f2:
    f2.readline()
    data = f2.readlines()
    for line in data:  
        Clusters.append([ int(line.split()[0]), 
                          float(line.split()[1]), 
                          float(line.split()[2]), 
                          float(line.split()[3]),
                          float(line.split()[4]) ])

with open('PES','r') as f3:
    data = f3.readlines()
    for line in data:   
        PES.append([ str(line.split()[0]), 
                     float(line.split()[1]), 
                     float(line.split()[2]), 
                     float(line.split()[3]) ])

# METROPOLIS LOOP  

np.random.seed()  # seed for random number generator
start_time_MC   = timeit.default_timer()
OUTPUT_data  = copy.deepcopy(INIT_data)
PES_copy     = copy.deepcopy(PES)

with open('LOG', 'w') as f5:
    f5.write('%s\n' % ('**********LOG info**********'))

    
    
LCG = 0
for cluster in OUTPUT_data:
    if LCG <= cluster[1]:
        LCG =  cluster[1]
        

allPossibleLocation = []

# generating all the grid point and store in allPossibleLocation
for k in range(0, len(PES_copy)): 
    for i in range(0, param.num_primcellInOneDirection):
        for j in range(0, param.num_primcellInOneDirection):
            yIncreaseFactor =  param.primcell_b *j
            y = PES[k][2] + yIncreaseFactor * np.sqrt(3)/2
            x = PES[k][1] + param.primcell_a * i -  yIncreaseFactor* 1 /2
            allPossibleLocation.append([x, y])                                     
        
with open('metropolis','w') as f4:
    for step in range(Metro_Max+1):
        #check overlap
        
        overlapNumberCount = 0
        overlap = False
        indexListAll = []
        overlapCheck = True # overlapCheck is to check whether not there is more overlap. It is true as long as ovlp is true. It is false once ovlp is false.
        while (overlapCheck and (overlapNumberCount < param. LimitForOverlap)):
            
            OUTPUT_data, ovlp, index_list, LCG = overlap_check(Clusters, OUTPUT_data, LCG)
            totalAtoms1 = 0
            for tar in OUTPUT_data:
                totalAtoms1 = totalAtoms1 + tar[0]
           # print(str(overlapNumberCount)+ ": " + str(totalAtoms1))
            if (ovlp == True):
                overlap = True
            overlapCheck = ovlp            
            Ncluster = len(OUTPUT_data)
            if (ovlp == True and step == 0):
                raise ValueError('Overlapping clusters found in the initial setup!')
            elif (ovlp == False and step == 0):
                print('No overlapping clusters found in the initial setp!')
                break
                
            for index in index_list:
                indexListAll.append(index)

            overlapNumberCount += 1

        if (overlap == True and step > 0):
            with open('LOG', 'a') as f5:
                f5.write('%5s  %12i\n' % ('step =',step))
                f5.write('%27s \n' %  ('**********OVERLAP**********'))
                f5.write('%27s \n' %  ('Overlapping clusters found!'))
                for lst in indexListAll:
                    f5.write('%s' % (lst))
                f5.write('\n')



        # writing output  
        if ( (step % write_step) == 0 ): 
            f4.write('%5s  %10i %16s %3i \n' % ('step =',step,'numclusters =',Ncluster))
            f4.write('%4s  %14s  %14s  %14s  %14s\n' %  ('Pt', 'R', 'X','Y','E'))
            for i in range(Ncluster):
                f4.write('%3i  %16.8f  %16.8f  %16.8f  %16.8f\n' % 
                        (OUTPUT_data[i][0], OUTPUT_data[i][1], OUTPUT_data[i][2], OUTPUT_data[i][3], OUTPUT_data[i][4]))
            f4.write('\n')
        # choosing a cluster  
        indx = int(np.random.rand() * Ncluster) 
        num_atm_temp  = OUTPUT_data[indx][0] 
        R_temp        = OUTPUT_data[indx][1]
        X_temp        = OUTPUT_data[indx][2]
        Y_temp        = OUTPUT_data[indx][3]
        E_temp        = OUTPUT_data[indx][4]
        # choosing step size  
        if ( num_atm_temp == 1 ): #For monomer case.  
            pxy1 = int( 4.0*(2.0*np.random.rand() - 1.0) )   # -4 < px < 4   , scalar for the vector below
            pxy2 = int( 4.0*(2.0*np.random.rand() - 1.0) )   # -4 < px < 4   , scalar for the vector below
            pxy3 = int( 4.0*(2.0*np.random.rand() - 1.0) )   # -4 < px < 4   , scalar for the vector below
            #py = int( 4.0*(2.0*np.random.rand() - 1.0) )   # -4 < py < 4   
        else:  #If part of cluster we want bigger step, otherwise never breaks. We ensure that 
            pxy1 = math.ceil( ( 2.0*R_temp)  * ( 2.0*np.random.rand() - 1.0 ) ) # the step-size is cluster size dependent
            pxy2 = math.ceil( ( 2.0*R_temp)  * ( 2.0*np.random.rand() - 1.0 ) ) # the step-size is cluster size dependent
            pxy3 = math.ceil( ( 2.0*R_temp)  * ( 2.0*np.random.rand() - 1.0 ) ) # the step-size is cluster size dependent
            #py = math.ceil( ( 2.0*R_temp)  * ( 2.0*np.random.rand() - 1.0 ) ) 
            # we don't need two random px and py since we are scaling the direction, if we have two, it will not make any sense
            
        #this is where the new X and Y coords are generated- this is probably where to add the scaling vectors. 
        # we only need three directions to do the trick
        direction_choice1 = choose_dir[0] #this is choosing the direction 
        direction_choice2 = choose_dir[1] 
        direction_choice3 = choose_dir[2] 
    
        X_new = X_temp + direction_choice1[0] *pxy1 * param.xstep_max + direction_choice2[0] *pxy2 * param.xstep_max  + direction_choice3[0] * pxy3 * param.xstep_max#choosing the direction of the step
        Y_new = Y_temp + direction_choice1[1] *pxy1 * param.ystep_max + direction_choice2[1] *pxy2 * param.ystep_max + direction_choice3[1] * pxy3 * param.ystep_max

    
        x_find = X_finder(Y_new, param.primcell_a, param.maxx)       #to make the PBC easier (I think) 

        if ( X_new == X_temp and Y_new == Y_temp ): #no change
            #print("continue!!!")
            continue
            
        # Periodic Boundry Conditions - edited for Al2O3 lattice 
        #do Y PBC first
        if Y_new < param.miny:
            Y_new = param.maxy - param.miny + Y_new
            X_new = (X_new - x_find[0]) + X_finder(Y_new, param.primcell_a,param.maxx)[0]
        elif Y_new > param.maxy:
            Y_new = param.miny - param.maxy + Y_new
            X_new = (X_new - x_find[0]) + X_finder(Y_new, param.primcell_a,param.maxx)[0]

        # X BOUNDARY- need to account for the step size. 
        if X_new < X_finder(Y_new, param.primcell_a, param.maxx)[0]:
            X_new = X_new + X_finder(Y_new,param.primcell_a, param.maxx) [1] -  X_finder(Y_new,param.primcell_a, param.maxx) [0]
        elif X_new > X_finder(Y_new, param.primcell_a, param.maxx)[1]:
            X_new = X_new + X_finder(Y_new,param.primcell_a, param.maxx) [0] -  X_finder(Y_new,param.primcell_a, param.maxx) [1]

    
        # check for new clusters  
        inside_cluster = False
        for i in range (Ncluster):
            # if it is inside a cluster
            temp_dist = distance(X_new, Y_new, OUTPUT_data[i][2], OUTPUT_data[i][3])
            temp_r    = OUTPUT_data[i][1] + param.Ratom
            if ( (temp_dist <  temp_r) and (i !=indx) ):   
                Eold = E_temp + OUTPUT_data[i][4]       # total Eold
                if ( OUTPUT_data[indx][0] == 1): 
                    Enew_rmv = 0.0
                    Enew_add, Rnew_add = Cluster_finder(Clusters, OUTPUT_data, i, i, 'add')
                elif ( OUTPUT_data[indx][0] == 2): # Pt2 is special?
                    # CLUSTER FINDER FOR EADD
                    Enew_add, Rnew_add = Cluster_finder(Clusters, OUTPUT_data, i, i, 'add')
                    Enew_rmv = PES_finder(X_new, Y_new, PES_copy, N_mesh)
                    Rnew_rmv = param.Ratom 
                # choose a new cluster
                else:
                    Enew_add, Enew_rmv, Rnew_add, Rnew_rmv = Cluster_finder(Clusters, OUTPUT_data, indx, i, 'both')
                Enew = Enew_add + Enew_rmv
                if ( (Enew-Eold < 0.0 and E_temp != OUTPUT_data[i][4]) or 
                     (E_temp == OUTPUT_data[i][4] and np.exp( -beta*abs(Enew-Eold)) / (np.pi*OUTPUT_data[i][0]**2)  > np.random.rand()) or 
                     ((E_temp != OUTPUT_data[i][4] and np.exp( -beta*(Enew-Eold) ) > np.random.rand()))  ):  # accept the move
                # modify old clusters 
                    OUTPUT_data[i][0] = OUTPUT_data[i][0] + 1     # Natom
                    OUTPUT_data[i][1] = Rnew_add                     # R
                    OUTPUT_data[i][4] = Enew_add                     # E
                    if ( OUTPUT_data[indx][0] == 1 ):        #Natom = 1
                        # DELETE A CLUSTER (single atom)
                        del (OUTPUT_data[indx])
                        Ncluster -= 1
                    else:
                        OUTPUT_data[indx][0] = OUTPUT_data[indx][0] - 1     # Natom
                        OUTPUT_data[indx][1] = Rnew_rmv                        # R
                        OUTPUT_data[indx][4] = Enew_rmv                        # E

                else:
                    with open('LOG', 'a') as f5:
                        f5.write('%4s  %10i\n' % ('step=',step))
                        f5.write('%s \n' %  ('**********MOVE**********'))
                        f5.write('%s \n' %  ('METROPOLIS MOVE TO A NEW CLUSTER IS NOT FAVORABLE!'))
                inside_cluster = True
                break  
 
        if( inside_cluster == False ):  # we get a new single atom meaing that a new cluster forms      
            #print('inside cluster = false')
            if ( OUTPUT_data[indx][0] == 1):
                Enew_rmv = 0.0
            elif ( OUTPUT_data[indx][0] == 2):
                Enew_rmv = PES_finder(X_temp, Y_temp, PES_copy, N_mesh)
                Rnew_rmv = param.Ratom
            # choose a new cluster 
       
            else:
                Enew_rmv, Rnew_rmv = Cluster_finder(Clusters, OUTPUT_data, indx, indx, 'remove')
            E_atom = PES_finder(X_new, Y_new, PES_copy, N_mesh)
            Enew   = E_atom + Enew_rmv 
            Eold   = OUTPUT_data[indx][4]

            if ( (Enew-Eold < 0.0) or ( np.exp( -beta*(Enew-Eold) ) > np.random.rand()) ):  # accept the move
                with open('LOG', 'a') as f5:
                    f5.write('%5s  %12i\n' % ('step =',step))
                    f5.write('%s \n' %  ('**********MOVE**********'))
                    f5.write('%s \n' %  ('SINGLE ATOM MOVE OUTSIDE OF A CLUSTER ACCEPTED!'))
                # modify old clusters  
                if ( OUTPUT_data[indx][0] == 1 ):     # Note for a single atom move on surface
                    del (OUTPUT_data[indx])
                    Ncluster -= 1
                else:
                    OUTPUT_data[indx][0] = OUTPUT_data[indx][0] - 1     # Natom
                    OUTPUT_data[indx][1] = Rnew_rmv                        # R
                    OUTPUT_data[indx][4] = Enew_rmv                        # E
                OUTPUT_data.append([1, param.Ratom, X_new, Y_new, E_atom])             
                Ncluster += 1      
            else:
                with open('LOG', 'a') as f5:
                    f5.write('%5s  %12i\n' % ('step =',step))
                    f5.write('%s \n' %  ('**********MOVE**********'))
                    f5.write('%s \n' %  ('SINGLE ATOM MOVE OUTSIDE OF THE CLUSTER IS NOT ACCEPTED!'))

with open('LOG', 'a') as f5:
    f5.write('\n')
    f5.write('%s \n' %  ('**********CALCULATION IS DONE**********'))
    f5.write('%13s  %12i\n' % ('Total steps =',step))
    f5.write('%30s %4i \n' %  ('Final number of clusters =', Ncluster))
    f5.write('%30s %16.8f \n' %  ('Total MC time in seconds =', timeit.default_timer() - start_time_MC))
    f5.write('%5s \n' %  ('DONE!'))
    f5.write('%s \n' %  ('***************************************'))


if param.SinteringResultPlot: 
    # plot setting

    xcoords = []
    ycoords = []
    Rcoords = []
    types   = []
    for i in range(len(OUTPUT_data)):
        xcoords.append(OUTPUT_data[i][2])
        ycoords.append(OUTPUT_data[i][3])
        Rcoords.append(30*int(OUTPUT_data[i][1])**2)
        types.append(int(OUTPUT_data[i][0]))

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