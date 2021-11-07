#!/usr/bin/env python

'''
Sintering simulation via Ostwald Ripening of metalic clusters deposited on the surface 
based on NVT ensemble and metropolis moves.

Borna Zandkarimi 2020
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

# current_LCG is current largest cluster radius
# coords is OUTPUT_DATA
def boundaryOverlapCheck(coords, current_x, current_y, current_radius, current_LCG ):
    overlap = False
    targetPosition = -1

   #print("maxx: " + str(param.maxx) + " maxy: " + str(param.maxy))
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
    if current_x + current_radius >=  param.maxx:
        isRight = True
        overlapPotentialType.append( "right")
    if current_x - current_radius <= param.minx:
        isLeft = True
        overlapPotentialType.append("left" )
    if current_y + current_radius>= param.maxy:
        isUp = True
        overlapPotentialType.append("up")
    if current_y - current_radius<= param.miny:
        isDown = True
        overlapPotentialType.append("down")
    
    # no boundary overlap case
    if overlapPotentialType == [] :
        return overlap, targetPosition, current_LCG
    
    # newClusterList[][0] = x
    # newClusterList[][1] = y
    # newClusterList[][2] = radius
    # newClusterList[][3] = types
    # newClusterList[][4] = position
    newClusterList = []
    
    #print(overlapPotentialType)
    for types in overlapPotentialType:
        if types == "right" :
            # delta is the longest length that other potential overlap cluster can exist
            delta = current_x + current_radius - param.maxx + 2 * current_LCG
            #counting position
            n = 0
            for cluster in coords:
                if cluster[2] <=  param.minx + delta:
                    newClusterList.append([cluster[2] + param.maxx, cluster[3], cluster[1], cluster[0], n])       
                n += 1
        if types == "left" :
            # this delt should be a negative number
            delta = current_x - current_radius + param.minx - 2 * current_LCG
            n = 0
            for cluster in coords:
                if cluster[2] >= param.maxx + delta:
                    newClusterList.append([cluster[2] - param.maxx, cluster[3], cluster[1], cluster[0], n])
                n += 1
        if types == "up" :
            delta = current_y + current_radius - param.maxy + 2 * current_LCG
            n = 0
            for cluster in coords:
                if cluster[3] <= param.miny + delta:
                    newClusterList.append([cluster[2], cluster[3] +param.maxy, cluster[1], cluster[0], n])
                n += 1
        if types == "down" :
            delta = current_y - current_radius + param.miny - 2 * current_LCG
            n = 0
            for cluster in coords:
                if cluster[3] >= param.maxy +delta:
                    newClusterList.append([cluster[2], cluster[3] - param.maxy, cluster[1],cluster[0], n])
                n += 1
    
    # check for corner cases
    if isRight and isUp :
        n = 0
        for cluster in coords:
                if (cluster[2] <= current_LCG) and (cluster[3] <= current_LCG):
                    newClusterList.append([cluster[2] + param.maxx, cluster[3] + param.maxy, cluster[1], cluster[0], n])
                n += 1
    if isRight and isDown :
        n = 0
        for cluster in coords:
                if (cluster[2] <= current_LCG) and (cluster[3] >= param.maxy - current_LCG):
                    newClusterList.append([cluster[2] + param.maxx, cluster[3] - param.maxy, cluster[1], cluster[0], n])
                n += 1
    if isLeft and isUp :
        n = 0
        for cluster in coords:
                if (cluster[2] >= param.maxx - current_LCG) and (cluster[3] <= current_LCG):
                    newClusterList.append([cluster[2] - param.maxx, cluster[3] + param.maxy, cluster[1], cluster[0], n])
                n += 1
    if isLeft and isDown :
        n = 0
        for cluster in coords:
                if (cluster[2] >= param.maxx - current_LCG) and (cluster[3] >= param.maxy - current_LCG):
                    newClusterList.append([cluster[2] - param.maxx, cluster[3] - param.maxy, cluster[1], cluster[0], n])
                n += 1
    
    OVERLAPNUMBERCOUNT = 0
    # checking overlap
    for testCluster in newClusterList:
        #  print("radius: " + str(current_radius))
        #  print("radius: " + str(testCluster[2]))
        if distance(current_x,current_y, testCluster[0],testCluster[1]) <= current_radius + testCluster[2]:
            #print(str(current_x) + " " + str(current_y) + " " +str( testCluster[0]) + " " + str(testCluster[1]) + " " + str(current_radius) + " " + str(testCluster[2]))
            targetPosition = testCluster[4]
            overlap = True
            OVERLAPNUMBERCOUNT += 1 
            break
               
    #print(overlap)
    # print(len(newClusterList))
    #  print(OVERLAPNUMBERCOUNT) 
    # return whether overlap, position of overlap cluster
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
                    if (OUTPUT_data[i][0] > OUTPUT_data[j][0]):
                        xnew = OUTPUT_data[i][2]
                        ynew = OUTPUT_data[i][3]
                    else:
                        xnew = OUTPUT_data[j][2]
                        ynew = OUTPUT_data[j][3]
                    OUTPUT_data.append([numnew,Rnew,xnew,ynew,Enew])
                    break
            # overlap boundary check
            ovlp1 = False
            ovlp1, overlapPosition, current_LCG = boundaryOverlapCheck(OUTPUT_data, OUTPUT_data[i][2],  OUTPUT_data[i][3], OUTPUT_data[i][1], LCG)
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
                print("boundary overlap")
                print("LCG: " + str(LCG))
                print("x1: " + str(OUTPUT_data[i][2]) + " y1: " + str(OUTPUT_data[i][3])+ " r1: " + str(OUTPUT_data[i][1]) + " x2: " + str(OUTPUT_data[j1][2]) + " y2: " + str(OUTPUT_data[j1][3])+ " r2: " + str(OUTPUT_data[j1][1]))
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
                if (OUTPUT_data[i][0] > OUTPUT_data[j1][0]):
                    xnew = OUTPUT_data[i][2]
                    ynew = OUTPUT_data[i][3]
                else:
                    xnew = OUTPUT_data[j1][2]
                    ynew = OUTPUT_data[j1][3]
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
    for ii in range(N_mesh): 
        dx =  np.abs( param.xstep_max*( math.ceil(X_new/param.xstep_max) % N_meshx ) - PES_copy[ii][1] )
        dy =  np.abs( param.ystep_max*( math.ceil(Y_new/param.ystep_max) % N_meshy ) - PES_copy[ii][2] )
        if ( dx  < 0.1 and dy  < 0.1 ):    
            # PES found 
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
numprimcell   = param.numprimecellFactor*param.num_clust + param.num_single_atom
xdups = ydups = (int(numprimcell**0.5) + 2)

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

# Boltzmann population information, specific to TiO2
T=700
hartree_to_ev = 27.2114
k_J           = 1.38064852e-23
k_ev          = 8.6173303e-5
k_hartree     = k_ev/hartree_to_ev
kT            = k_hartree*T
beta          = 1.0/kT

Pt2=[1.09778751,1.09970170,1.10494276,1.10750419,1.10750419]
Pt2_index=[0,1,2,3,4]
Pt2_boltz_weight=[math.exp(-Pt2[i]/kT) for i in range(0,len(Pt2))]
Q2=sum(Pt2_boltz_weight)
Pt2_prob=[Pt2_boltz_weight[i]/Q2 for i in range(0,len(Pt2_boltz_weight))]
Pt2_data=[[3.0218605,1.09778751],[3.0283540,1.09970170],[3.0056904 ,1.10494276],[2.9766112 ,1.10750419],
          [2.8958918 ,1.10750419]]

Pt3=[0.90056183,0.91199504,0.91510903]
Pt3_index=[0,1,2]
Pt3_boltz_weight=[math.exp(-Pt3[i]/kT) for i in range(0,len(Pt3))]
Q3=sum(Pt3_boltz_weight)
Pt3_prob=[Pt3_boltz_weight[i]/Q3 for i in range(0,len(Pt3_boltz_weight))]
Pt3_data=[[3.0040290,0.90056183],[2.9851504,0.91199504],[2.8275859,0.91510903]]

Pt4=[0.72210192,0.72929952,0.73153563,0.73305648,0.73571588,0.73606268]
Pt4_index=[0,1,2,3,4,5]
Pt4_boltz_weight=[math.exp(-Pt4[i]/kT) for i in range(0,len(Pt4))]
Q4=sum(Pt4_boltz_weight)
Pt4_prob=[Pt4_boltz_weight[i]/Q4 for i in range(0,len(Pt4_boltz_weight))]
Pt4_data=[[ 3.6684145,0.72210192],[3.0623059,0.72929952],[3.6433934,0.73153563],[3.2190594,0.73305648],
         [3.0076024,0.73571588],[3.0904500,0.73606268]]

Pt5=[0.5417517,0.54831806]
Pt5_index=[0,1]
Pt5_boltz_weight=[math.exp(-Pt5[i]/kT) for i in range(0,len(Pt5))]
Q5=sum(Pt5_boltz_weight)
Pt5_prob=[Pt5_boltz_weight[i]/Q5 for i in range(0,len(Pt5_boltz_weight))]
Pt5_data=[[3.7748643 ,0.54175170],[3.6812070,0.54831806]]

Pt6=[0.36589531,0.37585195,0.37773049,0.37906140,0.38048965]
Pt6_index=[0,1,2,3,4]
Pt6_boltz_weight=[math.exp(-Pt6[i]/kT) for i in range(0,len(Pt6))]
Q6=sum(Pt6_boltz_weight)
Pt6_prob=[Pt6_boltz_weight[i]/Q6 for i in range(0,len(Pt6_boltz_weight))]
Pt6_data=[[4.1423044, 0.36589531],[4.2300866 ,0.37585195],[4.1883873 ,0.37773049],[3.9624653,0.37906140],
     [4.1217715,0.38048965]]

Pt7=[ 0.17865857,0.18249841,0.18304948,0.18388096,0.18630240,0.19132231,0.19230703]
Pt7_index=[0,1,2,3,4,5,6]
Pt7_boltz_weight=[math.exp(-Pt7[i]/kT) for i in range(0,len(Pt7))]
Q7=sum(Pt7_boltz_weight)
Pt7_prob=[Pt7_boltz_weight[i]/Q7 for i in range(0,len(Pt7_boltz_weight))]
Pt7_data=[ [4.2377875,0.17865857],[4.2400853,0.18249841],[4.9129259,0.18304948],[4.7864910,0.18388096],
          [4.4981330,0.18630240],[4.1891448,0.19132231],[4.9218588,0.19230703]]

Pt8=[0.00000000,0.00707852,0.01104738,0.01222103,0.01320434]
Pt8_index=[0,1,2,3,4]
Pt8_boltz_weight=[math.exp(-Pt8[i]/kT) for i in range(0,len(Pt8))]
Q8=sum(Pt8_boltz_weight)
Pt8_prob=[Pt8_boltz_weight[i]/Q8 for i in range(0,len(Pt8_boltz_weight))]
Pt8_data=[ [4.6229707,0.00000000],[4.5883157 ,0.00707852],[4.6106749,0.01104738],[4.2228037,0.01222103],
          [4.4044703,0.01320434]]

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
# METROPOLIS LOOP  

np.random.seed()  # seed for random number generator
start_time_MC   = timeit.default_timer()
OUTPUT_data  = copy.deepcopy(INIT_data)# ??? I change INIT_data to INIT
PES_copy     = copy.deepcopy(PES)

with open('LOG', 'w') as f5:
    f5.write('%s\n' % ('**********LOG info**********'))

    
    
LCG = 0
for cluster in OUTPUT_data:
    if LCG <= cluster[1]:
        LCG =  cluster[1]

with open('metropolis','w') as f4:
    for step in range(Metro_Max+1):
        totalAtoms = 0
        for tar in OUTPUT_data:
            totalAtoms = totalAtoms + tar[0]
       # print(str(step)+ ": " + str(totalAtoms))
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
                f4.write('%3i  %16.8f  %16.8f  %16.8f  %16.8f\n' % (OUTPUT_data[i][0], OUTPUT_data[i][1], OUTPUT_data[i][2], OUTPUT_data[i][3], OUTPUT_data[i][4]))
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
            px = int( 4.0*(2.0*np.random.rand() - 1.0) )   # -4 < px < 4   
            py = int( 4.0*(2.0*np.random.rand() - 1.0) )   # -4 < py < 4   
        else:  #If part of cluster we want bigger step, otherwise never breaks. We ensure that 
            px = math.ceil( ( 2.0*R_temp)  * ( 2.0*np.random.rand() - 1.0 ) ) # the step-size is cluster size dependent
            py = math.ceil( ( 2.0*R_temp)  * ( 2.0*np.random.rand() - 1.0 ) ) 
        X_new = X_temp + px * param.xstep_max
        Y_new = Y_temp + py * param.ystep_max
        if ( X_new == X_temp and Y_new == Y_temp ):
            continue
        # Periodic Boundry Conditions  
        if ( X_new > param.maxx ):   
            X_new = param.minx - param.maxx + X_new
        elif ( X_new < param.minx ):   
            X_new = param.maxx - param.minx + X_new
        if ( Y_new > param.maxy ):  
            Y_new = param.miny - param.maxy + Y_new
        elif ( Y_new < param.miny ): 
            Y_new = param.maxy - param.miny + Y_new
        # check for new clusters  
        inside_cluster = False
        for i in range (Ncluster):
            # if it is inside a cluster. the cluster can be just one atom
            temp_dist = distance(X_new, Y_new, OUTPUT_data[i][2], OUTPUT_data[i][3])
            temp_r    = OUTPUT_data[i][1] + param.Ratom
            if ( (temp_dist <  temp_r) and (i !=indx) ):# i != indx means that it cannot be itself.   
                Eold = E_temp + OUTPUT_data[i][4]       # total Eold
                if ( OUTPUT_data[indx][0] == 1):# to check whether it is an atom or not
                    Enew_rmv = 0.0
                    Enew_add, Rnew_add = Cluster_finder(Clusters, OUTPUT_data, i, i, 'add')
                elif ( OUTPUT_data[indx][0] == 2):# ??? why the case of two atoms is special?
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
    for i in range(len(OUTPUT_data)):
        xcoords.append(OUTPUT_data[i][2])
        ycoords.append(OUTPUT_data[i][3])
        Rcoords.append(30*int(OUTPUT_data[i][1])**2)
        types.append(int(OUTPUT_data[i][0]))

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
