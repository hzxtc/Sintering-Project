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

def overlap_check(Clusters, OUTPUT_data):
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
                    ovlp = True
                    idx.append(i)
                    idx.append(j)
                    idx_pair.append([i,j])
                    #merge
                    numnew = OUTPUT_data[i][0] + OUTPUT_data[j][0]
                    for k in range(len(Clusters)):
                        if (Clusters[k][0] == numnew):
                            Rnew = Clusters[k][1]
                            Enew = Clusters[k][4]
                    if (OUTPUT_data[i][0] > OUTPUT_data[j][0]):
                        xnew = OUTPUT_data[i][2]
                        ynew = OUTPUT_data[i][3]
                    else:
                        xnew = OUTPUT_data[j][2]
                        ynew = OUTPUT_data[j][3]
                    OUTPUT_data.append([numnew,Rnew,xnew,ynew,Enew])
                    break
    #removing duplicates from the list
    idx = list(dict.fromkeys(idx))
    idx = list(dict.fromkeys(idx))
    
    for i in sorted(idx, reverse=True):
        del (OUTPUT_data[i])

    return OUTPUT_data, ovlp, idx_pair

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

with open('metropolis','w') as f4:
    for step in range(Metro_Max+1):
        #check overlap
        numberOfOverlap = 0
        overlap = True
        
        while numberOfOverlap < param.LimitForOverlap and overlap: # while not reach LimitForOverlap and still overlap
            OUTPUT_data, ovlp, index_list = overlap_check(Clusters, OUTPUT_data)
            Ncluster = len(OUTPUT_data)
            overlap = ovlp
            numberOfOverlap = numberOfOverlap + 1
            if (ovlp == True and step == 0):
                print(index_list)
                raise ValueError('Overlapping clusters found in the initial setup!')
            elif (ovlp == False and step == 0):
                print('No overlapping clusters found in the initial setp!')
            elif (ovlp == True and step > 0):
                with open('LOG', 'a') as f5:
                    f5.write('%5s  %12i\n' % ('step =',step))
                    f5.write('%27s \n' %  ('**********OVERLAP**********'))
                    f5.write('%27s \n' %  ('Overlapping clusters found!'))
                    for lst in index_list:
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