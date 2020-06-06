# -*- coding: utf-8 -*-
"""
Created on Tue May 06 23:50:47 2014

@author: Harpal
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Current structure of program
# 1. Loads graph that has been previously constructed by makeCrystalGraph.py
# 2. Randomises initial structure
# 3. Runs IsingModel()
# 4. User can then choose to run additional functions, e.g.
#    graphSRO()
#    isingModel()

# If user doesn't have 'eta' model then comment out all lines containg 'eta'
# Should put 'eta' code in try/except blocks really... (or whatever is recommended)

import networkx as nx
import math
import random
import os
import numpy

#import cProfile
#import pstats
#import StringIO
#from eta import ETA

##############################################################################

####################
# GLOBAL VARIABLES #
####################

G=''
H=''
J=''
size ='32x32x5'
#filename = r'input/%s-quantum-well-graph-Ganchenkova.gpickle'%size
#filename = r'make_crystal/v2_crystal/%s-graph-Ganchenkova.gpickle'%size
filename = r'make_crystal/v1_crystal/%s-quantum-well-graph-Ganchenkova.gpickle'%size
lenJ=''
# for graphing SRO
SROx = []
SROy1 = []
# Scientific constants
kB = 8.6173324*(10e-5) # eV
# Simulation parameters
T = 773 # K
kT = kB*T
indiumComposition = 0.5
steps = 5000000
statsFrequency = 10000  # generate statistics every statsFrequency steps
saveFrequency = int(steps)/10
##############################################################################
edgesDict={}

def getNeighbours(i,k,P=None):
    """
    input your node and the weight of the neighbours you wanna find, i and
    k. Retrieve a list of the nodes of the kth neigbour.
    """

    global edgesDict
    key=(i,k)
    if key in edgesDict:
        return edgesDict[key]
    else:
        global H
        if P:
            H=P
        edges = H.edges(i,k)   
        edges_of_weight = []
        neighbours = []
        for edge in edges:
            if edge[2]['weight'] == k:
                edges_of_weight.append(edge)
        for anElement in edges_of_weight:
            neighbours.append(anElement[1])
        edgesDict[key]=neighbours
        return neighbours
       
def attemptAtomSwap():
    """
    See paper for definitions:
    
    N_A_12_In = The number of indium atoms within the 1st and 2nd shell of the 
    randomly selected indium atom
    
    N_A_34_Ga = The number of indium atoms within the 3rd and 4th shell of the
    randomly selected Gallium atom
    
    AB12_In = n(A,B,1) + n(A,B,2) for the randomly selected Indium atom
    
    Ein = The energy of the initial configuration
    
    Efi = The energy of the final configuration.
    
    Retrieve E12 and E34 from the functions E_aa, E_ac, E_ac2_and E_cc.

    E_xy selects the appropriate binding energies as a function of indiumComposition

    E_xy were produced by plotting trendlines from experimental data and interpolating
    

    What I do here is:
    
    Pick an In and a Ga atom at random (not at the surface)
    
    Check the number of In atoms around the In atom N_A_In_12
    
    Find the total number of such pairs at "4" distances e.g. nAB12_In+nAB34_In
    
    Multiply that by the appropriate energy of those configurations to find the
    equilibrium energy of the local environment at the BEFORE the swap.
    
    Do the same for the new position that the Indium would be AFTER the swap.So
    treat the neighbours around the Gallium atoms position AS IF it 
    were an Indium atom.
    
    Hence we +1 to account for the fact that 1 In atom in the crystal 
    means that get_N_A_12 and get_N_A_34 returns zero.
    
    Plug the value into the ProbabilityOfAcceptance and the rest is history.       
    """
    randomInNode=None
    randomGaNode=None
    # pick an In atom at random NOT in the surface

    #while not randomInNode or not randomGaNode:
    #    randomNodes = random.sample(H.node,2)
    #    for r in randomNodes:
    #        if H.node[r]['species']=='Ga':
    #            randomGaNode=r
    #        else:
    #            randomInNode=r

    while not randomInNode or not randomGaNode:
        r=J[random.randint(0,lenJ-1)]
        if H.node[r]['species']=='Ga':
            randomGaNode=r
        else:
            randomInNode=r
       
#+1 includes counts the Indium atom that we are looking for neighbours for in the count
#Otherwise this could count ZERO In-In pairs when get_N_A_XXX() returns zero In neighbours
#Must +1for ALL as well because there should be an In atom in the position of the Ga post swap
    
    NumberOfInAtoms_In_site = 1+get_N_A_1ac(randomInNode)+get_N_A_1aa(randomInNode)+get_N_A_2ac(randomInNode)+get_N_A_2cc(randomInNode)
    NumberOfInAtoms_Ga_site = 1+get_N_A_1ac(randomGaNode)+get_N_A_1aa(randomGaNode)+get_N_A_2ac(randomGaNode)+get_N_A_2cc(randomGaNode)     

    nconfig_InIn1ac_In = NumberOfInAtoms_In_site*3
    #nconfig_InIn1aa_In = NumberOfInAtoms_In_site*3
    #nconfig_InIn2ac_In = NumberOfInAtoms_In_site*3
    nconfig_InIn2cc_In = NumberOfInAtoms_In_site
    nconfig_InIn1ac_Ga = NumberOfInAtoms_Ga_site*3
    #nconfig_InIn1aa_Ga = NumberOfInAtoms_Ga_site*3
    #nconfig_InIn2ac_Ga = NumberOfInAtoms_Ga_site*3
    nconfig_InIn2cc_Ga = NumberOfInAtoms_Ga_site
    
    Eb_1ac_comp = Eb_1ac(indiumComposition)
    Eb_1aa_comp = Eb_1aa(indiumComposition)
    Eb_2ac_comp = Eb_2ac(indiumComposition)
    Eb_2cc_comp = Eb_2cc(indiumComposition)     
    
    #Ein = nconfig_InIn1ac_In*Eb_1ac_comp + nconfig_InIn1aa_In*Eb_1aa_comp + nconfig_InIn2ac_In*Eb_2ac_comp + nconfig_InIn2cc_In*Eb_2cc_comp
    #Efi = nconfig_InIn1ac_Ga*Eb_1ac_comp + nconfig_InIn1aa_Ga*Eb_1aa_comp + nconfig_InIn2ac_Ga*Eb_2ac_comp + nconfig_InIn2cc_Ga*Eb_2cc_comp
    Ein = nconfig_InIn1ac_In*Eb_1ac_comp + nconfig_InIn1ac_In*Eb_1aa_comp + nconfig_InIn1ac_In*Eb_2ac_comp + nconfig_InIn2cc_In*Eb_2cc_comp
    Efi = nconfig_InIn1ac_Ga*Eb_1ac_comp + nconfig_InIn1ac_Ga*Eb_1aa_comp + nconfig_InIn1ac_Ga*Eb_2ac_comp + nconfig_InIn2cc_Ga*Eb_2cc_comp
    probabilityOfAcceptance = math.exp(((Ein)-(Efi))/(kT))
    acceptSwap = False
    if Ein > Efi or random.random() < probabilityOfAcceptance:
        acceptSwap = True
        G.node[randomInNode]['species'],G.node[randomGaNode]['species'] = 'Ga','In'
    #DEBUG STRING ON NEXT LINE. UNCOMMENT TO PRINT INFO DURING CALCULATION.
    #print("InitGa: %i, FinalGa: %i, dH: %f, P: %f, Accepted: %r" % (initialNumberOfGaGaBonds,finalNumberOfGaGaBonds,deltaH,probabilityOfAcceptance,acceptSwap))
    
    return acceptSwap

def Eb_1ac(indiumComposition):
    E1 = 0.4945*(indiumComposition**4) - 0.727*(indiumComposition**3) + 0.3929*(indiumComposition**2) - 0.0925*indiumComposition - 0.072
    return E1

def Eb_1aa(indiumComposition):
    if indiumComposition == 0.1:
        E2=-0.065
    else:
        E2 = 0.0738*indiumComposition - 0.0771
    return E2

def Eb_2ac(indiumComposition):
    if indiumComposition == 0.1:
        E3 = 0.007
    if indiumComposition == 0.25:
        E3 = 0.008
    else:
        E3 = 0.0133*indiumComposition + 0.0052
    return E3
    
def Eb_2cc(indiumComposition):
    E4 = 0.1238*(indiumComposition**2)-0.1122*indiumComposition + 0.04
    return E4

def get_N_A_1ac(i):
    """
    input a node, return the number of In atoms of the 1ac neighbours(k=1)    
    """
    neighboursOfRandomNode = getNeighbours(i,1)
    numberOfIndiumNeighbours = 0
    for neighbour in neighboursOfRandomNode:
        if G.node[neighbour]['species'] == 'In':
            numberOfIndiumNeighbours += 1 
    return numberOfIndiumNeighbours

def get_N_A_1aa(i):
    """
    input a node, return the number of In atoms of the 1aa neighbours(k=2)    
    """
    neighboursOfRandomNode = getNeighbours(i,2)
    numberOfIndiumNeighbours = 0
    for neighbour in neighboursOfRandomNode:
        if G.node[neighbour]['species'] == 'In':
            numberOfIndiumNeighbours += 1 
    return numberOfIndiumNeighbours
    
def get_N_A_2ac(i):
    """
    input a node, return the number of In atoms of the 2ac neighbours(k=3)    
    """
    neighboursOfRandomNode = getNeighbours(i,3)
    numberOfIndiumNeighbours = 0
    for neighbour in neighboursOfRandomNode:
        if G.node[neighbour]['species'] == 'In':
            numberOfIndiumNeighbours += 1 
    return numberOfIndiumNeighbours

def get_N_A_2cc(i):
    """
    input a node, return the number of In atoms of the 2cc neighbours(k=4)     
    """
    neighboursOfRandomNode = getNeighbours(i,4)
    numberOfIndiumNeighbours = 0
    for neighbour in neighboursOfRandomNode:
        if G.node[neighbour]['species'] == 'In':
            numberOfIndiumNeighbours += 1 
    return numberOfIndiumNeighbours

def runIsingModel_SRO_1(numberOfSwaps):
    """ 
    Runs attemptAtomSwap() numberOfSwaps times.
    """
    ns=float(numberOfSwaps)
    
    #etaIsing = ETA(numberOfSwaps/statsFrequency) # set up a progress bar if 'eta' package is available
    print ("Running Ising Model to calculate SRO_1... Please wait.")    
    successfulSwaps = 0
    SwapsSoFar = 0
    for i in range(numberOfSwaps):
        SwapsSoFar += 1
        if attemptAtomSwap():
            successfulSwaps += 1
        if(i%statsFrequency==0):
            print float(i)/ns*100,'% Done'
            #etaIsing.print_status()
            SROx.append(i)
            SROy1.append(getSRO_1())
        if(i%saveFrequency==0 and i>2):
            #writeGraphToXYZ(size+'_xIn='+str(indiumComposition)+'_T='+str(T)+' K_SRO_1_'+str(SROy1[-1])+'SwapsSoFar'+str(SwapsSoFar)+'.xyz')
            pass
    print("\nOut of %i attempted atom swaps, %i were successful (ratio of %f)."%(numberOfSwaps,successfulSwaps,float(successfulSwaps)/numberOfSwaps))
    return successfulSwaps

def writeGraphToXYZ(myFilename):
    f = open(myFilename,'w')
    f.write(str(len(G))+'\n')
    f.write('Atoms. File created from networkx graph by '+os.path.basename(__file__)+'\n')
    for nodeindex in G:
        atom = G.node[nodeindex]
        f.write(str(atom['species'])+' '+str(atom['x'])+' '+str(atom['y'])+' '+str(atom['z'])+'\n')
    f.close()
    print("Graph exported as .xyz file.")

def getSRO_1():
    """Returns the Warren-Cowley short range order parameter for the first
    nearest neighbour (l=1). See Chan, Liu and Zunger paper for more
    information. DOI: 10.1103/PhysRevB.82.045112 """
    indiumCount = []
    global H
    for node in H:
        if H.node[node]['species'] == 'Ga':
            neighbours = getNeighbours(node,1,H)
            numberOfIndiumNeighbours = 0
            for neighbour in neighbours:
                if H.node[neighbour]['species'] == 'In':
                    numberOfIndiumNeighbours += 1
            indiumCount.append(numberOfIndiumNeighbours/6.0) # e.g. for indiumComposition=0.5, numberofIndiumNeighbours would be 3 on average for a random composition but might fluctuate from site to site, and therefore indiumCount would be appended by 0.5 on average
    indiumProbability = sum(indiumCount)/len(indiumCount) # this calculates that average!
    return (1 - (indiumProbability/indiumComposition))
    




def graphSRO_1():
    print ("Graphing SRO_1ac... Please wait")
    import matplotlib.pyplot as plt
    plt.plot(SROx,SROy1)
    plt.xlabel('Step')
    plt.ylabel('SRO_1_'+str(indiumComposition)+'_'+str(T)+'K')
    plt.show()

def printSRO_1():
    print('Steps:')
    print(SROx)
    print('SRO_1:')
    print(SROy1)

def main():
   #pr=cProfile.Profile()
    # this means it grabs the filename from first argument when you run the script from the command line
    # alternatively could change this to "filename = '64-64-10-graph.gpickle'" or similar
    #pr.enable()
    global G,H,lenJ,J
    # Loads our initial graph that was made from makeCrystalGraph.py
    print("Loading initial graph.")
    try:
        G = nx.read_gpickle(filename)
    except:
        try:
            G = nx.read_graphml(filename)
        except:
            raise Exception("Couldn't load graph! Generate a graph first using makeCrystalGraph.py")
    print("Graph loaded.")
    # Assuming our initial graph is pure GaN, we need to create a random initial
    # composition.
    J=[] #J will be our list that holds our nodes where surface == False 
    print 'size is',size.split('x')[-1]
    if size.split('x')[-1]=='1':
        for node in G:
            if G.node[node]['surface']!=False:
                J.append(node)
        H = G.subgraph(J) #H is the subgraph which contains only surface atoms.
    else:
        for node in G:
            if G.node[node]['surface']==False:
                J.append(node)
        H = G.subgraph(J) #H is the subgraph which contains only surface atoms.
    #This speeds up the code
    lenJ=len(J)

   
    ###############################################################################
    print("Randomising initial configuration to match an overall composition of InₓGa₁₋ₓN (x=%f)." % indiumComposition)
    totalNumberOfGaSites = 0
    currentNumberOfInAtoms = 0



    for node in H: # count up number of Ga sites in our graph
        if H.node[node]['species'] == 'Ga':
            totalNumberOfGaSites += 1
        if H.node[node]['species'] == 'In':
            currentNumberOfInAtoms += 1
    if currentNumberOfInAtoms != 0:
        raise Exception("Code assumes that there are no In atoms in the initial graph.")
    desiredNumberOfInAtoms = int(math.floor(indiumComposition*totalNumberOfGaSites))
    #etaRandomising = ETA(desiredNumberOfInAtoms)


    for i in range(desiredNumberOfInAtoms):
       #etaRandomising.print_status()
        randomNode=J[random.randint(0,lenJ-1)]# find a random node
        while H.node[randomNode]['species']!='Ga': # and ensure random node is Ga and not on the 'surface' of the quantum well
                randomNode=J[random.randint(0,lenJ-1)]
        H.node[randomNode]['species']='In' # change node to In
        currentNumberOfInAtoms += 1
    print("\nInitial composition randomised. Out of a total of %i Ga sites, %i are now occupied by In atoms for a overall composition of %f." % (totalNumberOfGaSites,currentNumberOfInAtoms,float(currentNumberOfInAtoms)/totalNumberOfGaSites))   
    ###############################################################################
    runIsingModel_SRO_1(steps)
    #graphSRO_1()
    printSRO_1()
    a = numpy.asarray([SROx,SROy1])
    numpy.savetxt('CSV_'+size+'_xIn='+str(indiumComposition)+'_T='+str(T)+' K_SRO_1'+'.csv', a, delimiter=",")
    #saveDoto()
    #writeGraphToXYZ(size+'_xIn='+str(indiumComposition)+'_T='+str(T)+' K_SRO_1_'+str(SROy1[-1])+'Final'+'.xyz')
    nx.write_gpickle(G,size+'_xIn='+str(indiumComposition)+'_T='+str(T)+' K_SRO_1_'+str(SROy1[-1])+'Final'+'.gpickle')
    #pr.disable()
    #s=StringIO.StringIO()
    #ps=pstats.Stats(pr,stream=s).sort_stats('cumulative')
    #ps.print_stats(0.05)
    #print s.getvalue()

if __name__=='__main__':
    main()
