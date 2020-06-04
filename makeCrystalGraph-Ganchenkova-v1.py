#!/usr/bin/env python
# -*- coding: utf-8 -*-
# M.K. Horton, m.horton11@imperial.ac.uk

import sys
import math
from ase import *
import networkx as nx
import numpy as np
from scipy import spatial
from collections import Counter

# constructs a GaN supercell from scratch as a networkx graph
# includes real-space co-ordinates for easier visualisation at a later date

# need to refactor this to avoid code reuse! and make it generally clearer what's being printed at the end

detectSurfaces = True
wrapZ = False # TRUE for a bulk graph (no surfaces, periodic in z), FALSE for a quantum well (surfaces, only periodic in x and y)
export = True

# radius = 3.183 # 3.3 for Z=12, 3.189 produces weird bug

radii =  [3.17902, 3.19, 4.49944, 5.176]  # 1ac, aa, 2ac, 2c, added 0.01 just in case [3.16902, 3.18, 4.48944, 5.166]
weights = [1,2,3,4]

print('Constructing supercell.')

na = 32
nb = 32
nc = 5

a = 3.18
c = 5.166
atoms = Atoms([Atom('Ga', (0,       0,    0)),
               Atom('Ga', (1/3., 1/3., 1/2.))]) # only looking at Ga atoms!
#               Atom('N',  (0,       0, 3/8.)),
#               Atom('N',  (1/3., 1/3., 7/8.))])
cell = [(a,             0, 0),
        (a/2, a*math.sqrt(3)/2, 0),
        (0,             0, c)]
atoms.set_cell(cell, scale_atoms=True)

qw = atoms * (na,nb,nc)

qwAbsolutePositions = qw.get_positions()
qwFractionalPositions = qw.get_scaled_positions()

G = nx.Graph()
for i in range(len(qw)):
    G.add_node(i)
    G.node[i]['species'] = 'Ga'
    G.node[i]['x'] = float(qwAbsolutePositions[i][0])
    G.node[i]['y'] = float(qwAbsolutePositions[i][1])
    G.node[i]['z'] = float(qwAbsolutePositions[i][2])
    G.node[i]['fx'] = float(qwFractionalPositions[i][0])
    G.node[i]['fy'] = float(qwFractionalPositions[i][1])
    G.node[i]['fz'] = float(qwFractionalPositions[i][2])
    if qwFractionalPositions[i][2] < (0.1/nc) or qwFractionalPositions[i][2] > (1-(0.6/nc)):
        if detectSurfaces:
            G.node[i]['surface'] = True
        else:
            G.node[i]['surface'] = False
    else:
        G.node[i]['surface'] = False

# Now for the cool stuff!
print('Constructing neighbour list.')
print('Counters show number of atoms with a given number of neighbours and are used for error checking.')

# documentation here:
# http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query_ball_point.html#scipy.spatial.cKDTree.query_ball_point

def shiftCell(axis,frac):
    # axis: 0 is x, 1 is y, 2 is z
    # frac: is a float showing how much the cell is shifted fractionally in that direction
    print('Constructing (shifting axis '+str(axis)+' by fraction '+str(frac)+') neighbour list.')
    for i in range(len(qwFractionalPositions)):
        qwFractionalPositions[i][axis] += frac # change fractional co-ordinate
        qwFractionalPositions[i][axis] -= math.floor(qwFractionalPositions[i][axis])
    qw.set_scaled_positions(qwFractionalPositions) # update the absolute co-ordinates with the new fractional co-ordinates

def findEdges(radius,weight,print_stats):
    qwAbsolutePositions = qw.get_positions()
    tree = spatial.cKDTree(qwAbsolutePositions)
    for node in G:
        x = qwAbsolutePositions[node][0]
        y = qwAbsolutePositions[node][1]
        z = qwAbsolutePositions[node][2]
        neighbours = tree.query_ball_point([x,y,z],radius) # where radius is next nearest distance; may be wrong!
        for neighbour in neighbours:
           if neighbour != node:
                if G.has_edge(node,neighbour) == False:
                    G.add_edge(node,neighbour,weight=weight)
    if print_stats:
        neighbourList = []
        for node in G:
            neighbourList.append(len(getEdgesWithWeight(node,weight)))
        print(Counter(neighbourList))

def getEdgesWithWeight(node,weight):
    edges = G.edges(node,weight)
    selected_edges = []
    for edge in edges:
        if edge[2]['weight'] == weight: # get the weight for each edge and see if it's the weight we want
            selected_edges.append(edge)
    return selected_edges


tree = spatial.cKDTree(qwAbsolutePositions)

def findNeighbours(radius,weight):
    findEdges(radius,weight,False) # not print stats
    shiftCell(0,0.5) # x
    findEdges(radius,weight,False)
    shiftCell(1,0.5) # y
    findEdges(radius,weight,False)
    shiftCell(0,0.5) # x again
    findEdges(radius,weight,False)

    if wrapZ:
        shiftCell(2,0.5) # z
        findEdges(radius,weight,False)
        shiftCell(0,0.5) # x again
        findEdges(radius,weight,False)
        shiftCell(1,0.5) # y again
        findEdges(radius,weight,True)

for i in range(len(radii)):
    print ('Finding neighbours for radius '+str(radii[i])+' and assigning weight '+str(weights[i])+'.')
    findNeighbours(radii[i],weights[i])

if export:
    if wrapZ:
        nx.write_gpickle(G,str(na)+'x'+str(nb)+'x'+str(nc)+'-bulk-graph-Ganchenkova.gpickle')
    else:
        nx.write_gpickle(G,str(na)+'x'+str(nb)+'x'+str(nc)+'-quantum-well-graph-Ganchenkova.gpickle')
    print ('Graph constructed and exported successfully.')
