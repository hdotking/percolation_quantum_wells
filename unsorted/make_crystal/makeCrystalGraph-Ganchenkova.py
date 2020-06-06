#!/usr/bin/env python
# -*- coding: utf-8 -*-
# M.K. Horton, m.horton11@imperial.ac.uk

import sys
import math
from ase import *
import networkx as nx
import numpy as np
from scipy import spatial

# constructs a GaN supercell from scratch as a networkx graph
# includes real-space co-ordinates for easier visualisation at a later date

# need to refactor this to avoid code reuse! and make it generally clearer what's being printed at the end

detectSurfaces = True
wrapZ = False # TRUE for a bulk graph (no surfaces, periodic in z), FALSE for a quantum well (surfaces, only periodic in x and y)
export = True

radius = 3.183 # 3.3 for Z=12, 3.189 produces weird bug

print('Constructing supercell.')

na = 32 
nb = 32
nc = 10 #This is the length I am varying.

a = 3.19
c = 5.189
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

tree = spatial.cKDTree(qwAbsolutePositions)

for node in G:
    x = qwAbsolutePositions[node][0]
    y = qwAbsolutePositions[node][1]
    z = qwAbsolutePositions[node][2]
    neighbours = tree.query_ball_point([x,y,z],radius) # where radius is next nearest distance; may be wrong!
    for neighbour in neighbours:
        if neighbour != node:
            G.add_edge(node,neighbour)

neighbourList = []
for node in G:
    neighbourList.append(len(G.edges(node)))
from collections import Counter
print(Counter(neighbourList))

# then shift boundaries to get periodic bonds ...
print('Constructing (shifted x + 0.5) neighbour list.')
for i in range(len(qwFractionalPositions)):
    qwFractionalPositions[i][0] += 0.5
    qwFractionalPositions[i][0] -= math.floor(qwFractionalPositions[i][0])
qw.set_scaled_positions(qwFractionalPositions)
qwAbsolutePositions = qw.get_positions()
tree = spatial.cKDTree(qwAbsolutePositions)

for node in G:
    if len(G.edges(node)) < 12:
        x = qwAbsolutePositions[node][0]
        y = qwAbsolutePositions[node][1]
        z = qwAbsolutePositions[node][2]
        neighbours = tree.query_ball_point([x,y,z],radius) # where radius is next nearest distance; may be wrong!
        for neighbour in neighbours:
            if neighbour != node:
                G.add_edge(node,neighbour)

# and repeat ...
print('Constructing (shifted y + 0.5) neighbour list.')
for i in range(len(qwFractionalPositions)):
    qwFractionalPositions[i][1] += 0.5
    qwFractionalPositions[i][1] -= math.floor(qwFractionalPositions[i][1])
qw.set_scaled_positions(qwFractionalPositions)
qwAbsolutePositions = qw.get_positions()
tree = spatial.cKDTree(qwAbsolutePositions)

for node in G:
    if len(G.edges(node)) < 12:
        x = qwAbsolutePositions[node][0]
        y = qwAbsolutePositions[node][1]
        z = qwAbsolutePositions[node][2]
        neighbours = tree.query_ball_point([x,y,z],radius) # where radius is next nearest distance; may be wrong!
        for neighbour in neighbours:
            if neighbour != node:
                G.add_edge(node,neighbour)

# and repeat ...
print('Constructing (shifted x + 0.5) neighbour list.')
for i in range(len(qwFractionalPositions)):
    qwFractionalPositions[i][0] += 0.5
    qwFractionalPositions[i][0] -= math.floor(qwFractionalPositions[i][0])
qw.set_scaled_positions(qwFractionalPositions)
qwAbsolutePositions = qw.get_positions()
tree = spatial.cKDTree(qwAbsolutePositions)

for node in G:
    if len(G.edges(node)) < 12:
        x = qwAbsolutePositions[node][0]
        y = qwAbsolutePositions[node][1]
        z = qwAbsolutePositions[node][2]
        neighbours = tree.query_ball_point([x,y,z],radius) # where radius is next nearest distance; may be wrong!
        for neighbour in neighbours:
            if neighbour != node:
                G.add_edge(node,neighbour)

neighbourList = []
for node in G:
    neighbourList.append(len(G.edges(node)))
print(Counter(neighbourList))

surfaceNeighbourList = []
bulkNeighbourList = []
for node in G:
	if G.node[node]['surface'] == True:
		surfaceNeighbourList.append(len(G.edges(node)))
	else:
		bulkNeighbourList.append(len(G.edges(node)))
print(Counter(surfaceNeighbourList))
print(Counter(bulkNeighbourList))

# and repeat ...

if wrapZ:
    print('Constructing (shifted z + 0.5) neighbour list.')
    for i in range(len(qwFractionalPositions)):
        qwFractionalPositions[i][2] += 0.5
        qwFractionalPositions[i][2] -= math.floor(qwFractionalPositions[i][2])
    qw.set_scaled_positions(qwFractionalPositions)
    qwAbsolutePositions = qw.get_positions()
    tree = spatial.cKDTree(qwAbsolutePositions)

    for node in G:
        if len(G.edges(node)) < 12:
            x = qwAbsolutePositions[node][0]
            y = qwAbsolutePositions[node][1]
            z = qwAbsolutePositions[node][2]
            neighbours = tree.query_ball_point([x,y,z],radius) # where radius is next nearest distance; may be wrong!
            for neighbour in neighbours:
                if neighbour != node:
                    G.add_edge(node,neighbour)

    neighbourList = []
    for node in G:
        neighbourList.append(len(G.edges(node)))
    print(Counter(neighbourList))

    surfaceNeighbourList = []
    bulkNeighbourList = []
    for node in G:
    	if G.node[node]['surface'] == True:
    		surfaceNeighbourList.append(len(G.edges(node)))
    	else:
    		bulkNeighbourList.append(len(G.edges(node)))
    print(Counter(surfaceNeighbourList))
    print(Counter(bulkNeighbourList))
    
    print('Constructing (shifted x + 0.5) neighbour list.')
    for i in range(len(qwFractionalPositions)):
        qwFractionalPositions[i][0] += 0.5
        qwFractionalPositions[i][0] -= math.floor(qwFractionalPositions[i][0])
    qw.set_scaled_positions(qwFractionalPositions)
    qwAbsolutePositions = qw.get_positions()
    tree = spatial.cKDTree(qwAbsolutePositions)

    for node in G:
        if len(G.edges(node)) < 12:
            x = qwAbsolutePositions[node][0]
            y = qwAbsolutePositions[node][1]
            z = qwAbsolutePositions[node][2]
            neighbours = tree.query_ball_point([x,y,z],radius) # where radius is next nearest distance; may be wrong!
            for neighbour in neighbours:
                if neighbour != node:
                    G.add_edge(node,neighbour)

    neighbourList = []
    for node in G:
        neighbourList.append(len(G.edges(node)))
    print(Counter(neighbourList))

    surfaceNeighbourList = []
    bulkNeighbourList = []
    for node in G:
    	if G.node[node]['surface'] == True:
    		surfaceNeighbourList.append(len(G.edges(node)))
    	else:
    		bulkNeighbourList.append(len(G.edges(node)))
    print(Counter(surfaceNeighbourList))
    print(Counter(bulkNeighbourList))
    
    print('Constructing (shifted y + 0.5) neighbour list.')
    for i in range(len(qwFractionalPositions)):
        qwFractionalPositions[i][1] += 0.5
        qwFractionalPositions[i][1] -= math.floor(qwFractionalPositions[i][1])
    qw.set_scaled_positions(qwFractionalPositions)
    qwAbsolutePositions = qw.get_positions()
    tree = spatial.cKDTree(qwAbsolutePositions)

    for node in G:
        if len(G.edges(node)) < 12:
            x = qwAbsolutePositions[node][0]
            y = qwAbsolutePositions[node][1]
            z = qwAbsolutePositions[node][2]
            neighbours = tree.query_ball_point([x,y,z],radius) # where radius is next nearest distance; may be wrong!
            for neighbour in neighbours:
                if neighbour != node:
                    G.add_edge(node,neighbour)

    neighbourList = []
    for node in G:
        neighbourList.append(len(G.edges(node)))
    print(Counter(neighbourList))

    surfaceNeighbourList = []
    bulkNeighbourList = []
    for node in G:
    	if G.node[node]['surface'] == True:
    		surfaceNeighbourList.append(len(G.edges(node)))
    	else:
    		bulkNeighbourList.append(len(G.edges(node)))
    print('Surface count:')
    print(Counter(surfaceNeighbourList))
    print('Bulk count:')
    print(Counter(bulkNeighbourList))

for node in G:
	if len(G.edges(node)) == 8:
		print G.node[node]

if export:
    if wrapZ:
        nx.write_gpickle(G,str(na)+'x'+str(nb)+'x'+str(nc)+'-bulk-graph-Ganchenkova-050514.gpickle')
    else:
        nx.write_gpickle(G,str(na)+'x'+str(nb)+'x'+str(nc)+'-graph-Ganchenkova-050514.gpickle')
        
print 'Graph constructed and exported successfully.'