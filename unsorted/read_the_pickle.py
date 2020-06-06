import networkx as nx


size = "32x32x5"
G = nx.read_gpickle("/home/hdot/py2/git_repos/percolation_quantum_wells/32x32x5_xIn=0.5_T=773 K_SRO_1_0.773274739583Final.gpickle")


J = []
if size.split('x')[-1] == '1':
    for node in G:
        if G.node[node]['surface']:
            print(type(G.node[node]['weight']))
            J.append(node)
    H = G.subgraph(J)  # H is the subgraph which contains only surface atoms.
else:
    for node in G:
        if not G.node[node]['surface']:
            J.append(node)
    H = G.subgraph(J)  # H is the subgraph which contains only surface atoms.
# This speeds up the code
lenJ = len(J)
#
#
# #for node in H:
#     #print(H.nodes[node]['species'])
#
# nx.draw(G)
# plt.show()
#


print "hello world"



#edgesDict = {}
#
#
#def getNeighbours(i, k, P=None):
#    """
#    input your node and the weight of the neighbours you wanna find, i and
#    k. Retrieve a list of the nodes of the kth neigbour.
#    """
#
#    global edgesDict
#    key = (i, k)
#    if key in edgesDict:
#        return edgesDict[key]
#    else:
#        global H
#        if P:
#            H = P
#        edges = H.edges(i, k)
#        edges_of_weight = []
#        neighbours = []
#        print(edges)
#        print(edgesDict)
#        for edge in edges:
#            print(edge)
#            if edge[2]['weight'] == k:
#                edges_of_weight.append(edge)
#        for anElement in edges_of_weight:
#            neighbours.append(anElement[1])
#        edgesDict[key] = neighbours
#        return neighbours
#
#
#getNeighbours(195, 1)
