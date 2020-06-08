import networkx as nx


def getEdgesWithWeight(node, weight):
    edges = G.edges(node, weight)
    selected_edges = []
    for edge in edges:
        # the line below used to be edge[2]['weight'] == weight. is this doing the same?

        if G[edge[0]][edge[1]]["weight"] == weight: # get the weight for each edge and see if it"s the weight we want
            selected_edges.append(edge)
    return selected_edges

J = nx.Graph()
J.add_nodes_from([1, 2, 3, 4, 5])
J.add_edges_from([(1, 2, {'weight': 1}), (1, 3, {'weight': 2}), (1, 4, {'weight': 3}), (1, 5, {'weight': 4})])

print(nx.get_node_attributes(J, 'weight'))

