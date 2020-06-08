import networkx as nx
import random
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


def generate_random_3Dgraph(n_nodes, radius, seed=None):

    """
    The inputs are the number of nodes in the network, the radius value and, obviously, the position of the nodes. The parameters radius is used to specify if two nodes are connected or not. Two nodes are connected by an edge if their distance is at most equal to radius. This will give the number of edges, that is the connections between the nodes.
    """

    if seed is not None:
        random.seed(seed)

    # Generate a dict of positions
    pos = {i: (random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)) for i in range(n_nodes)}
    print(pos)
    # Create random 3D network
    G = nx.random_geometric_graph(n_nodes, radius, pos=pos)

    return G

# visualise with size and color

def network_plot_3D(G, angle, save=False):

    """
    A couple of things to note here.

    To draw the nodes we use 3D scatter plot
    To draw the edges we use 3D line plot
    Note how to define the colour of the node: we get the value of the maximum number of edges in a single node, and use that value to define the colour scale to go from zero to such a maximum value.
    """

    # Get node positions
    pos = nx.get_node_attributes(G, 'pos')

    # Get number of nodes
    n = G.number_of_nodes()
    # Get the maximum number of edges adjacent to a single node
    edge_max = max([G.degree(i) for i in range(n)])
    # Define color range proportional to number of edges adjacent to a single node
    colors = [plt.cm.plasma(G.degree(i)/edge_max) for i in range(n)] 
    # 3D network plot
    with plt.style.context(('ggplot')):

        fig = plt.figure(figsize=(10,7))
        ax = Axes3D(fig)

        # Loop on the pos dictionary to extract the x,y,z coordinates of each node
        for key, value in pos.items():
            xi = value[0]
            yi = value[1]
            zi = value[2]

            # Scatter plot
            ax.scatter(xi, yi, zi, c=colors[key], s=20+20*G.degree(key), edgecolors='k', alpha=0.7)

        # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
        # Those two points are the extrema of the line to be plotted
        for i,j in enumerate(G.edges()):
            x = np.array((pos[j[0]][0], pos[j[1]][0]))
            y = np.array((pos[j[0]][1], pos[j[1]][1]))
            z = np.array((pos[j[0]][2], pos[j[1]][2]))

        # Plot the connecting lines
            ax.plot(x, y, z, c='black', alpha=0.5)

    # Set the initial view
    ax.view_init(30, angle)
    # Hide the axes
    ax.set_axis_off()
    if save is not False:
         plt.savefig("example_graph"+".png")
         plt.close('all')
    else:
          plt.show()

# use a simple linear scaling to scale the size of each node.
#s=20+20*G.degree(key)


# visualise with 200 nodes

n = 200
Gee = generate_random_3Dgraph(n_nodes=n, radius=0.25, seed=1)
network_plot_3D(Gee, 0, save=False)

# for k in range(20, 201, 1):
#     G = generate_random_3Dgraph(n_nodes=k, radius=0.25, seed=1)
#     angle = (k - 20) * 360 / (200 - 20)
#
#     network_plot_3D(G, angle, save=True)
#     print(angle)