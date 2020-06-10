from xyz2graph import MolGraph, to_networkx_graph, to_plotly_figure
from plotly.offline import offline
import time
import sys

start = time.time()

# Create the MolGraph object
mg = MolGraph()

# Read the data from the .xyz file
filename = sys.argv[1]
mg.read_xyz(filename)

if '/' in filename:
    print('input with forward slash')
    size = filename.split("/")[-1]
    size = size.split("-")[0]
elif "\\" in filename:
    print('input with back slash')
    size = filename.split("\\")[-1]
    size = size.split("-")[0]
else:
    size = filename.split("-")[0]

# Create the Plotly figure object
fig = to_plotly_figure(mg)

# Plot the figure
offline.plot(fig, filename=f"{size}-quantum-well")

end = time.time()

print(f"script took this much time to execute:\t{end - start}")

# Convert the molecular graph to the NetworkX graph
# G = to_networkx_graph(mg)