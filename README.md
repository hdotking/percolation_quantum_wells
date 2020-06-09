# How to analyse the percolation threshold of non-random alloys - InGaN

![Alt text](working_dir/visualising/xyz_to_graph/InGaNQuantumWell-3x3x3.png)

## Make your crystals
Run the makeQuantumWell.py script to create your own quantum well.
E.g) to get a quantum well which is 10x10x5 InGaN unit cells in size: run the code as follows:

``` 
python makeQuantumWellGraph.py 10 5
>>10x10x5-quantum-well-graph.gpickle
```
## Run the Ising model to simulate annealing
The output file can then then be used as an argument for getDecompositionCalculateSRO1.py 

Visualising the codes output can be done with any software that can view .xyz files - see the visualisation/xyz2graph directory for examples.

We know that decomposition has occured by observing the SRO parameter in the CSV.

IF 0 < SRO_X < 1 --> The system has formed clusters

## OS Compatibility

the output of makeQuantumWellGraph.py differs when build on windows and linux (ubuntu)

This should not cause problems as the graphs can be made on either OS as required.
