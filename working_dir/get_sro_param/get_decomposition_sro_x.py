import networkx as nx
import math
import random
import os
import numpy as np
import sys

import time

# Current structure of program
# 1. Loads graph made from makeCrystalGraph.py & Randomises initial structure
# 3. Runs IsingModel()
# 4. User can then choose to run additional functions, e.g.
#    graphSRO()
#    isingModel()

start = time.time()


class SimulateAnnealing:

    def __init__(self, filename):

        # For graphing short-range order parameters
        self.SROx = []
        self.SROy1 = []
        self.SROy2 = []
        self.SROy3 = []
        self.SROy4 = []

        # Scientific constants
        kB = 8.6173324 * 10e-5  # eV

        # Simulation parameters
        self.T = 773  # K
        self.kT = kB * self.T
        # self.steps = 5000000
        self.steps = 200000
        self.In_composition = 0.5
        self.edgesDict = {}  # hasmhmap to reduce timecomplexity
        self.statsFrequency = 10000  # generate statistics every statsFrequency steps
        self.saveFrequency = int(self.steps) / 1000

        # define binding energies for this specific composition: since it only needs to be determined once.
        self.Eb_1ac = self.get_Eb_1ac()
        self.Eb_1aa = self.get_Eb_1aa()
        self.Eb_2ac = self.get_Eb_2ac()
        self.Eb_2cc = self.get_Eb_2cc()

        # handles linux vs. windows issues with feeding input through command line
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

        # define variables for naming and graphing
        self.size = size
        self.J = []  # J will be our list that holds our nodes where surface == False
        self.H = None  # H will be the subgraph holding the nodes defined in J

        # Loads our initial graph that was made from makeCrystalGraph.py
        print("Loading initial graph.")
        try:
            self.G = nx.read_gpickle(filename)
        except:
            raise Exception("Couldn't load graph! Generate a graph first using makeCrystalGraph.py")
        print("Graph loaded.")

    def initialise_microstructure(self):

        print('size is', self.size.split('x')[-1])
        if self.size.split('x')[-1] == '1':  # if we are simulating a nxnx1 quantum well
            for node in self.G:
                if self.G.nodes[node]['surface']:
                    self.J.append(node)
            self.H = self.G.subgraph(self.J)  # H is the subgraph which contains only surface atoms.
        else:
            for node in self.G:
                if not self.G.nodes[node]['surface']:
                    self.J.append(node)
            self.H = self.G.subgraph(self.J)  # H is the subgraph which contains only BULK atoms.

        # This speeds up the code as we are NOT attempting to swap the fixed surface atoms
        print("Randomising initial configuration to match an overall composition of InₓGa₁₋ₓN "
              f"(x={self.In_composition})")
        tot_no_Ga_atms = 0
        curr_no_In_atms = 0

        for node in self.H:  # count up number of Ga sites in our non-surface subgraph
            if self.H.nodes[node]['species'] == 'Ga':
                tot_no_Ga_atms += 1
            if self.H.nodes[node]['species'] == 'In':
                curr_no_In_atms += 1
        if curr_no_In_atms != 0:
            raise Exception("Code assumes that there are no In atoms in the initial graph.")
        target_no_In_atoms = int(math.floor(self.In_composition * tot_no_Ga_atms))

        for i in range(target_no_In_atoms):
            randomNode = self.J[random.randint(0, len(self.J) - 1)]  # find a random node
            while self.H.nodes[randomNode]['species'] != 'Ga':
                randomNode = self.J[random.randint(0, len(self.J) - 1)]
            self.H.nodes[randomNode]['species'] = 'In'  # change node to In
            curr_no_In_atms += 1

        print(
            f"\nInitial composition randomised. Out of a total of {tot_no_Ga_atms} Ga sites,"
            f"{curr_no_In_atms} are now occupied by In atoms for"
            f"an overall composition of {float(curr_no_In_atms) / tot_no_Ga_atms}.")

    def attempt_atom_swap(self):
        """
        See thesis for definitions:

        What I do here is:

        Pick an In and a Ga atom at random (not at the surface)

        Get the TOTAL number of In atoms around the In atom @ all 4 weights with get_no_In_atms_in_1ac/1aa/2ac/2cc()

        Scale that by Z(k)/2 to account for a geometrical restriction - this massively simplifies our process

        Multiply that by the appropriate binding energy of those configurations to find the
        equilibrium energy of the local environment BEFORE the swap.

        Do the same for the new position that the Indium would be AFTER the swap.
        So treat the neighbours around the Gallium atoms position AS IF it were an Indium atom.

        Hence we +1 to account for the fact that 1 In atom is at the origin of the crystal.

        Plug the value into the ProbabilityOfAcceptance and the rest is history.

        Note:
         o Ein = pre-swap
         o Efi = post-swap
        """

        rnd_In_atm = None
        rnd_Ga_atm = None

        # pick an In atom at random NOT in the surface
        while not rnd_In_atm or not rnd_Ga_atm:
            rnd_atm = self.J[random.randint(0, len(self.J) - 1)]
            if self.H.nodes[rnd_atm]['species'] == 'Ga':
                rnd_Ga_atm = rnd_atm
            else:
                rnd_In_atm = rnd_atm

        # This could return a count of ZERO atoms when get_no_In_atms_in_XXX(rnd_In_atm) returns no In neighbours
        # So we add 1 to account for the Indium atom at the origin
        no_In_atms_around_rnd_In_site = (1
                                         + self.get_no_In_atms_in_1ac(rnd_In_atm)
                                         + self.get_no_In_atms_in_1aa(rnd_In_atm)
                                         + self.get_no_In_atms_in_2ac(rnd_In_atm)
                                         + self.get_no_In_atms_in_2cc(rnd_In_atm))

        # Must +1 for Ga as well because there should be an In atom in the position of the Ga post swap
        no_In_atms_around_rnd_Ga_site = (1
                                         + self.get_no_In_atms_in_1ac(rnd_Ga_atm)
                                         + self.get_no_In_atms_in_1aa(rnd_Ga_atm)
                                         + self.get_no_In_atms_in_2ac(rnd_Ga_atm)
                                         + self.get_no_In_atms_in_2cc(rnd_Ga_atm))

        # Multiplying by scaling factors to account for an additional geometrical restriction - Z(k)/2
        scaled_no_InIn1ac__In = no_In_atms_around_rnd_In_site * 3
        scaled_no_InIn1aa__In = no_In_atms_around_rnd_In_site * 3
        scaled_no_InIn2ac__In = no_In_atms_around_rnd_In_site * 3
        scaled_no_InIn2cc__In = no_In_atms_around_rnd_In_site

        scaled_no_InIn1ac__Ga = no_In_atms_around_rnd_Ga_site * 3
        scaled_no_InIn1aa__Ga = no_In_atms_around_rnd_Ga_site * 3
        scaled_no_InIn2ac__Ga = no_In_atms_around_rnd_Ga_site * 3
        scaled_no_InIn2cc__Ga = no_In_atms_around_rnd_Ga_site

        # Multiple scaled number of atoms in all configurations by their respective binding energy
        Ein = (scaled_no_InIn1ac__In * self.Eb_1ac
               + scaled_no_InIn1aa__In * self.Eb_1aa
               + scaled_no_InIn2ac__In * self.Eb_2ac
               + scaled_no_InIn2cc__In * self.Eb_2cc)

        Efi = (scaled_no_InIn1ac__Ga * self.Eb_1ac
               + scaled_no_InIn1aa__Ga * self.Eb_1aa
               + scaled_no_InIn2ac__Ga * self.Eb_2ac
               + scaled_no_InIn2cc__Ga * self.Eb_2cc)

        # Attempt the Monte-Carlo Step
        probabilityOfAcceptance = math.exp((Ein - Efi) / self.kT)
        acceptSwap = False
        if Ein > Efi or random.random() < probabilityOfAcceptance:  # if new config is energetically favourable
            acceptSwap = True
            self.G.nodes[rnd_In_atm]['species'], self.G.nodes[rnd_Ga_atm]['species'] = 'Ga', 'In'

        return acceptSwap

    def get_Eb_1ac(self):
        """
        Determines the binding energy at weight 1ac as a function of the indium composition
        :return: Binding Energy (eV)
        """

        E1 = (0.4945 * (self.In_composition ** 4)
              - 0.727 * (self.In_composition ** 3)
              + 0.3929 * (self.In_composition ** 2)
              - 0.0925 * self.In_composition - 0.072)

        return E1

    def get_Eb_1aa(self):
        """
        Determines the binding energy at weight 1aa as a function of the indium composition
        :return: Binding Energy (eV)
        """

        if self.In_composition == 0.1:
            E2 = -0.065
        else:
            E2 = 0.0738 * self.In_composition - 0.0771

        return E2

    def get_Eb_2ac(self):
        """
        Determines the binding energy at weight 2ac as a function of the indium composition
        :return: Binding Energy (eV)
        """

        if self.In_composition == 0.1:
            E3 = 0.007
        elif self.In_composition == 0.25:
            E3 = 0.008
        else:
            E3 = 0.0133 * self.In_composition + 0.0052

        return E3

    def get_Eb_2cc(self):
        """
        Determines the binding energy at weight 2cc as a function of the indium composition
        :return: Binding Energy (eV)
        """

        E4 = 0.1238 * (self.In_composition ** 2) - 0.1122 * self.In_composition + 0.04

        return E4

    def get_no_In_atms_in_1ac(self, node):
        """
        input a node, return the number of neighbouring In atoms in the 1ac configuration (weight=1)
        """

        neighboursOfRandomNode = self.getNeighbours(node, 1)
        numberOfIndiumNeighbours = 0
        for neighbour in neighboursOfRandomNode:
            if self.G.nodes[neighbour]['species'] == 'In':
                numberOfIndiumNeighbours += 1
        return numberOfIndiumNeighbours

    def get_no_In_atms_in_1aa(self, node):
        """
        input a node, return the number of neighbouring In atoms in the 1aa configuration (weight=2)
        """

        neighboursOfRandomNode = self.getNeighbours(node, 2)
        numberOfIndiumNeighbours = 0
        for neighbour in neighboursOfRandomNode:
            if self.G.nodes[neighbour]['species'] == 'In':
                numberOfIndiumNeighbours += 1
        return numberOfIndiumNeighbours

    def get_no_In_atms_in_2ac(self, node):
        """
        input a node, return the number of neighbouring In atoms in the 2ac configuration (weight=3)
        """

        neighboursOfRandomNode = self.getNeighbours(node, 3)
        numberOfIndiumNeighbours = 0
        for neighbour in neighboursOfRandomNode:
            if self.G.nodes[neighbour]['species'] == 'In':
                numberOfIndiumNeighbours += 1
        return numberOfIndiumNeighbours

    def get_no_In_atms_in_2cc(self, node):
        """
        input a node, return the number of neighbouring In atoms in the 2ccc configuration (weight=4)
        """
        neighboursOfRandomNode = self.getNeighbours(node, 4)
        numberOfIndiumNeighbours = 0
        for neighbour in neighboursOfRandomNode:
            if self.G.nodes[neighbour]['species'] == 'In':
                numberOfIndiumNeighbours += 1
        return numberOfIndiumNeighbours

    def run_metropolis_monte_carlo(self, num_steps):
        """
        Runs attempt_atom_swap() self.num_steps times.
        """

        ns = float(num_steps)

        print("Running Ising Model to calculate SRO parameters... Please wait.")
        successfulSwaps = 0
        SwapsSoFar = 0
        for i in range(num_steps):
            SwapsSoFar += 1
            if self.attempt_atom_swap():
                successfulSwaps += 1
            if (i % self.statsFrequency == 0) or i == 0 or i == 10 or i == 100 or i == 500:  # calculate SRO prameters
                print(round(float(i) / ns * 100, 3), '% Done')
                self.SROx.append(i)
                self.SROy1.append(self.get_sro_x(1))
                self.SROy2.append(self.get_sro_x(2))
                self.SROy3.append(self.get_sro_x(3))
                self.SROy4.append(self.get_sro_x(4))
            if i % self.saveFrequency == 0 and i > 2:  # uncomment to create multiple microstructures
                # self.writeGraphToXYZ(self.size + '-xIn=' + str(self.In_composition) + '_T=' + str(self.T)
                #                      + 'K_SRO_X=' + str(self.SROy1[-1]) + '.xyz')
                pass
        print(f"\nOut of {num_steps} attempted atom swaps, {successfulSwaps} were successful "
              f"(ratio of {successfulSwaps / num_steps}).")

        return successfulSwaps

    def get_sro_x(self, weight):
        """
        Returns the Warren-Cowley short range order parameter for the first
        nearest neighbour (x=weight). See Chan, Liu and Zunger paper for more
        information. DOI: 10.1103/PhysRevB.82.045112
        """

        indiumCount = []  # - O(n*len(neighbours)) - n is proportional to the volume of q_well but the latter is small
        for node in self.H:  # creating a Ga-only subgraph will reduce this to O(n)--> O(n*(1-Indium_composition))
            if self.H.nodes[node]['species'] == 'Ga':
                neighbours = self.getNeighbours(node, weight)
                numberOfIndiumNeighbours = 0
                for neighbour in neighbours:
                    if self.H.nodes[neighbour]['species'] == 'In':
                        numberOfIndiumNeighbours += 1
                # e.g. for indiumComposition=0.5, numberofIndiumNeighbours would be 3 on average for random compositions
                # but might fluctuate from site to site, and therefore indiumCount would be appended by 0.5 on average
                if weight != 4:
                    indiumCount.append(numberOfIndiumNeighbours / 6.0)
                else:
                    indiumCount.append(numberOfIndiumNeighbours / 2.0)

        indiumProbability = sum(indiumCount) / len(indiumCount)  # this calculates that average!

        return 1 - (indiumProbability / self.In_composition)

    def getNeighbours(self, node, weight):
        """
        input your node and the weight of the neighbours you wanna find
        Retrieve a list of the nodes of the kth neigbour.
        """

        # check if we already have this key (reduce time complexity to O(1))
        key = (node, weight)
        if key in self.edgesDict:
            return self.edgesDict[key]
        else:  # O(n) time complexity
            neighbours_with_wt = []
            all_neighbours = list(self.H[node].keys())
            for a_nbr in all_neighbours:
                if self.H[node][a_nbr]['weight'] == weight:
                    neighbours_with_wt.append(a_nbr)

            self.edgesDict[key] = neighbours_with_wt
            return neighbours_with_wt

    def writeGraphToXYZ(self, filename):
        f = open(filename, 'w')
        f.write(str(len(self.G)) + '\n')
        f.write('Atoms. File created from networkx graph by ' + os.path.basename(__file__) + '\n')
        for node_index in self.G:
            atom = self.G.nodes[node_index]
            f.write(str(atom['species']) + ' ' + str(atom['x']) + ' ' + str(atom['y']) + ' ' + str(atom['z']) + '\n')
        f.close()
        print("Graph exported as .xyz file.")

    def save_output(self):
        SROx_out = np.transpose([self.SROx, self.SROy1, self.SROy2, self.SROy3, self.SROy4])
        np.savetxt('CSV_' + self.size + '-xIn=' + str(self.In_composition) + '_T=' + str(self.T) + 'K_SRO_X' + '.csv',
                   SROx_out, delimiter=',')
        self.writeGraphToXYZ(self.size + '-xIn=' + str(self.In_composition) + '_T=' + str(self.T) + 'K_SRO_X=' +
                             str(self.SROy1[-1]) + '.xyz')
        nx.write_gpickle(self.G, self.size + '-xIn=' + str(self.In_composition) + '_T=' + str(self.T) + 'K_SRO_X=' +
                         str(self.SROy1[-1]) + '.gpickle')


########################################################################################################################
# Driver Code
if __name__ == "__main__":

    try:
        arg = sys.argv[1]
    except IndexError:
        arg_zero = sys.argv[0]
        raise SystemExit("Usage: " + str(arg_zero) + " <input_filename.gpickle>")

    if len(sys.argv) > 2:
        raise SystemExit("Please enter the .gpickle file of the crystal structure to decompose")

    sa = SimulateAnnealing(sys.argv[1])
    sa.initialise_microstructure()
    sa.run_metropolis_monte_carlo(sa.steps)
    sa.save_output()

    end = time.time()
    print(end - start)
