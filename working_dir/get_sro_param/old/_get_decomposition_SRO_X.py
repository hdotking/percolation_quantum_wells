import networkx as nx
import math
import random
import os
import numpy
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
        self.steps = 100000
        self.In_composition = 0.5
        self.edgesDict = {}
        self.sro_1_dict = {}
        self.sro_2_dict = {}
        self.sro_3_dict = {}
        self.sro_4_dict = {}
        self.statsFrequency = 10000  # generate statistics every statsFrequency steps
        self.saveFrequency = int(self.steps) / 1000

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

        self.size = size
        self.J = []  # J will be our list that holds our nodes where surface == False

        # Loads our initial graph that was made from makeCrystalGraph.py
        print("Loading initial graph.")
        try:
            self.G = nx.read_gpickle(filename)
        except:
            raise Exception("Couldn't load graph! Generate a graph first using makeCrystalGraph.py")
        print("Graph loaded.")

        self.H = None

    def initialise_microstructure(self):

        print('size is', self.size.split('x')[-1])
        if self.size.split('x')[-1] == '1':
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
        rnd_In_atm = None
        rnd_Ga_atm = None
        # pick an In atom at random NOT in the surface

        while not rnd_In_atm or not rnd_Ga_atm:
            rnd_atm = self.J[random.randint(0, len(self.J) - 1)]
            if self.H.nodes[rnd_atm]['species'] == 'Ga':
                rnd_Ga_atm = rnd_atm
            else:
                rnd_In_atm = rnd_atm

        # +1 includes counts the Indium atom that we are looking for neighbours for in the count
        # Otherwise this could count ZERO In-In pairs when get_N_A_XXX() returns zero In neighbours

        no_In_neighbours_of_In_site = (1
                                       + self.get_no_atms_in_1ac(rnd_In_atm)
                                       + self.get_no_atms_in_1aa(rnd_In_atm)
                                       + self.get_no_atms_in_2ac(rnd_In_atm)
                                       + self.get_no_atms_in_2cc(rnd_In_atm))

        # Must +1 for Ga as well because there should be an In atom in the position of the Ga post swap
        no_In_neighbours_of_Ga_site = (1
                                       + self.get_no_atms_in_1ac(rnd_Ga_atm)
                                       + self.get_no_atms_in_1aa(rnd_Ga_atm)
                                       + self.get_no_atms_in_2ac(rnd_Ga_atm)
                                       + self.get_no_atms_in_2cc(rnd_Ga_atm))

        # multiplying by 3 to satisfy the symmetry of the wurtzite structure
        nconfig_InIn1ac_In = no_In_neighbours_of_In_site * 3
        nconfig_InIn2cc_In = no_In_neighbours_of_In_site
        nconfig_InIn1ac_Ga = no_In_neighbours_of_Ga_site * 3
        nconfig_InIn2cc_Ga = no_In_neighbours_of_Ga_site

        Eb_1ac_comp = self.Eb_1ac()
        Eb_1aa_comp = self.Eb_1aa()
        Eb_2ac_comp = self.Eb_2ac()
        Eb_2cc_comp = self.Eb_2cc()

        Ein = (nconfig_InIn1ac_In * Eb_1ac_comp
               + nconfig_InIn1ac_In * Eb_1aa_comp
               + nconfig_InIn1ac_In * Eb_2ac_comp
               + nconfig_InIn2cc_In * Eb_2cc_comp)

        Efi = (nconfig_InIn1ac_Ga * Eb_1ac_comp
               + nconfig_InIn1ac_Ga * Eb_1aa_comp
               + nconfig_InIn1ac_Ga * Eb_2ac_comp
               + nconfig_InIn2cc_Ga * Eb_2cc_comp)

        probabilityOfAcceptance = math.exp((Ein - Efi) / self.kT)
        acceptSwap = False
        if Ein > Efi or random.random() < probabilityOfAcceptance:
            acceptSwap = True
            self.G.nodes[rnd_In_atm]['species'], self.G.nodes[rnd_Ga_atm]['species'] = 'Ga', 'In'

        return acceptSwap

    def Eb_1ac(self):
        E1 = 0.4945 * (self.In_composition ** 4) - 0.727 * (self.In_composition ** 3) + 0.3929 * (
                self.In_composition ** 2) - 0.0925 * self.In_composition - 0.072
        return E1

    def Eb_1aa(self):
        if self.In_composition == 0.1:
            E2 = -0.065
        else:
            E2 = 0.0738 * self.In_composition - 0.0771
        return E2

    def Eb_2ac(self):
        if self.In_composition == 0.1:
            E3 = 0.007
        elif self.In_composition == 0.25:
            E3 = 0.008
        else:
            E3 = 0.0133 * self.In_composition + 0.0052
        return E3

    def Eb_2cc(self):
        E4 = 0.1238 * (self.In_composition ** 2) - 0.1122 * self.In_composition + 0.04
        return E4

    def get_num_atms_in_1ac(self, node):
        """
        input a node, return the number of In atoms of the 1ac neighbours(weight=1)
        """
        neighboursOfRandomNode = self.getNeighbours(node, 1)
        numberOfIndiumNeighbours = 0
        for neighbour in neighboursOfRandomNode:
            if self.G.nodes[neighbour]['species'] == 'In':
                numberOfIndiumNeighbours += 1
        return numberOfIndiumNeighbours

    def get_num_atms_in_1aa(self, node):
        """
        input a node, return the number of In atoms of the 1aa neighbours(weight=2)
        """
        neighboursOfRandomNode = self.getNeighbours(node, 2)
        numberOfIndiumNeighbours = 0
        for neighbour in neighboursOfRandomNode:
            if self.G.nodes[neighbour]['species'] == 'In':
                numberOfIndiumNeighbours += 1
        return numberOfIndiumNeighbours

    def get_num_atms_in_2ac(self, node):
        """
        input a node, return the number of In atoms of the 2ac neighbours(weight=3)
        """
        neighboursOfRandomNode = self.getNeighbours(node, 3)
        numberOfIndiumNeighbours = 0
        for neighbour in neighboursOfRandomNode:
            if self.G.nodes[neighbour]['species'] == 'In':
                numberOfIndiumNeighbours += 1
        return numberOfIndiumNeighbours

    def get_num_atms_in_2cc(self, node):
        """
        input a node, return the number of In atoms of the 2cc neighbours(weight=4)
        """
        neighboursOfRandomNode = self.getNeighbours(node, 4)
        numberOfIndiumNeighbours = 0
        for neighbour in neighboursOfRandomNode:
            if self.G.nodes[neighbour]['species'] == 'In':
                numberOfIndiumNeighbours += 1
        return numberOfIndiumNeighbours

    def run_Ising_model_sro1(self, num_steps):
        """
        Runs attemptAtomSwap() numberOfSwaps times.
        """

        ns = float(num_steps)

        print("Running Ising Model to calculate SRO_1... Please wait.")
        successfulSwaps = 0
        SwapsSoFar = 0
        for i in range(num_steps):
            SwapsSoFar += 1
            if self.attempt_atom_swap():
                successfulSwaps += 1
            if (i % self.statsFrequency == 0) or i == 0 or i == 10 or i == 100 or i == 500:
                print(round(float(i), 3) / ns * 100, '% Done')
                self.SROx.append(i)
                self.SROy1.append(self.getSRO_X(1))
                self.SROy2.append(self.getSRO_X(2))
                self.SROy3.append(self.getSRO_X(3))
                self.SROy4.append(self.getSRO_X(4))
            if (i % self.saveFrequency == 0 and i > 2):
                # writeGraphToXYZ(size+'-xIn='+str(indiumComposition)+'_T='+str(T)+' K_SRO_1_'+str(SROy1[-1])+'SwapsSoFar'+str(SwapsSoFar)+'.xyz')
                pass
        print("\nOut of %i attempted atom swaps, %i were successful (ratio of %f)." % (
            num_steps, successfulSwaps, float(successfulSwaps) / num_steps))
        return successfulSwaps

    def getSRO_X(self, weight):
        """Returns the Warren-Cowley short range order parameter for the first
        nearest neighbour (x=weight). See Chan, Liu and Zunger paper for more
        information. DOI: 10.1103/PhysRevB.82.045112
        """

        indiumCount = []  # - O(n*len(neighbours)) - n is proportional to the volume of q_well but the latter is <6
        fail_counter = 0
        for node in self.H:  # creating a Ga-only subgraph will reduce this to O(n)--> O(n*(1-Indium_composition))
            if self.H.nodes[node]['species'] == 'Ga':
                neighbours = self.getNeighbours(node, weight)
                print("SRO_" + str(weight), neighbours, len(neighbours), fail_counter)
                numberOfIndiumNeighbours = 0
                for neighbour in neighbours:
                    if self.H.nodes[neighbour]['species'] == 'In':
                        numberOfIndiumNeighbours += 1
                # e.g. for indiumComposition=0.5, numberofIndiumNeighbours would be 3 on average for a random composition
                # but might fluctuate from site to site, and therefore indiumCount would be appended by 0.5 on average
                indiumCount.append(numberOfIndiumNeighbours / 6.0)
            else:
                fail_counter += 1
        indiumProbability = sum(indiumCount) / len(indiumCount)  # this calculates that average!

        print('This is fail counter:', fail_counter)
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
            # python3 implementation
            neighbours_with_wt = []
            all_nbrs = list(self.H[node].keys())
            for a_nbr in all_nbrs:
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
        SROx_out = numpy.asarray([self.SROx, self.SROy1, self.SROy2, self.SROy3, self.SROy4])
        numpy.savetxt('CSV_' + self.size + '-xIn=' + str(self.In_composition) + '_T=' + str(self.T) + 'K_SRO_X' + '.csv'
                      , SROx_out, delimiter=",")
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
    sa.run_Ising_model_sro1(sa.steps)
    sa.save_output()

    end = time.time()

    print(end - start)
    key_list = list(sa.edgesDict.keys())

    # print(sa.nc)
    # do we have all of the neighbours of all of the nodes here?
    # Yes
    # for a 3x3x3 crystal , we have 2 atoms per 'unit cell'
    # therefore 27*2 atoms in total
    # BUT 3*3*2 of those will be surface atoms - two sheets
    # These surface atoms will not be swapping, but will count towards neighbours)
    # Therefore edges dict will NOT use these as keys, but will use them in their values
    # hence the counter below will show 54 - 18 = 36 atoms. The number of keys in edgesDict

    # THIS MEANS THAT edgesDict contains ALL the atoms and weights that we need for getNeighbours()
