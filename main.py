import numpy as np
from Build_PE import Build_PE
from Build_randbox import Acceptability
from Rand_Perturbation import Randomise_PE, Randomise_Position
import ase, ase.io, ase.visualize
from Rotate_chain import Rotator
from ase.io.opls import OPLSStructure, OPLSff


def main(N, L, rot=False):

    #Input: N = Number of chains to be created
    #       L = C Length of the polyethylene chain
    #       rot = True/False for random rotation of each polyethylene chain, Defaults to False

    #Output variables
    all_chain_L = []
    num_C = []
    AP_allchains = [] #Atom_positions of all N chains, sorted by generation history
    AT_allchains = [] #Atom_types of all N chains, sorted by generation history
    total_atoms = 0


    #Build N normal polyethylene chains
    for i in range(N):
        AP_allchains.append(Build_PE(L)[0])
        AT_allchains.append(Build_PE(L)[1])
        all_chain_L.append(AP_allchains[-1][-3][0])
        num_C.append(L)
        total_atoms += L*3
    print(str(N) + " standardised polyethylene generated.")
    print("Total expected atoms = " + str(total_atoms))

    # then build random box to contain chains
    flag = True
    rand_box_dim = [] #format is a,b,c,alp,bet,gam
    while flag:
        acceptable = Acceptability(all_chain_L, num_C)

        if acceptable[0]:
            flag = False
            for i in range(6):
                rand_box_dim.append(acceptable[i+1])
            print("Acceptable box configuration found.")

    new_AP_allchains = []

    rand_box_dim[0]= 2.5
    rand_box_dim[1] = 2.5
    rand_box_dim[2] = 2.5
    rand_box_dim[3] = 90
    rand_box_dim[4] = 90
    rand_box_dim[5] = 90


    if rot:
        for i in AP_allchains:
            new_AP_allchains.append(Rotator(i))


    else:
        new_AP_allchains = AP_allchains


    """
    # then perturb each chain
    AMP = 1
    New_AP_allchains = Randomise_PE(AP_allchains, AMP)
    print(New_AP_allchains)
    print("All polyethylene molecules randomised.")
    
"""
    # then place each chain randomly in box # account for PBC
    # Can be done by making sure there is a minimum 1.5A distance from each chain,
    # and all chains stay within the confines of the box dimensions in b and c direction.

    a, b, c = rand_box_dim[0], rand_box_dim[1], rand_box_dim[2]
    new_AP_allchains = Randomise_Position(new_AP_allchains, a, b, c)
    print("All polyethylene fragments randomly placed within box.")

    #Now, flatten all lists for output.
    # Use ASE to do this easily.
    #We have both rand_box_dim (dimensions of random box), and New_AP_allchains (coordiantes of randomised polethylene fragments)

    print ("Flattening elements to write...")
    elements = ''
    for s in AT_allchains:
        for t in s:
            elements += t
    print ("Success.")

    print("Flattening atom positions to write...")
    flat_AP = []
    for s in new_AP_allchains: #New_ap_allchains = list of np.arrays of 3 elements
        for t in s:
            flat_AP.append(t)
    print ("Success.")


    types = []
    atoms = ase.Atoms(elements, np.array(flat_AP), cell=rand_box_dim, pbc=False)



    atoms.arrays['molid'] = np.ones(len(atoms),dtype='int')

    [types.append("CT") if x%3 == 0 else types.append("HC") for x in range(len(atoms))]
    types = np.array(types)
    atoms.arrays['type'] = types


    ase.visualize.view(atoms)
    ase.io.write(r'C:\Users\daehu\OneDrive\Desktop\polyethylene\test.xyz', atoms, format='extxyz')

    print("Output file has been made. Please run on Visualisation software.")


    structure = OPLSStructure(r'C:\Users\daehu\OneDrive\Desktop\polyethylene\test.xyz')
    ase.visualize.view(structure)



    with open(r'C:\Users\daehu\OneDrive\Desktop\polyethylene\FF_defs.par') as fd:
        opls = OPLSff(fd)
    opls.write_lammps(structure, prefix='lmp')


    return -1



main(1, 2, rot=True) # 3 chains of L = 5


