import numpy as np
from Build_PE import Build_PE
from Build_randbox import Acceptability
from Rand_Perturbation import Randomise_PE, Randomise_Position
import ase, ase.io

def main(N): # Number of chains to create inside box.
    all_chain_L = []
    num_C = []
    AP_allchains = [] #Atom_positions of all N chains, sorted by generation history
    AT_allchains = [] #Atom_types of all N chains, sorted by generation history
    #build N normal polyethylene chains
    for i in range(N):
        ran_num = np.random.randint(5, 20)
        AP_allchains.append(Build_PE(2)[0])
        AT_allchains.append(Build_PE(2)[1])
        all_chain_L.append(AP_allchains[-1][-3][0])
        num_C.append(ran_num)
    print(str(N) + " standardised polyethylene generated.")

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

    # then perturb each chain
    AMP = 2
    New_AP_allchains = Randomise_PE(AP_allchains, AMP)
    print(New_AP_allchains)
    print("All polyethylene molecules randomised.")

    # then place each chain randomly in box # account for PBC
    # Can be done by making sure there is a minimum 1.5A distance from each chain,
    # and all chains stay within the confines of the box dimensions in b and c direction.
    """AMP = 3
    New_AP_allchains = Randomise_Position(New_AP_allchains, AMP)
    print (New_AP_allchains)
    print ("Polyethylene randomly placed within box.")
"""

    #Now, flatten all lists for output.
    # Use ASE to do this easily.
    #We have both rand_box_dim (dimensions of random box), and New_AP_allchains (coordiantes of randomised polethylene fragments)

    elements = ''
    for s in AT_allchains:
        for t in s:
            elements += t
    print (len(elements))

    flat_AP = []
    print (len(New_AP_allchains))
    print (New_AP_allchains)
    for s in New_AP_allchains: #New_ap_allchains = list of np.arrays of 3 elements
        for t in s:
            flat_AP.append(t)
    print (flat_AP)


    atoms = ase.Atoms(elements, np.array(flat_AP), cell=rand_box_dim, pbc=True)
    ase.io.write(r'C:\Users\daehu\Desktop\polyethylene\test.cell', atoms, format='castep-cell')
    return -1



print(main(3))


