import numpy as np


def Randomise_PE(AP_all_chains, AMP):
    # Takes all the normalised Polyethylene from Build_PE and applies a random perturbation to all atoms.
    # Ensure the chain does not self-overlap... model each atom as hard spheres, and do not let their spheres overlap..
    # If they overlap, reperturb, or move slightly until its ok.

    New_AP_all_chains = [] # Perturbed version of N polymer chains

    for i in AP_all_chains: # sequentially move through list of N polymer fragment positions
        New_AP = []  # separate A_P perturbations for i polymer.

        for j in i: #Sequentially move through list of atoms for i polymer.
            # Checks for C or H (1, 4, 7 .. = C) (2,3, 5,6 .. = H)
            # Perturb j by AMP
            flag = True  #True, until perturbed_j successfully appended to new_AP
            while flag:
                rand = np.random.uniform(-AMP, AMP, size=(3))
                perturbed_j = j+rand

                if len(New_AP) > 0:

                    for k in range(len(New_AP)):
                        #check if new atom (modelled as a point) is not within sphere of radius 0.75A of any atom in New_AP
                        if (perturbed_j[0] - New_AP[k][0])**2 + (perturbed_j[1] - New_AP[k][1])**2 + (perturbed_j[2] - New_AP[k][2])**2 < 0.75**2:
                            break

                    else:
                        flag = Ftalse
                        New_AP.append(perturbed_j)
                else:
                    flag = False
                    New_AP.append(perturbed_j)

        New_AP_all_chains.append(New_AP)

    return New_AP_all_chains

def Randomise_Position(AP_all_chains, AMP):
    #Functionally similar to Randomise_PE, but compares between each polyethylene fragment
    #Rather than self-comparison of each atom position with the same fragment atoms
    New_AP_all_chains = []
    for i in AP_all_chains:  # sequentially move through list of N polymer fragment positions
        for j in i:  # Sequentially move through list of atoms for i polymer.
            # Checks for C or H (1, 4, 7 .. = C) (2,3, 5,6 .. = H)
            # Perturb j by AMP
            flag = True  # True, until perturbed_j successfully appended to new_AP
            redo = True
            while flag:
                if redo:
                    rand = np.random.uniform(-AMP, AMP, size=(3))
                perturbed_j = j + rand

                redo = False
                current_chain = []

                if len(New_AP_all_chains) > 0:  # We do not perturb the position of the first chain.
                    # check through New_Pos (a list of arrays of atom position for each chain)
                    # for each member of New_Pos, find the highest and lowest values for x,y,z.
                    # then check that any member of current_chain is not within 1.5A of these values.

                    for k in New_AP_all_chains:
                        #Currently a list of atom positions in New_AP_all_chains
                        for l in k:
                            #Currently each atom position...
                            #Do a comparison with current atom position
                            if (perturbed_j[0] - l[0]) ** 2 + (perturbed_j[1] - l[1]) ** 2 + (
                                    perturbed_j[2] - l[2]) ** 2 < 1.5 ** 2:
                                redo = True
                                break
                        else:
                            current_chain.append(perturbed_j)
                            continue

                        break
                    flag = False


                else:
                    current_chain.append(j) # If nothing in allchain
                New_AP_all_chains.append(current_chain)

            #   print ("Successfully translated chain position without conflict.")

    return New_AP_all_chains


