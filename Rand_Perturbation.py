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
                        flag = False
                        New_AP.append(perturbed_j)
                else:
                    flag = False
                    New_AP.append(perturbed_j)

        New_AP_all_chains.append(New_AP)

    return New_AP_all_chains

def Randomise_Position(AP_all_chains, a, b, c):
    #Functionally similar to Randomise_PE, but compares between each polyethylene fragment
    #Rather than self-comparison of each atom position with the same fragment atoms
    New_AP_all_chains = []
    count=0

    for i in AP_all_chains:  # sequentially move through list of N polymer fragment positions
        current_chain = []

        flag = True
        while flag:  # turned false only when length of current_chain = length of i
            rand_x = 0.0
            rand_y = np.random.uniform(0.0, float(b))
            rand_z = np.random.uniform(0.0, float(c))
            redo = False

            for j in i: #iterate through all atom positions of a single chain
                rand_vector = np.array([rand_x, rand_y, rand_z])
                perturbed_j = j + rand_vector

                # check if NewAPallchains is empty - if so, fill current chain and add to newallchains
                if len(New_AP_all_chains) > 0:

                    # check, that current_chain for EACH CHAIN, does not overlap its atoms with any of their atoms.
                    for chain in New_AP_all_chains:
                        for atom in chain:
                            if (perturbed_j[0] - atom[0])**2 + (perturbed_j[1] - atom[1])**2 + (perturbed_j[2] - atom[2])**2 < 2**2:
                                current_chain = []

                                redo = True
                                # if overlaps, then break out of for loops before satisfying the below if statement.
                                break

                        else:
                            continue

                        break

                    if redo == True:
                        break

                    else:
                        current_chain.append(perturbed_j)

                else:
                    # fill current chain with all j, non-perturbed. should repeat until end of i=0 (first chain).
                    current_chain.append(j)




            if len(current_chain) == len(i):
                count +=1
                New_AP_all_chains.append(current_chain)
                current_chain = [] #resets current_chain
                flag = False
                print ("Chain " + str(count) + " placed in box.")


    return New_AP_all_chains





