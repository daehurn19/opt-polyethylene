import numpy as np


def Randomise_Position(AP_all_chains, b, c):
    #Functionally similar to Randomise_PE, but compares between each polyethylene fragment
    #Rather than self-comparison of each atom position with the same fragment atoms
    New_AP_all_chains = []

    #All mirrors are: +0, +b, +c, +bc
    mirror = []
    count = 0

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
                            if (perturbed_j[0] - atom[0])**2 + (perturbed_j[1] - atom[1])**2 + (perturbed_j[2] - atom[2])**2 < 1.5**2 \
                                    or (perturbed_j[0] - atom[0])**2 + (perturbed_j[1] - (atom[1]+b))**2 + (perturbed_j[2] - atom[2])**2 < 1.5**2 \
                                    or (perturbed_j[0] - atom[0]) ** 2 + (perturbed_j[1] - (atom[1])) ** 2 + (perturbed_j[2] - (atom[2]+c)) ** 2 < 1.5 ** 2 \
                                    or (perturbed_j[0] - atom[0]) ** 2 + (perturbed_j[1] - (atom[1] + b)) ** 2 + (perturbed_j[2] - (atom[2]+c)) ** 2 < 1.5 ** 2:


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





