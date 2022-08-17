import numpy as np


def Randomise_Position(AP_all, b, c):
    """
    Checks for a viable location to move a chain. This position must not be near any atom of another chain
    (atleast 1.5A away from each other chain member). When found - moves the chain to this position.

    :param AP_all: Atom positions of all the atoms in each chain. All chains should be centered around the origin.
    :param b: Lattice parameter in the y-axis. Set-up by Build_randbox function.
    :param c: Lattice parameter in the z-axis. Set-up by Build_randbox function.

    :return: Atom positions of all atoms in each chain. All chains should now be seperated and placed randomly in box.
    """


    n_AP_all = []

    count = 0

    for i in AP_all:  # sequentially move through list of N polymer fragment positions
        current_chain = []

        flag = True
        while flag:  # turned false only when length of current_chain = length of i
            rand_x = 0.0
            rand_y = np.random.uniform(0.0, float(b))
            rand_z = np.random.uniform(0.0, float(c))
            redo = False

            for j in i:  # iterate through all atom positions of a single chain
                rand_vector = np.array([rand_x, rand_y, rand_z])
                perturbed_j = j + rand_vector

                # check if n_AP_all is empty - if so, fill current chain and add to n_AP_all.
                if len(n_AP_all) > 0:

                    # check, that current_chain for EACH CHAIN, does not overlap its atoms with any of their atoms.
                    for chain in n_AP_all:
                        for atom in chain:
                            if (perturbed_j[0] - atom[0])**2 + (perturbed_j[1] - atom[1])**2 + (perturbed_j[2] - atom[2])**2 < 2**2 \
                                    or (perturbed_j[0] - atom[0])**2 + (perturbed_j[1] - (atom[1]+b))**2 + (perturbed_j[2] - atom[2])**2 < 2**2 \
                                    or (perturbed_j[0] - atom[0]) ** 2 + (perturbed_j[1] - (atom[1])) ** 2 + (perturbed_j[2] - (atom[2]+c)) ** 2 < 2 ** 2 \
                                    or (perturbed_j[0] - atom[0]) ** 2 + (perturbed_j[1] - (atom[1] + b)) ** 2 + (perturbed_j[2] - (atom[2]+c)) ** 2 < 2 ** 2:

                                current_chain = []

                                redo = True
                                # if overlaps, then break out of for loops before satisfying the below if statement.
                                break

                        else:
                            continue

                        break

                    if redo:
                        break

                    else:
                        current_chain.append(perturbed_j)

                else:
                    # fill current chain with all j, non-perturbed. should repeat until end of i=0 (first chain).
                    current_chain.append(j)

            if len(current_chain) == len(i):
                count += 1
                n_AP_all.append(current_chain)
                current_chain = []  # resets current_chain
                flag = False
                print("Chain " + str(count) + " placed in box.")

    return n_AP_all





