import numpy as np


test = [np.array([[0,0,0],[0,0.45,0.8],[0,0.45,-0.8]]),np.array([[0,0,0],[0,0.45,0.8],[0,0.45,-0.8]])]
test2 = [['C','H','H'],['C,H,H']]
def Randomise_PE(AP_all_chains, AMP):
    # Takes a normalised Polyethylene and applies a random perturbation to all atoms in the chain.
    # Ensure the chain does not self-overlap, save each new position in a list, and re-perturb if new member is within
    #  0.75A same spatial coordinates of C or 0.5A same spatial coordinates of H

    New_AP_all_chains = [] # Perturbed version of N chains


    for i in AP_all_chains: # sequentially move through list of N polymer fragment positions

        New_AP = []  # seperate A_P perturbations for i polymer.

        for j in i: #Sequentially move through list of A_P for i polymer.
            c = len(New_AP)  # Checks for C or H (1, 4, 7 .. = C) (2,3, 5,6 .. = H)
            print (str(j) + " current j")
            # Perturb j by AMP
            flag = True #True, until perturbed_j successfully appended to new_AP
            while flag:
                rand = np.random.uniform(-AMP, AMP, size=(3))
                perturbed_j = j+rand
                abs_j = np.sqrt(perturbed_j[0]**2 + perturbed_j[1]**2 + perturbed_j[2]**2)

                if len(New_AP) > 0:

                    for k in range(c):
                        print (str(k) + " current k")
                        if (k+1) % 3 == 1: #If it is a carbon comparison
                            if abs(abs_j - np.sqrt(New_AP[k][0]**2 + New_AP[k][1]**2 + New_AP[k][2]**2)) >= 0.75:
                                print ("success - the atom does not overlap with an existing C")
                            else:
                                print ("fail - the atom overlaps with an existing C")
                                break
                        else: #If it is hydrogen
                            if abs(abs_j - np.sqrt(New_AP[k][0]**2 + New_AP[k][1]**2 + New_AP[k][2]**2)) >= 0.5:
                                print("success - the atom does not overlap with an existing H")
                            else:
                                print ("fail - the atom overlaps with an existing H")
                                break
                        if (k+1) == c:
                            flag = False
                            New_AP.append(perturbed_j)
                else:
                    flag = False
                    New_AP.append(perturbed_j)
                    print ("Success, first C added")
        print ("Polymer fragment finished.")
        New_AP_all_chains.append(New_AP)

    return New_AP_all_chains

print (Randomise_PE(test, 1))




    
