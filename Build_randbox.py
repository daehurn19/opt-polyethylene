import numpy as np



def Build_Randombox(all_chain_L, num_C, CCL):
    #all_chain_L, a list containing the x-lengths of all chains to be placed into box.
    #Function: Builds a fixed box with a=length of chain,b,c,alpha=90,beta=90,gamma=90
    #The box has a boundary condition in all 3 axis. Volume of box scales linearly with N.
    #Note that the scaling of a is independent of N - this is because all C backbones are defined to the x-axis.

    volume =  len(all_chain_L) * 18 * sum(num_C) #General rule of thumb for organic molecules is 18 x no non-H atoms -  estimated
                                       # volume in Angstrom3

    #set a vector as equal to average chain length of each carbon polymer fragment
    a = sum(all_chain_L)/len(all_chain_L) + 0.816*CCL # L carbon propogated in x axis



    #Randomly perturbs the box length in b and c direction by up to 50% in either direction
    b = np.random.uniform(0.8,1.25) * np.sqrt((volume/a))
    c = np.random.uniform(0.8,1.25) * np.sqrt((volume/a))

    #Randomly perturbs the box angles by up to 50% high or low
    alpha = np.random.uniform(0.8, 1.25) * 90
    beta = np.random.uniform(0.8, 1.25) * 90
    gamma = np.random.uniform(0.8, 1.25) * 90

    pseudo_V = a * b * c * np.sqrt(1 + (2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)) - (np.cos(alpha) ** 2) - (np.cos(beta) ** 2) - (
                    np.cos(gamma) ** 2))


    if 1.1 * volume > pseudo_V > 0.90 * volume:
        return True, a, b, c, alpha, beta, gamma
    else:

        return False, a, b, c, alpha, beta, gamma







