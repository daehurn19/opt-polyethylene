import numpy as np



def Build_Randombox(all_chain_L, num_C):
    #all_chain_L, a list containing the x-lengths of all chains to be placed into box.
    #Function: Builds a fixed box with a=length of chain,b,c,alpha=90,beta=90,gamma=90
    #The box has a boundary condition in all 3 axis. Volume of box scales linearly with N.
    #Note that the scaling of a is independent of N - this is because all C backbones are defined to the x-axis.

    volume =  len(all_chain_L) * 18 * sum(num_C) #General rule of thumb for organic molecules is 18 x no non-H atoms -  estimated
                                       # volume in Angstrom3

    #set a vector as equal to average chain length of each carbon polymer fragment
    a = sum(all_chain_L)/len(all_chain_L) + 1.258 # L carbon propogated in x axis



    #Randomly perturbs the box length in b and c direction by up to 50% in either direction
    b = np.sqrt((volume/a)) + 1.5
    c = np.sqrt((volume/a)) + 1.5

    #Randomly perturbs the box angles by up to 50% high or low
    alpha = np.random.uniform(0.5, 1.5) * 90
    beta = np.random.uniform(0.5, 1.5) * 90
    gamma = np.random.uniform(0.5, 1.5) * 90


    return a,b,c,alpha,beta,gamma, volume


def Acceptability(all_chain_L, num_C):
    #Function: The generated lattice/box is acceptable if this is passed.
    #Checks if scalar triple product of the vectors is roughly equal to
    #the volume of the cell (~10% acceptable deviation).
    vals = Build_Randombox(all_chain_L, num_C)
    a,b,c,alpha,beta,gamma, volume = vals[0],vals[1],vals[2],vals[3],vals[4],vals[5], vals[6]
    pseudo_V = a*b*c*np.sqrt(1+(2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))-(np.cos(alpha)**2) - (np.cos(beta)**2)-(np.cos(gamma)**2))
    if pseudo_V < 1.1*volume and pseudo_V > 0.90*volume:
        return True, a, b, c, alpha, beta, gamma
    else:
        
        return False, a, b, c, alpha, beta, gamma





