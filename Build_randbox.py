import numpy as np

def Build_Randombox(chain_x, N_carbons, CCL):
    """
    Builds a random box that scales with the amount of carbons in the box. Checks that this random box has
    acceptable parameters (will not lead to errors/problems due to assertion/logic).

    :param chain_x: x-values of last carbon atom in each chain.
    :param N_carbons: Number of carbons total of all chains.
    :param CCL: C-C bond length in angstroms.
    
    :return: True, and all lattice parameters IF the box is an acceptable size. Otherwise, returns False.
    """

    # General rule of thumb for organic molecules is 18 x no non-H atoms -  estimated volume in Angstrom3
    volume = len(chain_x) * 18 * sum(N_carbons) 

    # set a vector as equal to average chain length of each carbon polymer fragment
    a = sum(chain_x)/len(chain_x) + 0.943*CCL # L carbon propagated in x-axis



    # Randomly perturbs the box length in b and c direction by up to 25% in either direction
    b = np.random.uniform(0.75, 1.25) * np.sqrt((volume/a))
    c = np.random.uniform(0.75, 1.25) * np.sqrt((volume/a))

    # Randomly perturbs the box angles by up to 25% high or low
    alpha = np.random.uniform(0.75, 1.25) * 90
    beta = np.random.uniform(0.75, 1.25) * 90
    gamma = np.random.uniform(0.75, 1.25) * 90

    pseudo_V = a * b * c * np.sqrt(1 + (2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)) - (np.cos(alpha) ** 2) - (np.cos(beta) ** 2) - (
                    np.cos(gamma) ** 2))

    if 1.1 * volume > pseudo_V > 0.90 * volume:
        return True, a, b, c, alpha, beta, gamma
    else:

        return False, a, b, c, alpha, beta, gamma







