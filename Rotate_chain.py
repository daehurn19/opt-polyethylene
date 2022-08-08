import numpy as np

def Rotator(AP_chain):
    """
    Randomly rotates a Polyethylene chain to increase variation and non-uniformity of initial structure.

    :param AP_chain: Single PE chain within box.

    :return: Randomly rotated Polyethylene chain - 0 to 360 degrees.
    """

    n_AP_chain = []

    # random angle
    rand_angle = np.random.uniform(0, 2*np.pi, size=1)

    R_yz = np.array([[1,            0,                 0],
                   [0, np.cos(rand_angle), np.sin(rand_angle)],
                   [0, -np.sin(rand_angle), np.cos(rand_angle)]])

    # For each atom's atomic position in the chain, apply the rotation matrix to the coordinates
    for i in AP_chain:
        n_AP_chain.append(np.dot(i, R_yz))

    return n_AP_chain




