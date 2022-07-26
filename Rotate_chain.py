import numpy as np

def Rotator(AP_chain):
    #Function: Rotates a chain a random amount between 0 and 360 degrees.
    #The chain propagates in the x-axis, the chain should be rotated along the y-z plane.

    #Output:
    new_AP_chain = []

    rand_angle = np.random.uniform(0, 2*np.pi, size=1)


    R_yz = np.array([[1,            0,                 0],
                   [0, np.cos(rand_angle), np.sin(rand_angle)],
                   [0, -np.sin(rand_angle), np.cos(rand_angle)]])

    #For each atom's atomic position in the chain, apply the rotation matrix to the coordinates
    for i in AP_chain:
        new_AP_chain.append(np.dot(i, R_yz))

    return new_AP_chain




