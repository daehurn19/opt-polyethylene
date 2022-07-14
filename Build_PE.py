import numpy as np


def Build_PE(L): # L = length of chain (to be random by nature)
    #Function: Builds a model Polyethylene in the x-direction. Returns atomic position (XYZ) and element.

    #c-c-c angle
    cca = (180-109.5)*np.pi/180 /2

    #Define rotation matrix where normal plane lies in the z-axis. The polyethylene chain extends via x-axis.

    R = np.array([[np.cos(cca), np.sin(cca), 0],
                   [-np.sin(cca), np.cos(cca), 0],
                   [0,            0,           1]])

    #Define mirror matrix - H positions flip per chain basis.

    R90 = np.array([[0, 1, 0],
                   [-1, 0, 0],
                   [0, 0, 1]])

    #H rotation out of plane 60 degrees up/down.

    Rh = np.array([[1,            0,                 0],
                   [0, np.cos(np.pi/3), np.sin(np.pi/3)],
                   [0, -np.sin(np.pi/3), np.cos(np.pi/3)]])

    #Bond lengths
    CCL = 1.33
    CHL = 0.99

    #Output data
    atom_position = []
    atom_type = []


    #Set chain in x direction:
    chainvec = np.array([1, 0, 0])

    #Set initial -CH2 block:
    atom_position.append(np.array([0,0,0]))
    atom_type.append('C') #C = 1, H = 2

    atom_position.append(np.array([0,CHL*np.cos(np.pi/3), CHL*np.sin(np.pi/3)]))
    atom_type.append('H')

    atom_position.append(np.array([0, CHL*np.cos(np.pi/3), -CHL*np.sin(np.pi/3)]))
    atom_type.append('H')

    i = 3 #Counter for number of atoms
    while len(atom_position) < L*3:
        i += 1
        R = np.transpose(R) #fluctuating rotation matrix per CH2 subgroup.
        R90 = np.transpose(R90)
        atom_position.append(np.array(atom_position[-3] + CCL * np.dot(R, chainvec)))
        atom_type.append('C') #C

        i +=1 #H1
        R90_chainvec = np.dot(R90,chainvec)
        atom_position.append(np.array(atom_position[-1] + CHL * np.dot(Rh, R90_chainvec)))
        atom_type.append('H')

        i+=1 #H2
        R90_chainvec = np.dot(R90, chainvec)
        atom_position.append(np.array(atom_position[-2] + CHL * np.dot(np.transpose(Rh), R90_chainvec)))
        atom_type.append('H')

    return [atom_position, atom_type]


