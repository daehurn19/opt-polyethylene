import numpy as np


def Build_PE(L, CCL, CHL): # L = length of chain (to be random by nature)
    """
    Builds an approximate low-energy polyethylene structure with conventional bonding angles, using matrix calculations
    propagated along the x-axis.

    :param L: Length of each chain - expected values should be even to avoid PBC problems.
    :param CCL: C-C bond length in Angstroms - defaults to 1.54A
    :param CHL: C-H bond length in Angstroms - defaults to 1.10A

    :return: Atomic positions of a typical Polyethylene chain, centered at the origin; chain heads towards x-axis.
    """

    # c-c-c angle
    cca = (180-109.5)*np.pi/180 / 2

    # Define rotation matrix where normal plane lies in the z-axis. The polyethylene chain extends via x-axis #
    # and displaced on the y-axis.
    R_C = np.array([[np.cos(cca), np.sin(cca), 0],
                   [-np.sin(cca), np.cos(cca), 0],
                   [0,            0,           1]])

    # Define mirror matrix - H positions flip per chain basis.
    R_mirrorH = np.array([[0, 1, 0],
                         [-1, 0, 0],
                         [0, 0, 1]])

    # H rotation out of plane 60 degrees up/down.
    R_planeH = np.array([[1,            0,                 0],
                         [0, np.cos(np.pi/3), np.sin(np.pi/3)],
                         [0, -np.sin(np.pi/3), np.cos(np.pi/3)]])

    # Output data
    atom_position = []
    atom_type = []

    # Set chain in x direction:
    chainvec = np.array([1.0, 0.0, 0.0])

    # Set initial -CH2 block:
    atom_position.append(np.array([0.0,0.0,0.0]))
    atom_type.append('C')  # C = 1, H = 2

    atom_position.append(np.array([0.0,CHL*np.cos(np.pi/3), CHL*np.sin(np.pi/3)]))
    atom_type.append('H')

    atom_position.append(np.array([0.0, CHL*np.cos(np.pi/3), -CHL*np.sin(np.pi/3)]))
    atom_type.append('H')

    i = 3  # Counter for number of atoms
    while len(atom_position) < L*3:
        i += 1

        atom_position.append(np.array(atom_position[-3] + CCL * np.dot(R_C, chainvec)))
        atom_type.append('C')

        i += 1  # H1
        R_mH_chainvec = np.dot(R_mirrorH, chainvec)
        atom_position.append(np.array(atom_position[-1] + CHL * np.dot(R_planeH, R_mH_chainvec)))
        atom_type.append('H')

        i += 1  # H2
        R_mH_chainvec = np.dot(R_mirrorH, chainvec)
        atom_position.append(np.array(atom_position[-2] + CHL * np.dot(np.transpose(R_planeH), R_mH_chainvec)))
        atom_type.append('H')

        R_C = np.transpose(R_C)  # fluctuating rotation matrix per CH2 subgroup.
        R_mirrorH = np.transpose(R_mirrorH)

    return [atom_position, atom_type]


