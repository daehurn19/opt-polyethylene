import numpy as np
from Build_PE import Build_PE
from Build_randbox import Build_Randombox
from Rand_Perturbation import Randomise_Position
import ase, ase.io, ase.visualize
from ase.calculators.lammpslib import LAMMPSlib
from ase.constraints import ExpCellFilter
from ase.optimize import LBFGSLineSearch
from Rotate_chain import Rotator
import os
from os.path import exists


def generate(N, L, runs, CCL=1.54, CHL=1.10, rot=True, pressure=0):
    """
    Generates a random polyethylene initial structure in a random, scaling box, and relaxes it using 
    an AIREBO potential. Periodic boundary conditions are active.
    
    :param N: Number of individual chains to be generated.
    :param L: Length of each chain - expected values should be even to avoid PBC problems.
    :param runs: Number of runs to attempt.
    :param CCL: C-C bond length in Angstroms - defaults to 1.54A
    :param CHL: C-H bond length in Angstroms - defaults to 1.10A
    :param rot: Applies random rotation matrix to each chain - defaults to True
    :param pressure: Pressure of system in GPa - defaults to 0 GPa
     
    :return: A polymer .res file to use with AIRSS or equivalent.
    """
    
    for run in range(runs):
        
        chain_x = []
        AP_all = []  # Atom_positions of all N chains, sorted by generation history
        AT_all = []  # Atom_types of all N chains, sorted by generation history
        N_carbons = []

        # Build N normal polyethylene chains
        for i in range(N):
            AP_all.append(Build_PE(L, CCL, CHL)[0])
            AT_all.append(Build_PE(L, CCL, CHL)[1])
            chain_x.append(AP_all[-1][-3][0])
            N_carbons.append(L)
            
        print(str(N) + " standardised polyethylene generated.")

        # format is a,b,c,alp,bet,gam
        box_dimensions = []

        # build random box to contain chains
        flag = True
        while flag:
            trial_box = Build_Randombox(chain_x, N_carbons, CCL)

            if trial_box[0]:
                flag = False
                for i in range(6):
                    box_dimensions.append(trial_box[i+1])
                print("Acceptable box configuration found.")

        n_AP_all = []

        print (box_dimensions)

        # Check for rotation of polymer
        if rot:
            for i in AP_all:
                n_AP_all.append(Rotator(i))

        else:
            n_AP_all = AP_all

        # then place each chain randomly in box # account for PBC
        # Can be done by making sure there is a minimum 1.5A distance from each chain,
        # and all chains stay within the confines of the box dimensions in b and c direction.

        b, c = box_dimensions[1], box_dimensions[2]
        n_AP_all = Randomise_Position(n_AP_all, b, c)
        print("All polyethylene fragments randomly placed within box.")

        # Now, flatten all lists for output.
        # Use ASE to do this easily.
        # We have both box_dimensions (dimensions of random box), and n_AP_all
        # (coordinates of randomised polyethylene fragments).

        print("Flattening elements to write...")
        elements = ''
        for s in AT_all:
            for t in s:
                elements += t
        print("Success.")

        print("Flattening atom positions to write...")
        flat_AP = []
        for s in n_AP_all:  # n_AP_all = list of np.arrays of 3 elements
            for t in s:
                flat_AP.append(t)
        print ("Success.")

        atoms = ase.Atoms(elements, np.array(flat_AP), cell=box_dimensions, pbc=True)
        atoms.set_cell(box_dimensions, scale_atoms=True)

        print(atoms.get_cell_lengths_and_angles)

        ase.visualize.view(atoms)
        print("Visualising...")

        # Optimisation step
        to_eV_per_A3 = 1.0 / 160.2177  # from GPa
        # Relax cell using GAP

        cmds = ["pair_style airebo 3.0 1 1",
                "pair_coeff * * CH.airebo C H"]

        lammps = LAMMPSlib(lmpcmds=cmds, log_file='test.log', keep_alive=True)

        atoms.calc = lammps
        ecf = ExpCellFilter(atoms, scalar_pressure=pressure)
        optimizer = LBFGSLineSearch(ecf)
        try:
            optimizer.run(1e-4, 1000)
        except:
            print("Failed to find minimum structure.")

        # If converged, output final structure and fake castep file

        if optimizer.converged():
            print("Converged.")
            c = 1
            flag = True
            file_name = r'poly{}-out.xyz'.format(c)
            while flag:
                if exists(file_name):
                    c += 1
                    file_name = r'poly{}-out.xyz'.format(c)
                else:
                    flag = False
                    ase.io.write(r'poly{}-out.xyz'.format(c), atoms, 'extxyz')
            volume = atoms.get_volume()
            pv = pressure * to_eV_per_A3 * volume
            enthalpy = atoms.get_potential_energy() + pv
            with open(r'poly{}.castep'.format(c), 'w') as f:
                f.write("Current cell volume = {:25.15f} A**3\n".format(volume))
                f.write("*  Pressure:   {:25.15f}\n".format(pressure))
                f.write("Python: Final Enthalpy     = {:25.15f} eV\n".format(enthalpy))

            os.system('cabal xyze cell < poly{}-out.xyz > poly{}-out.cell'.format(c, c))
            os.system('castep2res poly{} > poly{}.res;'.format(c, c))

    return 0



rand_struct = generate(1, 10, 2) # 3 chains of L = 5
