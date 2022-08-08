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


def generate(N, L, runs, CCL=1.54, CHL=1.10, rot=False, pressure=0):

    # Input: N = Number of chains to be created
    #        L = C Length of the polyethylene chain - must be multiples of 2.
    #        runs = Number of generation + relaxation attempts.
    #        CCL = C-C bond length (Angstroms) - defaults to 1.54
    #        CHL = C-H bond length (Angstroms) - defaults to 1.10
    #        rot = True/False for random rotation of each polyethylene chain, Defaults to False
    #        pressure = Pressure of the system - defaults to 0 (in GPa)

    for run in range(runs):

        #Output variables
        all_chain_L = []
        AP_allchains = [] #Atom_positions of all N chains, sorted by generation history
        AT_allchains = [] #Atom_types of all N chains, sorted by generation history
        num_C = []

        #Build N normal polyethylene chains
        for i in range(N):
            AP_allchains.append(Build_PE(L, CCL, CHL)[0])
            AT_allchains.append(Build_PE(L, CCL, CHL)[1])
            all_chain_L.append(AP_allchains[-1][-3][0])
            num_C.append(L)
        print(str(N) + " standardised polyethylene generated.")

        # then build random box to contain chains
        flag = True
        rand_box_dim = [] #format is a,b,c,alp,bet,gam
        while flag:
            acceptable = Build_Randombox(all_chain_L, num_C, CCL)

            if acceptable[0]:
                flag = False
                for i in range(6):
                    rand_box_dim.append(acceptable[i+1])
                print("Acceptable box configuration found.")

        new_AP_allchains = []

        print (rand_box_dim)
        norm_box_dim = rand_box_dim
        print (norm_box_dim)
        norm_box_dim[3] = 90
        norm_box_dim[4] = 90
        norm_box_dim[5] = 90

        #Check for rotation of polymer
        if rot:
            for i in AP_allchains:
                new_AP_allchains.append(Rotator(i))


        else:
            new_AP_allchains = AP_allchains


        # then place each chain randomly in box # account for PBC
        # Can be done by making sure there is a minimum 1.5A distance from each chain,
        # and all chains stay within the confines of the box dimensions in b and c direction.

        b, c = norm_box_dim[1], norm_box_dim[2]
        new_AP_allchains = Randomise_Position(new_AP_allchains, b, c)
        print("All polyethylene fragments randomly placed within box.")

        #Now, flatten all lists for output.
        # Use ASE to do this easily.
        #We have both rand_box_dim (dimensions of random box), and New_AP_allchains (coordiantes of randomised polethylene fragments)
        print ("Flattening elements to write...")
        elements = ''
        for s in AT_allchains:
            for t in s:
                elements += t
        print ("Success.")

        print("Flattening atom positions to write...")
        flat_AP = []
        for s in new_AP_allchains: #New_ap_allchains = list of np.arrays of 3 elements
            for t in s:
                flat_AP.append(t)
        print ("Success.")


        atoms = ase.Atoms(elements, np.array(flat_AP), cell=norm_box_dim, pbc=True)
        atoms.set_cell(rand_box_dim, scale_atoms=True)

        ase.visualize.view(atoms)
        print("Visualising...")

        # Optimisation step
        to_eV_per_A3 = 1.0 / 160.2177  # from GPa
        # Relax cell using GAP

        cmds = ["pair_style airebo 3.0 1 1",
                "pair_coeff * * CH.airebo C H"]

        lammps = LAMMPSlib(lmpcmds=cmds, log_file='test.log', keep_alive=True)

        atoms.calc = lammps
        ecf = ExpCellFilter(atoms, scalar_pressure=0)
        optimizer = LBFGSLineSearch(ecf)
        try:
            optimizer.run(1e-4, 1000)
        except:
            print("fails")
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



rand_struct = generate(1, 4, 2) # 3 chains of L = 5
