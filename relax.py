from ase import Atoms
from ase.constraints import ExpCellFilter
import ase.io
from ase.optimize import LBFGSLineSearch
import numpy as np
from ase.calculators.lammpslib import LAMMPSlib
import sys, os
from os.path import exists


#------------------------------------------------

# Collect command line input and load random structure
def relaxation(atoms, pressure=0):
        to_eV_per_A3 = 1.0 / 160.2177  # from GPa
        # Relax cell using GAP

        cmds = ["pair_style airebo 3.0 1 1",
                "pair_coeff * * CH.airebo C H"]

        lammps = LAMMPSlib(lmpcmds=cmds, log_file='test.log', keep_alive=True)

        atoms.calc = lammps
        ecf = ExpCellFilter(atoms, scalar_pressure=0)
        optimizer = LBFGSLineSearch(ecf)
        try:
            optimizer.run(1e-3, 500)
        except:
            print("fails")
        # If converged, output final structure and fake castep file

        if optimizer.converged():
            print ("Converged.")
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
            pv = pressure*to_eV_per_A3 * volume
            enthalpy = atoms.get_potential_energy() + pv
            with open(r'poly{}.castep'.format(c), 'w') as f:
                f.write("Current cell volume = {:25.15f} A**3\n".format(volume))
                f.write("*  Pressure:   {:25.15f}\n".format(pressure))
                f.write("Python: Final Enthalpy     = {:25.15f} eV\n".format(enthalpy))

            os.system('cabal xyze cell < poly{}-out.xyz > poly{}-out.cell'.format(c, c))
            os.system('castep2res poly{} > poly{}.res;'.format(c,c))
