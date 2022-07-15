import random as rn
import numpy as np
import Build_PE as BP

def Build_Fixedbox(N, atom_position): #N = number of PE molecules to be placed inside.
    #Function: Builds a fixed box with a=length of chain,b,c,alpha=90,beta=90,gamma=90
    #The box has a boundary condition in all 3 axis.
    print ("sdfs")

    volume = len(atom_position)/3 * 20 #General rule of thumb for organic molecules is 20 x no non-H atoms -  estimated
                                       # volume in Angstrom3

    #set a vector as equal to x-position of last carbon
    a = atom_position[-3][0]
    b = np.sqrt(volume/a) # b = c
    alpha = 90 #alpha = beta = gamma

    return a,b,alpha

print (Build_Fixedbox(1, BP.Build_PE(30)[0]))




