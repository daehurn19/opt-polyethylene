import numpy as np
from Build_PE import Build_PE
from Build_randbox import Acceptability

def main(N): # Number of chains to create inside box.
    all_chain_L = []
    num_C = []
    AP_allchains = []
    AT_allchains = []
    #build N normal polyethylene chains
    for i in range(N):
        ran_num = np.random.randint(5,20)
        AP_allchains.append(Build_PE(ran_num)[0])
        AT_allchains.append(Build_PE(ran_num)[1])
        all_chain_L.append(AP_allchains[-1][-3][0])
        num_C.append(ran_num)
    print(str(N) + " standardised polyethylene generated.")

    # then build random box to contain chains (with PBC active)
    flag = True
    rand_box_dim = [] #format is a,b,c,alp,bet,gam
    count = 0
    while flag:
        acceptable = Acceptability(all_chain_L, num_C)

        if acceptable[0]:
            flag = False
            for i in range(6):
                rand_box_dim.append(acceptable[i+1])
            print("Acceptable box configuration found.")

    return (rand_box_dim)


    # then perturb each chain
    # then place each chain randomly in box


print(main(10))

