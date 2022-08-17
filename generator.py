from main import generate
import numpy as np
from sys import argv

pressure = float(argv[1])
runs = int(argv[2])

for i in range(runs):
    N = np.random.randint(1,5)
    
    if N == 1:
        L = np.random.randint(1,5)*2
    elif N == 2:
        L = np.random.randint(1,3)*2
    elif N == 3 or N == 4:
        L = 2

    if __name__ == "__main__":
        generate(N, L, pressure=pressure)

