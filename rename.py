import os

for i in range(400):
    os.rename(r'relaxation\poly-out{}.xyz'.format(i+1), r'relaxation\poly{}-out.xyz'.format(i+1))