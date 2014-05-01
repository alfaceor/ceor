import os
filename=os.path.splitext(os.path.basename('66_22_1000100010001_0.04_20000000_0.001_1.0_0.13_-1.0_100.dat'))[0]
print filename
print '-'*10
rev, njobs, chain, temperature, ttime, dt, epsi, q, Ec, print_each =filename.split('_')

print chain
