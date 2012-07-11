import os
filename=os.path.splitext(os.path.basename('/home/alfaceor/rev_chain_temp.dat'))[0]
print filename
filename.split('_')
for i in filename.split('_'):
  print i
