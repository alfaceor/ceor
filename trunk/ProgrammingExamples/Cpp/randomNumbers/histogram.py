from numpy.random import normal
from numpy import loadtxt
# print gaussian_numbers
import matplotlib.pyplot as plt
from numpy.random import normal
import sys

if (len(sys.argv) > 1):
  datafilename=sys.argv[1]  
else:
  print 'usage: python histogram <datafilename>'
  sys.exit(1)

data = loadtxt(datafilename)
plt.hist(data)
plt.title("Gaussian Histogram")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.savefig('image.png')
plt.show()
