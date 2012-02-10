from numpy.random import normal
from numpy import loadtxt
gaussian_numbers = normal(size=1000)
# print gaussian_numbers
import matplotlib.pyplot as plt
from numpy.random import normal
gaussian_numbers = normal(size=1000)
data = loadtxt("data.txt")
plt.hist(data)
plt.title("Gaussian Histogram")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()
