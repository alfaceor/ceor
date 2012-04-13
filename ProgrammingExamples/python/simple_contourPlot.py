import numpy as np
import matplotlib.pyplot as plt
X = np.arange(-3, 3, 1)
Y = np.arange(-3, 3, 1)
X,Y = np.meshgrid(X,Y)
R = np.sqrt(X**2+Y**2)
# plt.contour(X,Y,R)
plt.contourf(X,Y,R)
plt.show()
