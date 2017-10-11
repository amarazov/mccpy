import numpy as np
from matplotlib import pyplot as plt

x, y = np.loadtxt('exampleFile.csv',
                  unpack=True,
                  delimiter=',')

plt.plot(x, y)

plt.title('Epic Info $Omega$')
plt.ylabel('Y axis')
plt.xlabel('X axis')

plt.show()
