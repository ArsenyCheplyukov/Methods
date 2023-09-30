import numpy as np
import pandas
import matplotlib
import matplotlib.pyplot as plt

x = 1/np.array([0.164, 0.328, 0.656, 0.984, 1.312, 1.640])
y = np.array([0.448, 0.432, 0.421, 0.417, 0.414, 0.412])

a = 0.409974
b = 0.00644424

plt.plot(x, b*x+a)
plt.plot(x, y, 'ro')
plt.ylabel('some numbers')
plt.show()