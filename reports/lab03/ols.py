import numpy as np
import matplotlib.pyplot as plt

x = np.array([-1.0, 0.0, 1.0, 2.0, 3.0, 4.0])
y = np.array([0.86603, 1.0, 0.86603, 0.5, 0.0, -0.5])
t = np.arange(-1.0, 4.1, 0.1)

plt.plot(x, y, label='function')
plt.plot(x, 0.89232248 - 0.29131943 * x, label='lls')
plt.plot(t, 0.94748879 - 0.043071036 * t - 0.082749464 * t * t, label='sls')
plt.legend()
plt.show()
