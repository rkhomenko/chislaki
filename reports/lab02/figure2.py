import numpy as np
import matplotlib.pyplot as plt

a = 3

t = np.arange(0, 2 * np.pi + 0.1, 0.1)

x = np.arange(-2, 5.1, 0.1)

plt.plot(a * np.cos(t) + a / 2, a * np.sin(t) + a / 2)
plt.plot(x, a ** 3 / (x ** 2 + a ** 2))

plt.xlim((-2, 5))
plt.ylim((-2, 5))

plt.savefig('figure2.png')
