import numpy as np
import matplotlib.pyplot as plt

t = np.arange(-1.5, 1.6, 0.1)

plt.plot(t, t * t)
plt.plot(t, np.log(t + 2))
plt.savefig('figure1.png')
