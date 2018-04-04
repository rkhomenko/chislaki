from scipy.optimize import fsolve
import numpy as np

a = 3

def f(x):
    x1, x2 = x
    return ((x1 ** 2 + a ** 2) * x2 - a ** 3, \
            (x1 - a / 2) ** 2 + (x2 - a / 2) ** 2 - a ** 2)

x1, x2 = fsolve(f, (1, 1))
print(x1, x2, f((x1, x2)))
