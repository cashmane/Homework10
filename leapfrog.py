import numpy as np
import matplotlib.pyplot as plt

def f(r):
    return r[1]**2 - r[0] - 5

def leapfrog(initial_r, initial_t, h, tmax):
    r1 = initial_r
    r2 = initial_r +0.5*h*f(initial_r)
    xList = []
    while initial_t < tmax-h:
        r1 = r1 + h*f(r2)
        r2 = r2 + h*f(r1)
        xList.append(r1)
        initial_t += h
    return np.array(xList)

if __name__ == '__main__':
    tmax = 50
    initial_r = np.array([1, 0])
    times = np.arange(0, 50, 0.001)
    xs = leapfrog(initial_r, 0, 0.001, tmax)
    plt.plot(times, xs[:,0])
    plt.show()
