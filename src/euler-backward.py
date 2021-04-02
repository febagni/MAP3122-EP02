import matplotlib.pyplot as plt
import numpy as np


def backwardEuler(f, yinit, x_range, h):
    m = len(yinit)
    n = int((x_range[-1] - x_range[0])/h)

    x = x_range[0]
    y = yinit

    xsol = np.empty(0)
    xsol = np.append(xsol, x)

    ysol = np.empty(0)
    ysol = np.append(ysol, y)

    for i in range(n):
        yprime = f(x+h, y)/(1+h)

        for j in range(m):
            y[j] = y[j] + h*yprime[j]

        x += h
        xsol = np.append(xsol, x)

        for r in range(len(y)):
            ysol = np.append(ysol, y[r])  # Saves all new y's

    return [xsol, ysol]


def myFunc(x, y):
    '''
    We define our ODEs in this function.
    '''
    dy = np.zeros((len(y)))
    dy[0] = 2*x + (y-(x*x))*(y-(x*x))
    return dy

n = 5000
h = (3-1.1)/n
x = np.array([1.1, 3.0])
yinit = np.array([-8.79])

[ts, ys] = backwardEuler(f=myFunc, yinit=yinit, x_range=x, h=h)

#--- Calculate the exact solution, for comparison ---#
t = [x[0]+i*h for i in range(n+1)]
yexact = []
for i in range(n+1):
    ye = t[i]*t[i] + 1/(1 - t[i])
    yexact.append(ye)

plt.plot(ts, ys, 'r')
plt.plot(t, yexact, 'b')
plt.xlim(x[0], x[1])
plt.legend(["Backward Euler method",
            "Exact solution"], loc=2)
plt.xlabel('x', fontsize=17)
plt.ylabel('y', fontsize=17)
plt.tight_layout()
plt.show()