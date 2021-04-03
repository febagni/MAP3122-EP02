#!/usr/bin/env python3

##############################################################
# File      :   euler_forward.py
# Project   :   MAP3122 - EP02 - Metodos numericos para 
#                                resolucao de EDOs
# Date      :   April/2021#!/usr/bin/env python3
##############################################################
# Author    :   Felipe Bagni
# Email     :   febagni@usp.br
# NUSP      :   11257571
##############################################################
# Author    :   Gabriel Yugo Kishida
# Email     :   gabriel.kishida@usp.br
# NUSP      :   11257647
##############################################################

import numpy as np
from plotter import distance_graph

'''
The Euler's forward method or explicit Euler's method
or Euler-Cauchy method is formulated as:
    y[i+1] = y[i] + h * f(x[i], y[i]),
    
where f(x[i], y[i]) is the differential equation evaluated
at x[i] and y[i].
'''


def forwardEuler(f, yinit, x_range, h):
    '''
    This function/module performs the forward Euler method steps.
    '''
    m = len(yinit) # Number of ODEs
    n = int((x_range[-1] - x_range[0])/h) # Number of sub-intervals
    
    x = x_range[0] # Initializes variable x
    y = yinit # Initializes variable y
    
    xsol = np.empty(0) # Creates an empty array for x
    xsol = np.append(xsol, x) # Fills in the first element of xsol

    ysol = []
    ysol.append(y) # Fills in the initial conditions


    for _ in range(n):

        for j in range(m):
            y[j] = y[j] + h*f[j](x, y) # Eq. (8.2)
            
        x += h # Increase x-step
        xsol = np.append(xsol, x) # Saves it in the xsol array
        
        for r in range(len(y)):
            ysol.append(y[r]) # Saves all new y's 
            
    return [xsol, ysol]


#########################################################################################

# passar para exercicio-#.py
'''
def myFunc(x, y):
'''
    #We define our ODEs in this function
'''
    dy = np.zeros((len(y)))
    dy[0] = 3*(1+x) - y[0]

    return dy


f = []
f.append(lambda t,u : 3*(1 + t) - u[0])  #=dx 

h = 0.2
x = np.array([1.0, 2.0])
yinit = np.array([4.0]) #np.array([4.0, 1.0])


[ts, ys] = forwardEuler(f=f, yinit=yinit, x_range=x, h=h)


#--- Calculates the exact solution, for comparison ---#
dt = int((x[-1] - x[0]) / h)
t = [x[0]+i*h for i in range(dt+1)]
yexact = []
for i in range(dt+1):
    ye = 3*t[i] + np.exp(1-t[i])
    yexact.append(ye)

distance_graph(ts, ys, t, yexact, x, "sol implementada", "sol exata")
'''