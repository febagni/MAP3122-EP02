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

def forwardEuler(f, yinit, I, n):
    '''
    This function/module performs the forward Euler method steps.
    '''
    m = len(yinit) # Number of ODEs

    h = (I[1] - I[0])/n
    
    x = I[0] # Initializes variable x
    y = yinit.copy() # Initializes variable y
    
    xsol = np.empty(0) # Creates an empty array for x
    xsol = np.append(xsol, x) # Fills in the first element of xsol

    ysol = []
    ysol.append(y.copy()) # Fills in the initial conditions


    for _ in range(n):

        for j in range(m):
            y[j] = y[j] + h*f[j](x, y) # Euler-forward eq.
            
        x += h # Increase x-step
        xsol = np.append(xsol, x) # Saves it in the xsol array
        ysol.append(y.copy()) # Saves all new y's 
    return [xsol, ysol]