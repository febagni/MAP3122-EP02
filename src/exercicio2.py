#!/usr/bin/env python3

##############################################################
# File      :   exercicio2.py
# Project   :   MAP3122 - EP02 - Metodos numericos para 
#                                resolucao de EDOs
# Date      :   April/2021
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
from plotter import plot_2d_1f
from euler_backward import implicit_euler_system

def init():
    x_0 = 1.5
    y_0 = 1.5
    I = np.array([0, 10.0])
    return x_0, y_0, I

def myFunc(x, y):
        '''
        We define our ODEs in this function
        '''
        lmbda = 2/3
        alfa = 4/3
        beta = 1
        gama = 1
        dy = np.zeros((len(y)))
        dy[0] = 3*(1+x) - y[0]

        return dy

def ex2_2():
    print("Exercício 2.2: ")
    x_0, y_0, I = init()
    u_0 = np.zeros(2)
    u_0[0] = x_0
    u_0[1] = y_0
    f = []
    n = 500
    newton_iter_num = 7
    h = (I[1] - I[0])/n
    # u[0] == x e u[1] == y
    f.append(lambda t,u : ((2/3)*u[0]) - ((4/3)*u[0]*u[1])) #=dx 
    f.append(lambda t,u : u[0]*u[1] - u[1]) #=dy

    [ts, resposta_euler_implicito] = implicit_euler_system(u_0, f, I[0], I[1], n, newton_iter_num)
    
    plot_2d_1f(ts, resposta_euler_implicito[0], str("Gráfico de Solução por Euler Implicito"), 'r', "solução pelo método")
    

    
def main ():
    ex2_2()


main()