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
from plotter import plot_2d_1f, distance_graph
from euler_backward import implicit_euler_system
from euler_forward import forwardEuler
from rk4 import rk4system

def init():
    x_0 = 1.5
    y_0 = 1.5
    I = np.array([0, 10.0])
    return x_0, y_0, I

def ex2_1():
    print("Exercício 2.1: ")
    x_0, y_0, I = init()
    u_0 = np.array([x_0, y_0])
    f = []
    n = 5000
    h = (I[1] - I[0])/n
    # u[0] == x e u[1] == y
    f.append(lambda t,u : ((2/3)*u[0]) - ((4/3)*u[0]*u[1])) #=dx 
    f.append(lambda t,u : u[0]*u[1] - u[1]) #=dy

    [ts, ys] = forwardEuler(f, u_0, I, h)

    print(ys)
    
    x_sol = ys[0]
    y_sol = ys[1]

    #plot_2d_1f(ts, x_sol, str("Gráfico de Solução por Euler Implicito"), 'b', "solução pelo método")
    #plot_2d_1f(ts, y_sol, str("Gráfico de Solução por Euler Implicito"), 'r', "solução pelo método")

    distance_graph(ts, y_sol, ts, x_sol, I, "raposas", "coelhos")

    plot_2d_1f(x_sol, y_sol, str("Gráfico de Retrato de fase Coelhos X Raposas"), 'g', "retrato de fase")

def ex2_2():
    print("Exercício 2.2: ")
    x_0, y_0, I = init()
    u_0 = np.array([x_0, y_0])
    f = []
    n = 1000
    newton_iter_num = 7
    #h = (I[1] - I[0])/n
    # u[0] == x e u[1] == y
    f.append(lambda t,u : ((2/3)*u[0]) - ((4/3)*u[0]*u[1])) #=dx 
    f.append(lambda t,u : u[0]*u[1] - u[1]) #=dy

    [ts, resposta_euler_implicito] = implicit_euler_system(u_0, f, I[0], I[1], n, newton_iter_num)

    print(resposta_euler_implicito)
    
    x_sol = np.transpose(resposta_euler_implicito)[0]
    y_sol = np.transpose(resposta_euler_implicito)[1]

    #plot_2d_1f(ts, x_sol, str("Gráfico de Solução por Euler Implicito"), 'b', "solução pelo método")
    #plot_2d_1f(ts, y_sol, str("Gráfico de Solução por Euler Implicito"), 'r', "solução pelo método")

    distance_graph(ts, y_sol, ts, x_sol, I, "raposas", "coelhos")

    plot_2d_1f(x_sol, y_sol, str("Gráfico de Retrato de fase Coelhos X Raposas"), 'g', "retrato de fase")

def ex2_4():
    x_0, y_0, I = init()
    u_0 = np.array([x_0, y_0])
    f = []
    n = 500
    f.append(lambda t,u : ((2/3)*u[0]) - ((4/3)*u[0]*u[1])) #=dx 
    f.append(lambda t,u : u[0]*u[1] - u[1]) #=dy
    [ts, resposta_rk4] = rk4system(u_0,f,I[0],I[1],n)
    x_sol = np.transpose(resposta_rk4)[0]
    y_sol = np.transpose(resposta_rk4)[1]

    distance_graph(ts, y_sol, ts, x_sol, I, "raposas", "coelhos")

    plot_2d_1f(x_sol, y_sol, str("Gráfico de Retrato de fase Coelhos X Raposas"), 'g', "retrato de fase")

def main ():
    ex2_2()

main()