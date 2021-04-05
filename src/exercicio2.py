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
    '''
    Brief :     Funcao inicializadora que fornece valores comuns às análises do exercício
    Returns:    Valores iniciais u_0 / Intervalo I = [t0, tf]
    '''
    x_0 = 1.5
    y_0 = 1.5
    I = np.array([0., 10.0])
    return x_0, y_0, I

def ODE_solver(n, method):
    '''
    Brief : Essa funçao resolve a EDO para diferentes métodos, aplicando "n" iteracoes.
    Parameters: n - Número de iteracoes a serem chamadas,
                method - Nome do método a ser utilizado.
    Returns:    ts - Valores de tempo onde a solucao é analisada;
                ys - Valores da soluçao da EDO calculados
    '''
    x_0, y_0, I = init()
    u_0 = np.array([x_0, y_0])
    f = []

    f.append(lambda t,u : ((2/3)*u[0]) - ((4/3)*u[0]*u[1])) #=dx 
    f.append(lambda t,u : u[0]*u[1] - u[1]) #=dy

    if method=="euler_forward":
        [ts, ys] = forwardEuler(f, u_0, I, n)
    elif method=="euler_backward":
        [ts, ys] = implicit_euler_system(u_0, f, I[0], I[1], n, newton_iter_num=7)
    elif method=="rk4":
        [ts, ys] = rk4system(u_0, f, I[0], I[1], n)
    else:
        print("Error: Method not found")

    return [ts, ys]

def plot_graphs(ts, ys, method):

    '''
    Brief : Essa funçao plota os gráficos em funcao do tempo ts, e também cria o retrato
            de fase dos valores de ys.
    Parameters: ts - Valores de tempo onde a solucao é analisada;
                ys - Valores da soluçao da EDO calculados
    '''
    
    _, _, I = init()
    
    if method=="euler_forward":
        method_name = "Euler Explícito"
    elif method=="euler_backward":
        method_name = "Euler Implícito"
    elif method=="rk4":
        method_name = "Runge-Kutta 4"
    else:
        print("Error: Method not found")
    
    x_sol = np.transpose(ys)[0]
    y_sol = np.transpose(ys)[1]

    distance_graph(ts, y_sol, ts, x_sol, I, "Raposas", "Coelhos", "tempo", "população", str("Tamanho da população ao longo do tempo, por " + method_name))

    plot_2d_1f(x_sol, y_sol, str("Gráfico de Retrato de fase Coelhos X Raposas, por " + method_name), 'g', "coelhos", "raposas", "retrato de fase")

def ex2_1():
    '''
    Brief : Essa funçao resolve o que é requerido no exercício 2.1:
            Resoluçao da EDO utilizando o método de Euler explícito
    '''
    print("Exercício 2.1: ")

    [ts, ys] = ODE_solver(n=500, method="euler_forward")

    plot_graphs(ts, ys, method="euler_forward")

def ex2_2():
    '''
    Brief : Essa funçao resolve o que é requerido no exercício 2.2
            Resoluçao utilizando o método de Euler Implícito
    '''
    print("Exercício 2.2: ")

    [ts, ys] = ODE_solver(n=1000, method="euler_backward")

    plot_graphs(ts, ys, method="euler_backward")
    

def ex2_3():
    '''
    Brief : Essa funçao resolve o que é requerido no exercício 2.3
            Cálculo das diferenças entre o método Implícito e Explícito
    '''
    print("Exercício 2.3: ")
    _, _, I = init()
    n = [250, 500, 1000, 2000, 4000]
    for n_i in n:
        [ts_forward, ys_forward] = ODE_solver(n=n_i, method="euler_forward")
        [ts_backward, ys_backward] = ODE_solver(n=n_i, method="euler_backward")

        E = [np.asarray(j) - np.asarray(ys_forward[i]) for i,j in enumerate(ys_backward)] 

        E_x_sol = np.transpose(E)[0]
        E_y_sol = np.transpose(E)[1]

        distance_graph(ts_backward, E_y_sol, ts_forward, E_x_sol, I, "E_y", "E_x", "t", "E", str("Erro em função do tempo para n = " + str(n_i)))
    
def ex2_4():
    '''
    Brief : Essa funçao resolve o que é requerido no exercício 2.4
            Resoluçao da EDO utilizando o método de Runge-Kutta de ordem 4.
    '''
    print("Exercício 2.4: ")

    [ts, ys] = ODE_solver(n=1000, method="rk4")

    plot_graphs(ts, ys, method="rk4")

def main():
    '''
    Brief : Essa funçao eh a main do exercicio 2
    '''
    ex2_1()
    ex2_2()
    ex2_3()
    ex2_4()

main()