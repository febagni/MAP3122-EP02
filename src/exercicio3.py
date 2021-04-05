#!/usr/bin/env python3

##############################################################
# File      :   exercicio3.py
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
from plotter import plot_2d_1f, distance_graph, plot_3d_graph, distance_graph_3
from euler_backward import implicit_euler_system
from euler_forward import forwardEuler
from rk4 import rk4system

def init():
    '''
    Brief :     Funcao inicializadora que fornece valores comuns às análises do exercício
    Returns:    Número de iteracoes "n",
                Tempo inicial "t0",
                Vetor de constantes "B" (para o sistema de EDOs),
                Matriz de Constantes "A" (para o sistema de EDOs).
    '''
    n = 5000
    t0 = 0
    B = np.array([1.0, 1.0, -1.0])
    A = np.array([[0.0010, 0.001, 0.015],
                  [0.0015, 0.001, 0.001],
                  [0.000, -0.0005, 0.0]])
    return n, t0, B, A

def get_f(A,B):
    '''
    Brief :     Funcao que fornece as funcoes f do sistema de EDOs, dadas a matriz A
                e o vetor B.
    Parameters: Matriz de Constantes "A",
                Vetor de Constantes "B".
    Returns:    Vetor com as funçoes do Sistema.
    '''
    f = []
    #f.append(lambda t, u : u[0]*(1.0 - 0.001*u[0] - 0.001*u[1] - 0.0015*u[2]))
    #f.append(lambda t, u : u[1]*(1.0 - 0.0015*u[0] - 0.001*u[1] - 0.001*u[2]))
    #f.append(lambda t, u : u[2]*(-1.0 + 0.0033*u[0] - (-0.0005)*u[1] - 0*u[2]))
    #return f
    for i in range(len(A)):
        f.append(lambda t, u, i=i: u[i]*(B[i] - A[i][0]*u[0] - A[i][1]*u[1] - A[i][2]*u[2]))
    return f

def divide_sol(sol):
    '''
    Brief :     Funcao que recebe as respostas de um sistema de EDOs com 3 variáveis
                e devolve os valores em três vetores distintos x, y e z.
    Parameters: Matriz de solucoes "sol".
    Returns:    Vetores com respostas x, y e z
    '''
    xs = np.transpose(sol)[0]
    ys = np.transpose(sol)[1]
    zs = np.transpose(sol)[2]
    return xs, ys, zs

def ex3_methods(method):
    '''
    Brief :     Funcao que calcula a solucao para o problema proposto pelo exercício 3,
                para diferentes tipos de métodos de resoluçao de EDO. Além de calcular 
                a solucao, também plota os gráficos para a devida análise.
    Parameters: method - Nome do método a ser utilizado
    '''
    n, t0, B, A = init()

    u_0 = np.array([500., 500., 10.])
    tf = np.array([100., 500., 2000.])
    alpha = np.array([0.001, 0.002, 0.0033, 0.0036, 0.005, 0.0055]) 

    for k in range(len(alpha)):
        A[2][0] = -alpha[k]
        tf_iter = tf[int(k/2)]
        f = get_f(A, B)
        I = [t0, tf_iter]

        if method == "rk4":
            [ts, sol] = rk4system(u_0, f, t0, tf_iter, n)
            method_name = "Runge-Kutta 4"
        elif method == "euler_forward":
            [ts, sol] = forwardEuler(f, u_0, I, n)
            method_name = "Euler Explícito"
        else : print("Error: Method not found")
        xs, ys, zs = divide_sol(sol)

        time_title = "T_f = " + str(tf_iter)
        alpha_title = "alpha = " + str(alpha[k])
        title_complement = str(time_title + " e " + alpha_title + " com " + method_name + "\n")

        plot_3d_graph(xs, ys, zs, "Retrato de fase para " + title_complement, "Coelhos", "Lebres", "Raposas", "Retrato de Fase")

        plot_2d_1f(xs, ys, "Retrato de fase : " + title_complement + "Coelhos e Lebres", "red", "Coelhos", "Lebres", "Retrato de Fase")
        plot_2d_1f(xs, zs, "Retrato de fase : " + title_complement + "Coelhos e Raposas", "red", "Coelhos", "Raposas", "Retrato de Fase")
        plot_2d_1f(ys, zs, "Retrato de fase : " + title_complement + "Lebres e Raposas", "red", "Lebres", "Raposas", "Retrato de Fase")
            
        distance_graph_3(ts, xs, ys, zs, I, "Coelhos", "Lebres,", "Raposas", "População ao longo do tempo: \n" + title_complement)


def ex3_1():
    '''
    Brief :     Funcao que Realiza os cálculos pedidos para a funçao de Runge Kutta 4
                e para o Euler explícito.
    '''
    print("Exercício 3.1")
    ex3_methods(method="rk4")
    ex3_methods(method="euler_forward")

def ex3_sensibility_test(method):
    '''
    Brief :     Funcao que realiza o teste de sensibilidade às diferencas das condicoes
                iniciais de um sistema de EDOs. O teste pode ser feito em diferentes 
                métodos. Além disso, plota os gráficos a serem analisados.
    Parameters: method - Nome do método de resolucao de EDO a ser utilizado.
    '''
    print("Exercício 3 - Teste de Sensibilidade")

    caso1 = np.array([37. ,75. ,137.])
    caso2 = np.array([37. ,74. ,137.])

    n, t0, B, A = init()
    tf = np.array([400])
    alpha = np.array([0.005])

    A[2][0] = -alpha[0]
    tf_iter = tf[0]
    f = get_f(A, B)
    I = [t0, tf_iter]
    
    if method == "rk4":
        [ts1, sol1] = rk4system(caso1, f, t0, tf_iter, n)
        [ts2, sol2] = rk4system(caso2, f, t0, tf_iter, n)
        method_name = "Runge-Kutta 4"
    elif method == "euler_forward":
        [ts1, sol1] = forwardEuler(f, caso1, I, n)
        [ts2, sol2] = forwardEuler(f, caso2, I, n)
        method_name = "Euler Explícito"
    else : print("Error: Method not found")
    print("| Caso: 1 | Método: " + method_name + " | Raposas: " + str(sol1[n][2]) + " | Lebres: " + str(sol1[n][1]) + " | Coelhos: " + str(sol1[n][0]) + " |")
    print("| Caso: 2 | Método: " + method_name + " | Raposas: " + str(sol2[n][2]) + " | Lebres: " + str(sol2[n][1]) + " | Coelhos: " + str(sol2[n][0]) + " |")
    print()
    print("As raposas cresceram: " + str(((sol2[n][2]-sol1[n][2])/sol1[n][2])*100) + "% , as Lebres cresceram: " + str(((sol2[n][1]-sol1[n][1])/sol1[n][1])*100) + "%, e os coelhos cresceram: " + str(((sol2[n][0]-sol1[n][0])/sol1[n][0])*100) + "%.")
    print()
    xs1, ys1, zs1 = divide_sol(sol1)
    xs2, ys2, zs2 = divide_sol(sol2)

    ts = [ts1, ts2]
    xs = [xs1, xs2]
    ys = [ys1, ys2]
    zs = [zs1, zs2]

    title_complement = " para " + method_name
    titles = ["Caso com 75 lebres" + title_complement, "Caso com 74 lebres" + title_complement]

    for i in range(2):
        distance_graph_3(ts[i], xs[i], ys[i], zs[i], I, "Coelhos", "Lebres,", "Raposas", titles[i])

def ex3_2():
    '''
    Brief :     Funcao que realiza o teste de sensibilidade para os dois casos
                propostos pelo exercício.
    '''
    ex3_sensibility_test("rk4")
    ex3_sensibility_test("euler_forward")

def main():
    '''
    Brief : Essa funcao eh a main do exercicio 3
    '''
    ex3_1()
    print()
    ex3_2()

main()