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
from plotter import plot_2d_1f, distance_graph, plot_3d_graph, distance_graph_3, pair_distance_graph_3
from euler_backward import implicit_euler_system
from euler_forward import forwardEuler
from rk4 import rk4system

def init():
    n = 5000
    t0 = 0
    B = np.array([1.0, 1.0, -1.0])
    A = np.array([[0.0010, 0.001, 0.015],
                  [0.0015, 0.001, 0.001],
                  [0.000, -0.0005, 0.0]])
    return n, t0, B, A

def get_f(A,B):
    f = []
    #f.append(lambda t, u : u[0]*(1.0 - 0.001*u[0] - 0.001*u[1] - 0.0015*u[2]))
    #f.append(lambda t, u : u[1]*(1.0 - 0.0015*u[0] - 0.001*u[1] - 0.001*u[2]))
    #f.append(lambda t, u : u[2]*(-1.0 + 0.0033*u[0] - (-0.0005)*u[1] - 0*u[2]))
    #return f
    for i in range(len(A)):
        f.append(lambda t, u, i=i: u[i]*(B[i] - A[i][0]*u[0] - A[i][1]*u[1] - A[i][2]*u[2]))
    return f

def ex3_1():
    ex3_methods(method="rk4")
    ex3_methods(method="euler_forward")

def divide_sol(sol):
    xs = np.transpose(sol)[0]
    ys = np.transpose(sol)[1]
    zs = np.transpose(sol)[2]
    return xs, ys, zs

def ex3_methods(method):
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
        elif method == "euler_forward":
            [ts, sol] = forwardEuler(f, u_0, I, n)
        else : print("Error: Method not found")
        xs, ys, zs = divide_sol(sol)

        plot_3d_graph(xs, ys, zs)

        plot_2d_1f(xs, ys, "Retrato de fase: Coelhos e Lebres", "red")
        plot_2d_1f(xs, zs, "Retrato de fase: Coelhos e Raposas", "red")
        plot_2d_1f(ys, zs, "Retrato de fase: Lebres e Raposas", "red")
            
        distance_graph_3(ts, xs, ys, zs, I, "Coelhos", "Lebres,", "Raposas")

def ex3_sensibility_test(method):

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
    elif method == "euler_forward":
        [ts1, sol1] = forwardEuler(f, caso1, I, n)
        [ts2, sol2] = forwardEuler(f, caso2, I, n)
    else : print("Error: Method not found")
    print("| Caso: 1 | Método: " + method + " | Raposas: " + str(sol1[n][2]) + " | Lebres: " + str(sol1[n][1]) + " | Coelhos: " + str(sol1[n][0]) + " |")
    print("| Caso: 2 | Método: " + method + " | Raposas: " + str(sol2[n][2]) + " | Lebres: " + str(sol2[n][1]) + " | Coelhos: " + str(sol2[n][0]) + " |")
    print()
    print("As raposas cresceram: " + str(((sol2[n][2]-sol1[n][2])/sol1[n][2])*100) + "% , as Lebres cresceram: " + str(((sol2[n][1]-sol1[n][1])/sol1[n][1])*100) + "%, e os coelhos cresceram: " + str(((sol2[n][0]-sol1[n][0])/sol1[n][0])*100))
    print()
    xs1, ys1, zs1 = divide_sol(sol1)
    xs2, ys2, zs2 = divide_sol(sol2)

    ts = [ts1, ts2]
    xs = [xs1, xs2]
    ys = [ys1, ys2]
    zs = [zs1, zs2]
    titles = ["Caso com 75 lebres", "Caso com 74 lebres"]

    pair_distance_graph_3(ts,xs,ys,zs,I,"Coelhos","Lebres","Raposas", titles)
    #distance_graph_3(ts, xs, ys, zs, I, "Coelhos", "Lebres,", "Raposas")

def ex3_2():
    ex3_sensibility_test("rk4")
    ex3_sensibility_test("euler_forward")

def main():
    ex3_1()
    #ex3_2()

main()