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
from plotter import plot_2d_1f, distance_graph, plot_3d_graph
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
    n, t0, B, A = init()
    u_0 = np.array([500, 500, 10])
    tf = np.array([100, 500, 2000])
    
    alpha = np.array([0.001, 0.002, 0.0033, 0.0036, 0.005, 0.0055])

    for k in range(len(alpha)):
        f = []
        A[2][0] = -alpha[k]
        tf_iter = tf[int(k/2)]
        h = (tf_iter - t0)/n
        f = get_f(A, B)
        I = [t0, tf_iter]
        #[ts_euler, sol_euler] = forwardEuler(f,u,I,h)
        [ts_rk, sol_rk] = rk4system(u_0,f,t0,tf_iter,n)
        
        xs = np.transpose(sol_rk)[0]
        ys = np.transpose(sol_rk)[1]
        zs = np.transpose(sol_rk)[2]
        #plot_2d_1f(ts_rk,xs,"xs","red","a")
        #plot_2d_1f(ts_rk,ys,"ys","red","a")
        #plot_2d_1f(ts_rk,zs,"zs","red","a")
        plot_3d_graph(xs,ys,zs)

def main():
    ex3_1()

main()