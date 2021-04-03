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
from plotter import plot_2d_1f, distance_graph
from euler_backward import implicit_euler_system
from euler_forward import forwardEuler
from rk4 import rk4system

def init():
    n = 5000
    t0 = 0
    B = [1.0,1.0,-1.0]
    A = [[0.0010, 0.001, 0.015],
         [0.0015, 0.001, 0.001],
         [0.000, -0.0005, 0.0]]
    return n, t0, B, A

def ex3_1():
    n,t0,B,A = init()
    u = [500,500,10]
    tf = [100,500,2000]
    alpha = [0.001, 0.002, 0.0033, 0.0036, 0.005, 0.0055]
    f = []
    for k in range(len(alpha)):
        A[2][0] = alpha[k]
        tf_iter = tf[int(k/2)]
        h = (tf_iter - t0)/n
        for i in range(len(A)):
            f.append(lambda t,u : u[i]*(B[0] - A[i][0]*u[0] - A[i][1]*u[1] - A[i][2]*u[2]))
        I = [t0, tf_iter]
        [ts_euler, sol_euler] = forwardEuler(f,u,I,h)
        [ts_rk, sol_rk] = rk4system(u,f,t0,tf_iter,n)

def main():
    ex3_1()

main()