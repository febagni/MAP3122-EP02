#!/usr/bin/env python3

##############################################################
# File      :   exercicio1.py
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

import matplotlib.pyplot as plt
 
def get_time_array(t0,tf,n):
    time_array = []
    h = (tf-t0)/n
    t = t0
    for _ in range(1,n+1):
        time_array.append(t)
        t += h
    return time_array

def plot_scatter(x, y, title):
    plt.scatter(x, y, label='erros', color='r') #diminuir tamanho dos pontos
    plt.xlabel('t')
    plt.ylabel('E')
    plt.title(title)
    plt.legend()
    plt.show()

#3d graph method

#line diff graph method

def distance_graph(ts, ys, t, yexact, x):
    plt.plot(ts, ys, 'r')
    plt.plot(t, yexact, 'b')
    plt.xlim(x[0], x[1])
    plt.legend(["Method solution", 
                "Exact solution"], loc=2)
    plt.xlabel('x', fontsize=17)
    plt.ylabel('y', fontsize=17)
    plt.tight_layout()
    plt.show()