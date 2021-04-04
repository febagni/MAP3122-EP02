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
from mpl_toolkits import mplot3d
 
def get_time_array(t0,tf,n):
    time_array = []
    h = (tf-t0)/n
    t = t0
    for _ in range(1,n+1):
        time_array.append(t)
        t += h
    return time_array

def plot_2d_1f(x, y, title, color_letter, label_msg):
    plt.plot(x, y, label=label_msg, color=color_letter) 
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.legend()
    plt.show()

def plot_multiple_graphs(x, y, title, n_rows, n_columns):
    for i in range(len(y)):
        plt.subplot(n_rows, n_columns, i+1)
        plt.plot(x,y[i])
        plt.title("Aproximacao para x" + str(i) + "(t)")
        plt.subplots_adjust(hspace=0.4, wspace = 0.5)
    plt.suptitle(title)
    plt.show()

def plot_multiple_distance_graphs(x, y, exact_y, title, n_rows, n_columns, legend1, legend2):
    for i in range(len(y)):
        plt.subplot(n_rows, n_columns, i+1)
        plt.plot(x,y[i])
        plt.plot(x,exact_y[i])
        plt.legend([legend1, legend2], loc=2, fontsize=6)
        plt.title("Aproximacao para x" + str(i) + "(t)")
        plt.subplots_adjust(hspace=0.4, wspace = 0.5)
    plt.suptitle(title)
    plt.show()

#3d graph method

#line diff graph method

def distance_graph(ts, ys, t, yexact, I, legend1, legend2):
    plt.plot(ts, ys, 'r')
    plt.plot(t, yexact, 'b')
    plt.xlim(I[0], I[1])
    plt.legend([legend1, 
                legend2], loc=2)
    plt.xlabel('x', fontsize=17)
    plt.ylabel('y', fontsize=17)
    plt.tight_layout()
    plt.show()

def plot_3d_graph(xs,ys,zs):
    ax = plt.axes(projection="3d")
    ax.plot3D(xs,ys,zs)
    plt.show()