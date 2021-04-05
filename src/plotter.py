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
    '''
    Brief :     Funcao que recebe o tempo inicial e o final, o numero de iteracoes
                e devolve o vetor com todos os pontos do tempo que serao analisados.
    Parameters: t0 - Tempo inicial,
                tf - Tempo final,
                n - Número de iteracoes.
    Returns:    time_array - Vetor com os valores de tempo a ser analisados.
    '''
    time_array = []
    h = (tf-t0)/n
    t = t0
    for _ in range(1,n+1):
        time_array.append(t)
        t += h
    return time_array

def plot_2d_1f(x, y, title, color_letter, labelx, labely, label_msg=""):
    '''
    Brief :     Funcao que plota um gráfico 2D com uma funcao linear y em funcao de x.
    Parameters: x - Vetor com valores do eixo x,
                y - Vetor com valores do eixo y,
                title - Título do gráfico,
                color_letter - Cor do gráfico,
                labelx - Nome do Eixo X
                labely - Nome do Eixo Y,
                label_msg - Legenda do gráfico (pode ser nao existente).
    '''
    if label_msg=="":
        plt.plot(x, y, color=color_letter) 
    else:
        plt.plot(x, y, label=label_msg, color=color_letter) 
        plt.legend()
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(title)
    plt.show()

def plot_multiple_graphs(x, y, title, n_rows, n_columns, legend, labelx, labely):
    '''
    Brief :     Funcao que plota múltiplos gráfico 2D com funcoes lineares
                y em funcao de x (com diferentes y e um mesmo x).
    Parameters: x - Vetor com valores do eixo x,
                y - Vetor com vetores de valores do eixo y,
                title - Título do gráfico,
                n_rows - número de linhas de gráficos,
                n_columns - número de colunas de gráficos,
                legend - Legenda para a funcao
                labelx - Nome do Eixo X
                labely - Nome do Eixo Y,
    '''
    for i in range(len(y)):
        plt.subplot(n_rows, n_columns, i+1)
        plt.plot(x,y[i])
        plt.legend([legend], loc=2, fontsize=6)
        plt.title("Aproximacao para x" + str(i) + "(t)")
        plt.xlabel(labelx)
        plt.ylabel(labely)
        plt.subplots_adjust(hspace=0.4, wspace = 0.5)
    plt.suptitle(title)
    plt.show()


def distance_graph(ts, ys, t, yexact, I, legend1, legend2, labelx, labely, title):
    '''
    Brief :     Funcao que plota dois gráficos em um plano 2D. Chama-se distance_graph
                pois apresenta visualmente a distância entre duas funcoes.
                Pode ser usado para comparar um y calculado e um y exato.
    Parameters: ts - Vetor com valores do tempo (eixo x),
                ys - Vetor com valores calculados do eixo y,
                yexact - Vetor com valores exatos do eixo y,
                I - Intervalo limite a ser analisado,
                legend1 - Legenda para a funcao ys,
                legend2 - Legenda para a funcao yexact,
                labelx - Nome do Eixo X,
                labely - Nome do Eixo Y,
                title - Título do gráfico,
    '''
    plt.plot(ts, ys, 'r')
    plt.plot(t, yexact, 'b')
    plt.xlim(I[0], I[1])
    plt.legend([legend1, 
                legend2], loc=2)
    plt.xlabel(labelx, fontsize=17)
    plt.ylabel(labely, fontsize=17)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plot_multiple_distance_graphs(x, y, exact_y, title, n_rows, n_columns, legend1, legend2, xlabel, ylabel):
    '''
    Brief :     Funcao que plota diversos gráficos de distância em só uma página.
                Pode ser usado para comparar um y calculado e um y exato.
    Parameters: x - Vetor com valores do eixo x,
                y - Vetor com vetores de valores do eixo y,
                exact_y - Vetor com vetores de valores exatos do eixo y,
                title - Título do gráfico,
                n_rows - número de linhas de gráficos,
                n_columns - número de colunas de gráficos,
                legend1 - Legenda para a funcao y,
                legend2 - Legenda para a funcao exact_y,
                labelx - Nome do Eixo X,
                labely - Nome do Eixo Y,
    '''
    for i in range(len(y)):
        plt.subplot(n_rows, n_columns, i+1)
        plt.plot(x,y[i])
        plt.plot(x,exact_y[i])
        plt.legend([legend1, legend2], loc=2, fontsize=6)
        plt.title("Aproximacao para x" + str(i) + "(t)")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.subplots_adjust(hspace=0.4, wspace = 0.5)
    plt.suptitle(title)
    plt.show()

def distance_graph_3(ts, xs, ys, zs, I, legend1, legend2, legend3, title):
    '''
    Brief :     Funcao que plota três funcoes em só um gráfico.
    Parameters: ts - Vetor com valores do tempo (eixo x),
                xs - Vetor com valores calculados de x para o eixo y,
                ys - Vetor com valores calculados de y para o eixo y,
                zs - Vetor com valores calculados de z para o eixo y,
                I - Intervalo limite a ser analisado,
                legend1 - Legenda para a funcao xs,
                legend2 - Legenda para a funcao ys,
                legend3 - Legenda para a funcao zs,
                title - Título do gráfico,
    '''
    plt.plot(ts, xs, "blue")
    plt.plot(ts, ys, "black")
    plt.plot(ts, zs, "red")
    plt.xlim(I[0], I[1])
    plt.legend([legend1, 
                legend2,
                legend3], loc=2)
    plt.xlabel("Tempo", fontsize=17)
    plt.ylabel("População", fontsize=17)
    plt.title(title)
    plt.tight_layout()
    plt.show()

'''
def pair_distance_graph_3(ts, xs, ys, zs, I, legend1, legend2, legend3, labelx, labely, titles):
    for i in range(2):
        plt.subplot(1,2,i+1)
        plt.plot(ts[i], xs[i], "blue")
        plt.plot(ts[i], ys[i], "black")
        plt.plot(ts[i], zs[i], "red")
        plt.title(titles[i])
        plt.xlim(I[0], I[1])
        plt.legend([legend1, 
                legend2,
                legend3], loc=2)
        plt.xlabel(labelx, fontsize=17)
        plt.ylabel(labely, fontsize=17)
        plt.tight_layout() 
    plt.show()
'''

def plot_3d_graph(xs, ys, zs, title, labelx, labely, labelz, legend):
    '''
    Brief :     Funcao que plota uma curva paramétrica em 3D
    Parameters: xs - Vetor com valores x do ponto,
                ys - Vetor com valores y do ponto,
                zs - Vetor com valores z do ponto,
                title - Título do gráfico,
                labelx - Nome do Eixo X,
                labely - Nome do Eixo Y,
                labelz - Nome do Eixo Z,
                I - Intervalo limite a ser analisado,
                legend - Legenda para a curva.              
    '''
    ax = plt.axes(projection="3d")
    ax.plot3D(xs,ys,zs)
    plt.title(title)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_zlabel(labelz)
    plt.legend([legend])
    plt.show()