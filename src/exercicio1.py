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

import numpy as np
from plotter import plot_graph, get_time_array
from rk4 import rk4system, calc_error

x0 = [1,1,1,-1]
 
A = np.array([[-2,-1,-1,-2],
              [ 1,-2, 2,-1],
              [-1,-2,-2,-1],
              [ 2,-1, 1,-2]])

n = [20, 40, 80, 160, 320, 640]

xf = []
R = []

for i in range (len(n)-1):
    time_array = get_time_array(0,2,n[i])
    resultado = rk4system(x0,A,0,2,n[i])
    resultado_aux = rk4system(x0,A,0,2,n[i + 1])
    erro = calc_error(resultado,0,2,n[i])
    erro_aux = calc_error(resultado_aux,0,2,n[i + 1])
    xf.append(resultado[n[i]])
    R.append(max(erro)/max(erro_aux))
    plot_graph(time_array, erro, str("Gráfico de E(t) para n =" + str(n[i])))
    if i == len(n) - 2 : plot_graph(get_time_array(0,2,n[i+1]), erro_aux, str("Gráfico de E(t) para n =" + str(n[i+1])))

print(R)