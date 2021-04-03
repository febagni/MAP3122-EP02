
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

import numpy as np
from plotter import plot_2d_1f, get_time_array, distance_graph, plot_multiple_graphs, plot_multiple_distance_graphs
from rk4 import rk4system, calc_error
from euler_backward import implicit_euler_system

def explicit_solution(t):

    '''
    Brief : Essa funçao fornece a solucao explicita para o item 1 do primeiro exercicio,
            dado determinado valor de t.
    Parameters: t - valor de tempo.
    Returns:    Um vetor com os valores de explicitos de x(t) para t
    '''

    x_explicit = np.array([
                 np.exp(-t)*np.sin(t) + np.exp(-3*t)*np.cos(3*t),
                 np.exp(-t)*np.cos(t) + np.exp(-3*t)*np.sin(3*t),
                 -np.exp(-t)*np.sin(t) + np.exp(-3*t)*np.cos(3*t),
                 -np.exp(-t)*np.cos(t) + np.exp(-3*t)*np.sin(3*t)])
    return x_explicit

def get_exact_x(time_array):
    exact_x = []
    for t in range(len(time_array)):
        exact_x.append(explicit_solution(time_array[t]))
    return exact_x

def plot_x (t,x,title):
    x = np.transpose(x)
    exact_x = np.transpose(get_exact_x(t))
    plot_multiple_distance_graphs(t, x, exact_x, title, 2, 2 , "Valor calculado", "Valor Exato")

def exercise1_test1():
    '''
    Brief : Essa funçao representa o script para a solucao do exercicio 1 - teste 1
    '''

    print("Exercício 1 - Teste 1: ")

    x0 = [1,1,1,-1]

    f = []

    f.append(lambda t,x : -2*x[0] -1*x[1] -1*x[2] -2*x[3])
    f.append(lambda t,x : 1*x[0] -2*x[1] + 2*x[2] -1*x[3])
    f.append(lambda t,x : -1*x[0] -2*x[1] -2*x[2] -1*x[3])
    f.append(lambda t,x : 2*x[0] -1*x[1] + 1*x[2] -2*x[3])

    n = [20, 40, 80, 160, 320, 640]

    for i in range (len(n)-1):
        [t_sol,x_t] = rk4system(x0, f, 0, 2, n[i]) 
        plot_x(t_sol,x_t,"Gráficos de x(t) para n=" + str(n[i]))
        [t_aux,x_aux] = rk4system(x0, f, 0, 2, n[i + 1])
        erro = calc_error(t_sol, x_t, explicit_solution,)
        erro_aux = calc_error(t_aux, x_aux, explicit_solution)
        R = max(erro)/max(erro_aux)
        print("O R para a iteracao i = " + str(i) + " é : " + str(R))
        plot_2d_1f(t_sol, erro, str("Gráfico de E(t) para n =" + str(n[i])), 'r', "Erro")
        if i == len(n) - 2 : 
            plot_x(t_aux,x_aux,"Gráficos de x(t) para n=" + str(n[i+1]))
            plot_2d_1f(t_aux, erro_aux, str("Gráfico de E(t) para n =" + str(n[i+1])), 'r', "Erro")

def exercise1_test2():
    print("Exercício 1 - Teste 2: ")
    u_0 = np.array([-8.79])
    f = []
    I = np.array([1.1, 3.0])
    n = 5000
    newton_iter_num = 7
    h = (I[1] - I[0])/n

    f.append(lambda t,x : 2*t + (x-t*t)*(x-t*t))

    [ts, resposta_euler_implicito] = implicit_euler_system(u_0, f, I[0], I[1], n, newton_iter_num)

    #--- Calculates the exact solution, for comparison ---#
    dt = int((I[1] - I[0]) / h)
    t = [I[0]+i*h for i in range(dt+1)]
    yexact = []
    for i in range(dt+1):
        ye = t[i]*t[i] + 1/(1-t[i])
        yexact.append(ye)

    plot_2d_1f(ts, resposta_euler_implicito, str("Gráfico de Solução por Euler Implicito"), 'r', "solução pelo método")
    plot_2d_1f(t, yexact, str("Gráfico de Solução Explícita"), 'b', "solução exata explícita")

    distance_graph(ts, resposta_euler_implicito, t, yexact, I, "Valor Calculado", "Valor Exato")

    #Falta cáculo de E_2 e plots 

    E_2 = [abs(yexact[i] - resposta_euler_implicito[i]) for i in range(len(yexact))]

    plot_2d_1f(ts, E_2, str("Gráfico de E_2(t)"), 'r', "erros")

    print("A solução pelo método de Euler implícito é: ",resposta_euler_implicito[5000])

def main():
    '''
    Brief : Essa funçao eh a main do exercicio 1
    '''
    exercise1_test1()
    print()
    exercise1_test2()

#if __name__ == "__main__":
#    main()
main()