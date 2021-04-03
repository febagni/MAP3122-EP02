
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
from plotter import plot_scatter, get_time_array
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

def exercise1_test1():
    '''
    Brief : Essa funçao representa o script para a solucao do exercicio 1 - teste 1
    '''
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
        resultado = rk4system(x0, A, 0, 2, n[i])
        resultado_aux = rk4system(x0, A, 0, 2, n[i + 1])
        erro = calc_error(resultado, explicit_solution, 0, 2, n[i])
        erro_aux = calc_error(resultado_aux, explicit_solution, 0, 2, n[i + 1])
        xf.append(resultado[n[i]])
        R.append(max(erro)/max(erro_aux))
        plot_scatter(time_array, erro, str("Gráfico de E(t) para n =" + str(n[i])))
        if i == len(n) - 2 : plot_scatter(get_time_array(0,2,n[i+1]), erro_aux, str("Gráfico de E(t) para n =" + str(n[i+1])))
    print(R)

def exercise1_test2():

    u_0 = np.array([-8.79])
    f = []
    f.append(lambda t,x : 2*t + (x-t*t)*(x-t*t))

    resposta_euler_implicito = implicit_euler_system(u_0, f, 1.1, 3.0, 5000, 7)

    #h = 0.2
    #x = np.array([1.1, 3.0])

    #[ts, ys] = forwardEuler(f=myFunc, yinit=u_0, x_range=x, h=h)


    #--- Calculates the exact solution, for comparison ---#
    #dt = int((x[-1] - x[0]) / h)
    #t = [x[0]+i*h for i in range(dt+1)]
    #yexact = []
    #for i in range(dt+1):
    #    ye = t[i]*t[i] + 1/(1-t[i])
    #    yexact.append(ye)

    #distance_graph(ts, ys, t, yexact, x)

    #Falta cáculo de E_2 e plots 

    print("A solução pelo método de Euler implícito é: ",resposta_euler_implicito[5000])

def main():
    '''
    Brief : Essa funçao eh a main do exercicio 1
    '''
    exercise1_test1()
    exercise1_test2()

#if __name__ == "__main__":
#    main()
main()