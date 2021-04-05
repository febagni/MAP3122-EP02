
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
from euler_forward import forwardEuler

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

    '''
    Brief : Essa funçao devolve os valores exatos de x dado um vetor de tempo "time_array".
    Parameters: time_array - Um vetor com os valores de tempo.
    Returns:    Um vetor com os valores exatos de x(t) para cada tempo fornecido
    '''

    exact_x = []
    for t in range(len(time_array)):
        exact_x.append(explicit_solution(time_array[t]))
    return exact_x

def plot_x (t,x,title):

    '''
    Brief : Essa funcao recebe o tempo t e os valores de x e plota os gráficos 
            necessários para a análise dos resultados
    Parameters: t - Vetor com todos os valores de tempo analisados
                x - Matriz com todos os valores (x1,x2,x3... xn) de x analisados para cada t
                title - título do gráficoo
    '''

    x = np.transpose(x)
    exact_x = np.transpose(get_exact_x(t))
    plot_multiple_graphs(t, x, title, 2, 2, "Valor Calculado x(t)", "t", "x") # Plota somente o gráfico com a solucoes calculadas
    plot_multiple_distance_graphs(t, x, exact_x, title, 2, 2 , "Valor Calculado x(t)", "Valor Exato x*(t)", "t", "x") # Plota o gráfico comparativo com as solucoes calculadas e exatas

def exercise1_test1():
    '''
    Brief : Essa funçao representa o script para a solucao do exercicio 1 - teste 1
    '''

    print("Exercício 1 - Teste 1: ")

    x0 = [1,1,1,-1]

    f = []
    
    # Esses appends foram calculados por meio da matriz A que foi fornecida no enunciado do exercício
    f.append(lambda t,x : -2*x[0] -1*x[1] -1*x[2] -2*x[3])
    f.append(lambda t,x : 1*x[0] -2*x[1] + 2*x[2] -1*x[3])
    f.append(lambda t,x : -1*x[0] -2*x[1] -2*x[2] -1*x[3])
    f.append(lambda t,x : 2*x[0] -1*x[1] + 1*x[2] -2*x[3])
    
    # Valores de n para os quais os dados serao analisados
    n = [20, 40, 80, 160, 320, 640]

    # Iterando sobre o vetor n, analisando valores para cada elemento do vetor
    for i in range (len(n)-1):
        [t_sol,x_t] = rk4system(x0, f, 0, 2, n[i])                  # Cálculo de RK4 para n_i
        plot_x(t_sol,x_t,"Gráficos de x(t) para n=" + str(n[i]))
        [t_aux,x_aux] = rk4system(x0, f, 0, 2, n[i + 1])            # Cálculo de RK4 para n_i+1
        erro = calc_error(t_sol, x_t, explicit_solution,)           # Cálculo do erro para n_i
        erro_aux = calc_error(t_aux, x_aux, explicit_solution)      # Cálculo do erro para n_i+1
        R = max(erro)/max(erro_aux)                                 # Cálculo de R
        print("O R para a iteracao i = " + str(i+1) + " é : " + str(R))
        plot_2d_1f(t_sol, erro, str("Gráfico de E_1(t) para n =" + str(n[i])), 'r', "t", "E_1", "Erro E_1(t)")
        print("O valor calculado x(Tf) para n =  " + str(n[i]) + " é : " + str(x_t[n[i]]))
        print("O valor exato x*(Tf) para n = " + str(n[i]) + " é : " + str(explicit_solution(t_sol[n[i]])))
        print()
        if i == len(n) - 2 : # Último caso: impressao dos dados para n = 640
            print("O valor calculado x(Tf) para n =  " + str(n[i+1]) + " é : " + str(x_aux[n[i+1]]))
            print("O valor exato x*(Tf) para n = " + str(n[i+1]) + " é : " + str(explicit_solution(t_aux[n[i+1]])))
            plot_x(t_aux,x_aux,"Gráficos de x(t) para n=" + str(n[i+1]))
            plot_2d_1f(t_aux, erro_aux, str("Gráfico de E(t) para n =" + str(n[i+1])), 'r', "t", "E_1", "Erro E_1(t)")

def exercise1_test2():
    '''
    Brief : Essa funçao representa o script para a solucao do exercicio 1 - teste 2
    '''
    
    print("Exercício 1 - Teste 2: ")
    u_0 = np.array([-8.79])
    f = []
    I = np.array([1.1, 3.0])
    n = 5000
    newton_iter_num = 7
    h = (I[1] - I[0])/n

    # Funcao f(t,x) para a qual será aplicada o metodo.
    f.append(lambda t,x : 2*t + (x-t*t)*(x-t*t))

    [ts, resposta_euler_implicito] = implicit_euler_system(u_0, f, I[0], I[1], n, newton_iter_num)

    #--- Calculates the exact solution, for comparison ---#
    dt = int((I[1] - I[0]) / h)
    t = [I[0]+i*h for i in range(dt+1)]
    yexact = []
    for i in range(dt+1):
        ye = t[i]*t[i] + 1/(1-t[i])
        yexact.append(ye)

    plot_2d_1f(ts, resposta_euler_implicito, str("Gráfico de Solução por Euler Implicito"), 'r', "t", "x", "solução pelo método x(t)")
    plot_2d_1f(t, yexact, str("Gráfico de Solução Explícita"), 'b', "t", "x", "solução exata explícita x*(t)")

    distance_graph(ts, resposta_euler_implicito, t, yexact, I, "Valor Calculado", "Valor Exato", "t", "x")

    #Falta cáculo de E_2 e plots 

    E_2 = [abs(yexact[i] - resposta_euler_implicito[i]) for i in range(len(yexact))]

    plot_2d_1f(ts, E_2, str("Gráfico de E_2(t)"), 'r', "t", "E_2", "erro E_2(t)")

    print("A solução pelo método de Euler implícito é: ",resposta_euler_implicito[n])
    print("A solução explícita é: ", yexact[n])

def main():
    '''
    Brief : Essa funçao eh a main do exercicio 1
    '''
    exercise1_test1()
    print()
    exercise1_test2()

main()