#!/usr/bin/env python3

##############################################################
# File      :   metodos.py
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

def rk4iter(u, t, h, f) :
    '''
    Brief :Essa funçao aplica o algoritmo Runge Kutta em uma unica iteracao
    Parâmetros: x - valor a ser iterado,
                t - valor de t a iterar,
                h - passo da funçao,
                f - funcao f para aplicar Runge Kutta.
    '''
    K1 = np.zeros(len(f))
    K2 = np.zeros(len(f))
    K3 = np.zeros(len(f))
    K4 = np.zeros(len(f))

    for i in range(len(f)):
        K1[i] = h*f[i](t,u)
    for i in range(len(f)):
        K2[i] = h*f[i](t + h/2,u + K1/2)
    for i in range(len(f)):
        K3[i] = h*f[i](t + h/2,u + K2/2)
    for i in range(len(f)):
        K4[i] = h*f[i](t+h,u + K3)
    return u + (K1 + 2*K2 +2*K3 + K4)/6

def rk4system(u, f, t0, tf, n):
    '''
    Brief:  Essa funçao aplica o algoritmo Runge Kutta resolvendo um SPVI
            sendo o sistema de forma linear. O sistema pode ser fornecido
            por meio da matriz A.
    Parâmetros: u - valores iniciais, 
                A - matriz do sistema linear,
                t0 - valor inicial de t, 
                tf - valor final de t,
                n - número de divisoes entre t0 e tf.
    '''
    newx = np.zeros(len(f))
    x = u.copy()
    rk4values = []
    rk4values.append(x)
    tsol = []
    h = (tf-t0)/n
    t = t0
    tsol.append(t)
    for _ in range (1,n+1):
        newx = rk4iter(x,t,h,f)
        t = t + h
        x = newx.copy()
        tsol.append(t)
        rk4values.append(x)
    return [tsol, np.array(rk4values)]

def calc_error(t_solution, calc_solution, explicit_solution):
    '''
    Brief : Essa funçao é responsável por calcular o erro de cada iteraçao
            tendo as respostas explícitas.
    Parameters: calc_solution - Os valores de x(t) calculados por algum método iterativo,
                explicit_solution - Funçoes de t que sao a solucao explicita do problema,
                t0 / tf - tempo inicial e tempo final, respecitvamente,
                n - número de iteracoes realizadas.
    Returns:    Um vetor com os valores de erro para cada iteracao.
    '''
    errIter = []
    for k in range (len(calc_solution)):
        errIter.append(max([abs(explicit_solution(t_solution[k])[i] - calc_solution[k][i]) for i in range (len(calc_solution[k]))]))
    return errIter