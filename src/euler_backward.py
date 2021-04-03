#!/usr/bin/env python3

##############################################################
# File      :   euler_backward.py
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
import numpy as np


def partial_derivative(t, u, partial_in, f, step):
    '''
    Brief : Essa funçao calcula a derivada parcial de uma das funcoes de f (sendo
            f um vetor de funcoes f1, f2, ...). A derivada parcial é derivada na
            variavel de indice "partial_in" e é feita de forma numérica (para um step
            pequeno).
    Parâmetros: t - valor de tempo,
                u - variaveis do sistema,
                partial_in - indice (de "u") da variavel sobre a qual será derivada f
                f - vetor de funcoes, uma das quais será derivada,
                step - O tamanho do passo que define a precisao da derivada numerica.
    Retorna:    Valor da derivada parcial para os determinados valores de "input".
    '''
    step_u = u.copy()
    step_u[partial_in] += step
    return (f(t,step_u) - f(t,u))/step

def newton_iter(t, u, f, h, last_u):
    '''
    Brief : Essa funcao calcula uma iteracao do método de aproximacao de Newton para
            obter U_k+1 - com o objetivo de se apicar o método implícito de Euler.
    Parâmetros: t - valor de tempo,
                u - U_k+1 (antes da iteracao) variaveis do sistema,
                f - vetor de funcoes, uma das quais será derivada,
                h - tamanho do passo temporal percorrido a cada iteracao,
                last_u - U_k variáveis do sistema.
    Retorna:    Valor de U_k+1 após a iteracao de Newton
    '''
    n = len(u)
    jacobian = np.identity(n)
    G = np.zeros(len(u))
    new_t = t + h
    for i in range(n):
        for j in range(n):
            jacobian[i][j] -= h*partial_derivative(new_t, u, j, f[i], 0.001)
    inv_jacobian = np.linalg.inv(jacobian)
    for i in range(len(u)):
        G[i] = u[i] - h * f[i](t,u) - last_u[i]
    return u - np.matmul(inv_jacobian, G) # == new_u

def implicit_euler_iter(u, t, h, f, newton_iter_num):
    '''
    Brief : Essa funcao calcula uma iteracao do método de aproximacao de Newton para
            obter U_k+1 - com o objetivo de se apicar o método implícito de Euler.
    Parâmetros: u - U_k (antes da iteracao) variaveis do sistema,
                t - valor de tempo,
                h - tamanho do passo temporal percorrido a cada iteracao,
                f - vetor de funcoes, uma das quais será derivada,
                newton_iter - números de iteracao de newton que serao feitas.
    Retorna:    Valor de U_k+1 para a iteracao de Euler Implícito
    '''
    newton_u = u.copy()
    euler_u = u.copy()
    for _ in range(newton_iter_num) :
        newton_u = newton_iter(t, newton_u, f, h, u)
    for i in range(len(u)):
        euler_u[i] = u[i] + h * f[i](t, newton_u)
    return euler_u

def implicit_euler_system(u, f, t0, tf, n, newton_iter_num):
    '''
    Brief : Essa funcao calcula iteracoes do método de Euler implícito para um
            sistema de funcoes.
    Parâmetros: u - U_0 valores iniciais do sistema,
                f - vetor de funcoes que criam o sistema,
                t0 e tf - tempos iniciais e finais da iteracao,
                n - número de iteracoes do método de Euler a serem realizadas
                newton_iter_num - números de iteracao de newton que serao feitas.
    Retorna :   Valores de u para cada passo temporal calculado.
    '''
    u_values = []
    new_u = u.copy()
    u_values.append(new_u)
    h = (tf-t0)/n
    t = t0

    tsol = np.empty(0) # Creates an empty array for t
    tsol = np.append(tsol, t) # Fills in the first element of tsol

    for _ in range(1, n+1):
        new_u = implicit_euler_iter(new_u, t, h, f, newton_iter_num)
        u_values.append(new_u)
        t += h
        tsol = np.append(tsol, t) # Saves it in the tsol array

    return [tsol, u_values]
    
